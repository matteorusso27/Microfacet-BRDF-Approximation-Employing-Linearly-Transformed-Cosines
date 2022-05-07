#include <algorithm>
#include "brdf.h"

class BrdfGGXMultiscatter : public Brdf
{
private:

    struct RayInfo
	{
		// direction
		vec3 w;
		float theta;
		float cosTheta;
		float sinTheta;
		float tanTheta;
		float alpha;
		float Lambda;

		void updateDirection(const vec3& w, const float alpha_x, const float alpha_y)
		{
			this->w = w;
			theta = acosf(w.z);
			cosTheta = w.z;
			sinTheta = sinf(theta);
			tanTheta = sinTheta / cosTheta;
			const float invSinTheta2 = 1.0f / (1.0f - w.z*w.z);
			const float cosPhi2 = w.x*w.x*invSinTheta2;
			const float sinPhi2 = w.y*w.y*invSinTheta2;
			alpha = sqrtf( cosPhi2*alpha_x*alpha_x + sinPhi2*alpha_y*alpha_y ); 
			// Lambda
			if(w.z > 0.9999f)
				Lambda = 0.0f;
			else if(w.z < -0.9999f)
				Lambda = -1.0f;
			else
			{	
				const float a = 1.0f/tanTheta/alpha;
				Lambda = 0.5f*(-1.0f + ((a>0)?1.0f:-1.0f) * sqrtf(1 + 1/(a*a)));
			}
		}

		// height
		float h;
		float C1;
		float G1;

		void updateHeight(const float& h)
		{
			this->h = h;
			C1 = std::min(1.0f, std::max(0.0f, 0.5f*(h+1.0f)));

			if(this->w.z > 0.9999f)
				G1 = 1.0f;
			else if(this->w.z <= 0.0f)
				G1 = 0.0f;
			else
				G1 = powf(this->C1, this->Lambda);
		}
	};

	static inline float generateRandomNumber()
	{
		const float U = ((float)rand()) / (float)RAND_MAX;
		return U;
	}

	static inline bool IsFiniteNumber(float x)
	{
		return (x <= std::numeric_limits<float>::max() && x >= -std::numeric_limits<float>::max()); 
	}

	inline vec3 safe_sqrt(const vec3 &s) const 
	{
        vec3 value;
        for (int i=0; i<3; i++)
            value[i] = std::sqrt(std::max(0.0f, s[i]));
        return value;
    }

	vec2 sampleP22_11(const float theta_i, const float U, const float U_2, const float alpha_x, const float alpha_y) const
	{
		vec2 slope;

		if(theta_i < 0.0001f)
		{
			const float r = sqrtf(U/(1.0f-U));
			const float phi = 6.28318530718f * U_2;
			slope.x = r * cosf(phi);
			slope.y = r * sinf(phi);
			return slope;
		}

		// constant
		const float sin_theta_i = sinf(theta_i);
		const float cos_theta_i = cosf(theta_i);
		const float tan_theta_i = sin_theta_i/cos_theta_i;

		// slope associated to theta_i
		const float slope_i = cos_theta_i/sin_theta_i;

		// projected area
		const float projectedarea = 0.5f * (cos_theta_i + 1.0f);
		if(projectedarea < 0.0001f || projectedarea!=projectedarea)
			return vec2(0,0);
		// normalization coefficient
		const float c = 1.0f / projectedarea;

		const float A = 2.0f*U/cos_theta_i/c - 1.0f;
		const float B = tan_theta_i;
		const float tmp = 1.0f / (A*A-1.0f);

		const float D = sqrtf(std::max(0.0f, B*B*tmp*tmp - (A*A-B*B)*tmp));
		const float slope_x_1 = B*tmp - D;
		const float slope_x_2 = B*tmp + D;
		slope.x = (A < 0.0f || slope_x_2 > 1.0f/tan_theta_i) ? slope_x_1 : slope_x_2;

		float U2;
		float S;
		if(U_2 > 0.5f)
		{
		S = 1.0f;
		U2 = 2.0f*(U_2-0.5f);
		}
		else
		{
		S = -1.0f;
		U2 = 2.0f*(0.5f-U_2);
		}
		const float z = (U2*(U2*(U2*0.27385f-0.73369f)+0.46341f)) / (U2*(U2*(U2*0.093073f+0.309420f)-1.000000f)+0.597999f);
		slope.y = S * z * sqrtf(1.0f+slope.x*slope.x);

		return slope;
	}

    vec3 fresnelConductorExact(float cosThetaI, const vec3 &eta, const vec3 &k) const
	{
		/* Modified from "Optics" by K.D. Moeller, University Science Books, 1988 */

		float cosThetaI2 = cosThetaI*cosThetaI,
			sinThetaI2 = 1-cosThetaI2,
			sinThetaI4 = sinThetaI2*sinThetaI2;

		vec3 temp1 = eta*eta - k*k - vec3(sinThetaI2),
				a2pb2 = safe_sqrt((temp1*temp1 + k*k*eta*eta*4.0f)),
				a     = safe_sqrt(((a2pb2 + temp1) * 0.5f));

		vec3 term1 = a2pb2 + vec3(cosThetaI2),
				term2 = a*(2.0f*cosThetaI);

		vec3 Rs2 = (term1 - term2) / (term1 + term2);

		vec3 term3 = a2pb2*cosThetaI2 + vec3(sinThetaI4),
				term4 = term2*sinThetaI2;

		vec3 Rp2 = Rs2 * (term3 - term4) / (term3 + term4);

		return 0.5f * (Rp2 + Rs2);
	}

	vec3 samplePhaseFunction_conductor(const vec3& wi, const float alpha_x, const float alpha_y, const vec3& m_eta, const vec3& m_k, vec3& weight) const
	{
		const float U1 = generateRandomNumber();
		const float U2 = generateRandomNumber();

		// sample D_wi
		// stretch to match configuration with alpha=1.0	
		const vec3 wi_11 = normalize(vec3(alpha_x * wi.x, alpha_y * wi.y, wi.z));

		// sample visible slope with alpha=1.0
		vec2 slope_11 = sampleP22_11(acosf(wi_11.z), U1, U2, alpha_x, alpha_y);

		// align with view direction
		const float phi = atan2(wi_11.y, wi_11.x);
		vec2 slope(cosf(phi)*slope_11.x - sinf(phi)*slope_11.y, sinf(phi)*slope_11.x + cos(phi)*slope_11.y);

		// stretch back
		slope.x *= alpha_x;
		slope.y *= alpha_y;

		// compute normal
		vec3 wm;
		// if numerical instability
		if( (slope.x != slope.x) || !IsFiniteNumber(slope.x) ) 
		{
			if(wi.z > 0) wm = vec3(0.0f,0.0f,1.0f);
			else wm = normalize(vec3(wi.x, wi.y, 0.0f));
		}
		else
			wm = normalize(vec3(-slope.x, -slope.y, 1.0f));

		// reflect
		const vec3 wo = -wi + 2.0f * wm * dot(wi, wm);
		weight = fresnelConductorExact(dot(wi, wm), m_eta, m_k);

		return wo;
	}

	float D_ggx(const vec3& wm, const float alpha_x, const float alpha_y) const
	{
		if( wm.z <= 0.0f)
			return 0.0f;

		// slope of wm
		const float slope_x = -wm.x/wm.z;
		const float slope_y = -wm.y/wm.z;

		// P22
		const float tmp = 1.0f + slope_x*slope_x/(alpha_x*alpha_x) + slope_y*slope_y/(alpha_y*alpha_y);
		const float P22 = 1.0f / (M_PI * alpha_x * alpha_y) / (tmp * tmp);

		// value
		const float value = P22 / (wm.z*wm.z*wm.z*wm.z);
		return value;
	}

	vec3 evalPhaseFunction_conductor(const RayInfo& ray, const vec3& wo, const float alpha_x, const float alpha_y, const vec3& m_eta, const vec3& m_k) const
	{
		if(ray.w.z > 0.9999f)
			return vec3(0.0f);

		// half vector 
		const vec3 wh = normalize(-ray.w+wo);
		if(wh.z < 0.0f)
			return vec3(0.0f);
		
		// projected area
		float projectedArea;
		if(ray.w.z < -0.9999f)
			projectedArea = 1.0f;
		else 
			projectedArea = ray.Lambda * ray.w.z;

		// value
		const vec3 value = fresnelConductorExact(dot(-ray.w, wh), m_eta, m_k) * std::max(0.0f, dot(-ray.w, wh)) * D_ggx(wh, alpha_x, alpha_y) / 4.0f / projectedArea / dot(-ray.w, wh);
		return value;
	}

	// MIS weights for bidirectional path tracing on the microsurface
	float MISweight_conductor(const vec3& wi, const vec3& wo, const float alpha_x, const float alpha_y) const
	{
		if(wi.x == -wo.x && wi.y == -wo.y && wi.z == -wo.z)
			return 1.0f;
		const vec3 wh = normalize(wi+wo);
		const float value = D_ggx( (wh.z>0) ? wh : -wh , alpha_x, alpha_y);
		return value;
	}

	inline float invC1(const float U) const
	{
		const float h = std::max(-1.0f, std::min(1.0f, 2.0f*U-1.0f));
		return h;	
	}

	inline float sampleHeight(const RayInfo& ray, const float U) const
	{
		if(ray.w.z > 0.9999f)
			return std::numeric_limits<float>::max();
		if(ray.w.z < -0.9999f)
		{
			const float value = invC1(U*ray.C1);
			return value;
		}
		if(fabsf(ray.w.z) < 0.0001f)
			return ray.h;

		// probability of intersection
		if (U > 1.0f - ray.G1) // leave the microsurface
			return std::numeric_limits<float>::max();

		const float h = invC1( 
				ray.C1 / powf((1.0f-U),1.0f/ray.Lambda)
				);
		return h;
	}

	


public:

	virtual float eval(const vec3& V, const vec3& L, const float alpha, float& pdf) const
	{
		if(V.z <= 0)
		{
			pdf = 0;
			return 0;
		}

		// masking
		const float a_V = 1.0f / alpha / tanf(acosf(V.z));
		const float LambdaV = (V.z<1.0f) ? 0.5f * (-1.0f + sqrtf(1.0f + 1.0f/a_V/a_V)) : 0.0f;
		const float G1 = 1.0f / (1.0f + LambdaV);

		// shadowing
		float G2;
		if(L.z <= 0.0f)
			G2 = 0;
		else
		{
			const float a_L = 1.0f / alpha / tanf(acosf(L.z));
			const float LambdaL = (L.z<1.0f) ? 0.5f * (-1.0f + sqrtf(1.0f + 1.0f/a_L/a_L)) : 0.0f;
			G2 = 1.0f / (1.0f + LambdaV + LambdaL);
		}

		// D
		const vec3 H = normalize(V+L);
		const float slopex = H.x/H.z;
		const float slopey = H.y/H.z;
		float D = 1.0f / (1.0f + (slopex*slopex+slopey*slopey)/alpha/alpha);
		D = D*D;
		D = D / (3.14159f * alpha * alpha * H.z*H.z*H.z*H.z);

		pdf = fabsf(D * H.z / 4.0f / dot(V,H));
		float res = eval_conductor(L, V, alpha, alpha, {0, 0, 0}, {1, 1, 1}, 1024).x;

		return res;
	}

	virtual vec3 sample(const vec3& V, const float alpha, const float U1, const float U2) const
	{
		const float phi = 2.0f*3.14159f * U1;
		const float r = alpha*sqrtf(U2/(1.0f-U2));
		const vec3 N = normalize(vec3(r*cosf(phi), r*sinf(phi), 1.0f));
		const vec3 L = -V + 2.0f * N * dot(N, V);
		return L;
	}
	vec3 eval_conductor(const vec3& wi, const vec3& wo, const float alpha_x, const float alpha_y, const vec3& m_eta, const vec3& m_k, const int scatteringOrderMax) const
	{
		
		if(wi.z <= 0 || wo.z <= 0)
			return vec3(0.0f);
		// init
		RayInfo ray;
		ray.updateDirection(-wi, alpha_x, alpha_y);	
		ray.updateHeight(1.0f);
		vec3 energy(1.0f);

		RayInfo ray_shadowing;
		ray_shadowing.updateDirection(wo, alpha_x, alpha_y);

		// eval single scattering	
		// half-vector
		const vec3 wh = normalize(wi+wo);
		const float D = D_ggx(wh, alpha_x, alpha_y);
		const float G2 = 1.0f / (1.0f + (-ray.Lambda-1.0f) + ray_shadowing.Lambda);
		vec3 singleScattering = fresnelConductorExact(dot(-ray.w, wh), m_eta, m_k)  *  D * G2 / (4.0f * wi.z);
		
		// MIS weight 
		float wi_MISweight;

		// multiple scattering
		vec3 multipleScattering(0.0f);
		
		// random walk
		int current_scatteringOrder = 0;	
		while(current_scatteringOrder < scatteringOrderMax)
		{
			// next height
			float U = generateRandomNumber();
			ray.updateHeight( sampleHeight(ray, U) );		
					
			// leave the microsurface?
			if( ray.h == std::numeric_limits<float>::max() )
				break;
			else
				current_scatteringOrder++;

			// next event estimation 
			if( current_scatteringOrder > 1) // single scattering is already computed
			{
				vec3 phasefunction = evalPhaseFunction_conductor(ray, wo, alpha_x, alpha_y, m_eta, m_k); 
				ray_shadowing.updateHeight(ray.h);
				float shadowing = ray_shadowing.G1;
				vec3 I = energy * phasefunction * shadowing;

				// MIS
				const float MIS = wi_MISweight / ( wi_MISweight + MISweight_conductor(-ray.w, wo, alpha_x, alpha_y) );


				if ( IsFiniteNumber(I[0]) )
					multipleScattering += I * MIS;
			}

			// next direction
			vec3 weight;
			ray.updateDirection(samplePhaseFunction_conductor(-ray.w, alpha_x, alpha_y, m_eta, m_k, weight), alpha_x, alpha_y);
			energy = energy * weight;
			ray.updateHeight(ray.h);

			if(current_scatteringOrder == 1)
				wi_MISweight = MISweight_conductor(wi, ray.w, alpha_x, alpha_y);

			// if NaN (should not happen, just in case)
			if( (ray.h != ray.h) || (ray.w.x != ray.w.x)) 
				return vec3(0.0f);
		}
		// 0.5f = MIS weight of singleScattering
		// multipleScattering already weighted by MIS
		
		return singleScattering + 2.f*multipleScattering;
	}

	vec3 eval_conductor_ss(const vec3& wi, const vec3& wo, const float alpha_x, const float alpha_y, const vec3& m_eta, const vec3& m_k, const int scatteringOrderMax) const
	{
		if(wi.z <= 0 || wo.z <= 0)
			return vec3(0.0f);

		// init
		RayInfo ray;
		ray.updateDirection(-wi, alpha_x, alpha_y);	
		ray.updateHeight(1.0f);
		vec3 energy(1.0f);

		RayInfo ray_shadowing;
		ray_shadowing.updateDirection(wo, alpha_x, alpha_y);

		// eval single scattering	
		// half-vector
		const vec3 wh = normalize(wi+wo);
		const float D = D_ggx(wh, alpha_x, alpha_y);
		const float G2 = 1.0f / (1.0f + (-ray.Lambda-1.0f) + ray_shadowing.Lambda);
		vec3 singleScattering = fresnelConductorExact(dot(-ray.w, wh), m_eta, m_k)  *  D * G2 / (4.0f * wi.z);
		
		return singleScattering;
	}

	vec3 sample_conductor(const vec3& wi, const float alpha_x, const float alpha_y, const vec3& m_eta, const vec3& m_k, const int scatteringOrderMax)
{
		auto energy = vec3(1,1,1);

		// init
		RayInfo ray;
		ray.updateDirection(-wi, alpha_x, alpha_y);
		ray.updateHeight(1.0f);
			
		// random walk
		int current_scatteringOrder = 0;
		while(true)
		{
			// next height
			float U = generateRandomNumber();
			ray.updateHeight( sampleHeight(ray, U) );		

			// leave the microsurface?
			if( ray.h == std::numeric_limits<float>::max() )
				break;
			else
				current_scatteringOrder++;

			// next direction
			vec3 weight;
			ray.updateDirection(samplePhaseFunction_conductor(-ray.w, alpha_x, alpha_y, m_eta, m_k, weight), alpha_x, alpha_y);
			energy = energy * weight;
			ray.updateHeight(ray.h);

			// if NaN (should not happen, just in case)
			if( (ray.h != ray.h) || (ray.w.x != ray.w.x)) 
			{
				energy = vec3(0,0,0);
				return vec3(0,0,1);
			}

			if( current_scatteringOrder > scatteringOrderMax )
			{
				energy = vec3(0,0,0);
				return vec3(0,0,1);
			}
		}

		return ray.w;
	}

	float pdf(vec3& incoming,vec3& outgoing,float roughness) const {
		/*
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return 0.0f;
		*/
		// Calculate the reflection half-vector 
		vec3 wh = normalize(outgoing+incoming);

		RayInfo ray;
		const float alpha_x = roughness;
		const float alpha_y = roughness;
		ray.updateDirection(incoming, alpha_x, alpha_y);

		// single-scattering PDF + diffuse 
		// otherwise too many fireflies due to lack of multiple-scattering PDF
		// (MIS works even if the PDF is wrong and not normalized)
		return D_ggx(wh, alpha_x, alpha_y) / (1.0f + ray.Lambda) / (4.0f * incoming.z) + outgoing.z;
	}
};
