
#include <glm/glm.hpp>
using namespace glm;

struct LTC {
  // lobe amplitude
  float amplitude;

  // parametric representation
  float m11, m22, m13, m23;
  vec3  X, Y, Z;

  // matrix representation
  mat3  M;
  mat3  invM;
  float detM;

  LTC() {
    amplitude = 1;
    m11       = 1;
    m22       = 1;
    m13       = 0;
    m23       = 0;
    X         = vec3(1, 0, 0);
    Y         = vec3(0, 1, 0);
    Z         = vec3(0, 0, 1);
    update();
  }
/*
  void copy(const LTC& ltc) {
    this->amplitude = ltc.amplitude;
    this->m11       = ltc.m11;
    this->m22       = ltc.m22;
    this->m13       = ltc.m13;
    this->m23       = ltc.m23;
    this->X         = ltc.X;
    this->Y         = ltc.Y;
    this->Z         = ltc.Z;
    this->M         = ltc.M;
    this->invM      = ltc.invM;
    this->detM      = ltc.detM;
  }
*/
  inline void update()  // compute matrix from parameters
  {
    M    = mat3(X, Y, Z) * mat3(m11, 0, 0, 0, m22, 0, m13, m23, 1);
    invM = inverse(M);
    detM = abs(glm::determinant(M));
  }
/*
  void update2() {
    invM = inverse(M);
    detM = abs(glm::determinant(M));
  }
*/


  inline float eval(const vec3& L) const {
    
    vec3 Loriginal = glm::normalize(invM * L);
    vec3 L_        = M * Loriginal;
    
    float l        = glm::length(L_);
    float Jacobian = detM / (l * l * l);

    float D = 1.0f / 3.14159f * glm::max<float>(0.0f, Loriginal.z);

    float res = amplitude * D / Jacobian;
/*
    printf("Loriginal: %f,%f,%f\n", Loriginal.x,Loriginal.y,Loriginal.z);
    printf("L_: %f,%f,%f\n", L_.x,L_.y,L_.z);
    printf("l: %f\n", l);
    printf("Jacobian: %f\n", Jacobian);
    printf("D: %f\n", D);
    printf("amplitude: %f\n", amplitude);
    printf("res: %f\n***\n", res);
   */
    return res;
  }

  inline vec3 sample(const float U1, const float U2) const {
    const float theta = acosf(sqrtf(U1));
    const float phi   = 2.0f * 3.14159f * U2;
    const vec3  L     = normalize(M * vec3(sinf(theta) * cosf(phi),
                                          sinf(theta) * sinf(phi), cosf(theta)));
    return L;
  }
  
    inline float eval_pre(const yocto::vec3f& i){
    vec3 i_ = vec3(i.x,i.y,i.z);
    return eval(i_);
  }

 inline float computeMaxValue() {
  float max_value = 0.0;

  const int Nsample = 64;
  for (int j = 0; j < Nsample; ++j)
    for (int i = 0; i < Nsample; ++i) {
      const float U1 = (i + 0.5f) / (float)Nsample;
      const float U2 = (j + 0.5f) / (float)Nsample;
      max_value = std::max<float>(max_value, eval(sample(U1, U2)));
    }

  return max_value;
}
};


/*
float interpolate_4d(const float lut[], const vec4& uvwx) {
  const int LUT_SIZE = 64;

  // get coordinates normalized
  float s = uvwx.x * (LUT_SIZE - 1);
  float t = uvwx.y * (LUT_SIZE - 1);
  float r = uvwx.z * (LUT_SIZE - 1);
  float q = uvwx.w * (LUT_SIZE - 1);

  // get image coordinates and residuals
  int   i = (int)s, j = (int)t, k = (int)r, l = (int)q;
  int   ii = glm::min(i + 1, LUT_SIZE - 1);
  int   jj = glm::min(j + 1, LUT_SIZE - 1);
  int   kk = glm::min(k + 1, LUT_SIZE - 1);
  int   ll = glm::min(l + 1, LUT_SIZE - 1);
  float u = s - i, v = t - j, w = r - k, x = q - l;

  // trilinear interpolation
  int size2 = LUT_SIZE * LUT_SIZE;
  int size3 = size2 * LUT_SIZE;

  return lut[l * size3 + k * size2 + j * LUT_SIZE + i] * (1 - u) * (1 - v) *
             (1 - w) * (1 - x) +
         lut[l * size3 + k * size2 + j * LUT_SIZE + ii] * u * (1 - v) *
             (1 - w) * (1 - x) +
         lut[l * size3 + k * size2 + jj * LUT_SIZE + i] * (1 - u) * v *
             (1 - w) * (1 - x) +
         lut[l * size3 + kk * size2 + j * LUT_SIZE + i] * (1 - u) * (1 - v) *
             w * (1 - x) +
         lut[l * size3 + kk * size2 + jj * LUT_SIZE + i] * (1 - u) * v * w *
             (1 - x) +
         lut[l * size3 + kk * size2 + j * LUT_SIZE + ii] * u * (1 - v) * w *
             (1 - x) +
         lut[l * size3 + k * size2 + jj * LUT_SIZE + ii] * u * v * (1 - w) *
             (1 - x) +
         lut[l * size3 + kk * size2 + jj * LUT_SIZE + ii] * u * v * w *
             (1 - x) +
         lut[ll * size3 + k * size2 + j * LUT_SIZE + i] * (1 - u) * (1 - v) *
             (1 - w) * x +
         lut[ll * size3 + k * size2 + j * LUT_SIZE + ii] * u * (1 - v) *
             (1 - w) * x +
         lut[ll * size3 + k * size2 + jj * LUT_SIZE + i] * (1 - u) * v *
             (1 - w) * x +
         lut[ll * size3 + kk * size2 + j * LUT_SIZE + i] * (1 - u) * (1 - v) *
             w * x +
         lut[ll * size3 + kk * size2 + jj * LUT_SIZE + i] * (1 - u) * v * w *
             x +
         lut[ll * size3 + kk * size2 + j * LUT_SIZE + ii] * u * (1 - v) * w *
             x +
         lut[ll * size3 + k * size2 + jj * LUT_SIZE + ii] * u * v * (1 - w) *
             x +
         lut[ll * size3 + kk * size2 + jj * LUT_SIZE + ii] * u * v * w * x;
}

void create_tab_on_the_fly(float lut_fly[]){

	//Spherical Coordinates

	//This is fixed due to isotropic assumption
	float phiOutgoing = 0;
	float sinPhiOutgoing = 0;
	float cosPhiOutgoing = 1;

	const int grand = 64;

	float phiRange[grand];
	float thetaRange[grand];
	float roughnessRange[grand];

	float min_roughness = 0.0001f;

	for (int i = 0; i < grand; i++) {
		phiRange[i] = (2 * yocto::pif) * i / (grand - 1);
		thetaRange[i] = (yocto::pif / 2) * i / (grand - 1);
		if (i == 0) {  //avoid singularities
			roughnessRange[i] = min_roughness;
		}
		else {
			roughnessRange[i] = i / float(grand - 1);
			roughnessRange[i] = roughnessRange[i] * roughnessRange[i];
		}
	}
	int count = 0;
	for (float roughness : roughnessRange) {
		for (float phiIncoming : phiRange) {
			for (float thetaIncoming : thetaRange) {
				for (float thetaOutgoing : thetaRange) {


					float sinPhiIncoming = std::sin(phiIncoming);
					float cosPhiIncoming = std::cos(phiIncoming);

					float cosThetaIncoming = std::cos(thetaIncoming);
					float sinThetaIncoming = std::sin(thetaIncoming);

					float cosThetaOutgoing = std::cos(thetaOutgoing);
					float sinThetaOutgoing = std::sin(thetaOutgoing);


					glm::vec3 incoming = glm::vec3(sinThetaIncoming * cosPhiIncoming, sinThetaIncoming * sinPhiIncoming, cosThetaIncoming);
					glm::vec3 outgoing = glm::vec3(sinThetaOutgoing * cosPhiOutgoing, sinThetaOutgoing * sinPhiOutgoing, cosThetaOutgoing);
					
					/*EXTRACT LTC VALUES */
/*
					mat3 M = M_GGX(thetaOutgoing,std::sqrt(roughness));
					
					LTC ltc;
					ltc.M = M;
					ltc.invM = glm::inverse(M);
					ltc.detM = std::abs(glm::determinant(M));
					ltc.amplitude = amplitude_GGX(thetaOutgoing,std::sqrt(roughness));

					float max_value = computeMaxValue();
					float value = ltc.eval(incoming)/max_value;
					if (isnan(value)){
						lut_fly[count] = 0;
					}
					else{
						lut_fly[count] = value;
					}
					count++;
				}
			}
		}
	}
  printf("Tabella caricata\n");
}
*/