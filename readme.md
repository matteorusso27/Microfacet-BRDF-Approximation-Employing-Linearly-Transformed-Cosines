The BRDF function describes the reflection of light off a surface, it gives the fraction of reflected light in a given direction. The GGX Microfacet BRDF is currently considered the most realistic parametric BRDF: it gives excellent results for many surfaces and materials, but has a major drawback when materials have a high degree of roughness, making them appear darker, resulting in unnatural products. This is due to a shadowing term that cuts off some light paths at the microfacet level, mean- ing that energy conservation cannot be achieved with this approach

In 2016, Heitz et al. derived a stochastic method to simulate the initial interaction between light rays and the surface and their subsequent bounces. The main principle is to consider the surface as a volume in which one can perform a random walk to evaluate the BRDF more accurately: in this way, conservation of energy is preserved. This work was able to reproduce the albedo and the multiple scattering BRDF with a high degree of accuracy. The main drawback is the integration with modern path tracers, since it is a purely stochastic method. The algorithm can increase rendering time by a factor of 7 to 15, as pointed out by Kulla and Conty.

In ”Real-Time Polygonal-Light Shading with Linearly Transformed Cosines” Heitz et al.'s work, we found their approach could approximate various types of BRDFs using an efficient, fast, and strong approach given by Linearly Transformed Cosines. They approximated the GGX-BRDF as well as other distributions; the results seemed promising and inspired the work proposed in this thesis: the definition of a fit for the Multiple Scattering BRDF using Linearly Transformed Cosines.

The approximation proved convincing at high roughness values, with our model matching the behavior of Heitz et al's. algorithm. We compared our results with scenes rendered using the GGX-BRDF, and we found large improvements due to our energy conserving algorithm

The BRDF function describes the reflection of light off a surface, it gives the fraction of reflected light in a given direction. The GGX Microfacet BRDF is currently considered the most realistic parametric BRDF: it gives excellent results for many surfaces and materials, but has a major drawback when materials have a high degree of roughness, making them appear darker, resulting in unnatural products. This is due to a shadowing term that cuts off some light paths at the microfacet level, mean- ing that energy conservation cannot be achieved with this approach

In 2016, Heitz et al. derived a stochastic method to simulate the initial interaction between light rays and the surface and their subsequent bounces. The main principle is to consider the surface as a volume in which one can perform a random walk to evaluate the BRDF more accurately: in this way, conservation of energy is preserved. This work was able to reproduce the albedo and the multiple scattering BRDF with a high degree of accuracy. The main drawback is the integration with modern path tracers, since it is a purely stochastic method. The algorithm can increase rendering time by a factor of 7 to 15, as pointed out by Kulla and Conty.

In ”Real-Time Polygonal-Light Shading with Linearly Transformed Cosines” Heitz et al.'s work, we found their approach could approximate various types of BRDFs using an efficient, fast, and strong approach given by Linearly Transformed Cosines. They approximated the GGX-BRDF as well as other distributions; the results seemed promising and inspired the work proposed in this thesis: the definition of a fit for the Multiple Scattering BRDF using Linearly Transformed Cosines.

The approximation proved convincing at high roughness values, with our model matching the behavior of Heitz et al's. algorithm. We compared our results with scenes rendered using the GGX-BRDF, and we found large improvements due to our energy conserving algorithm

![alt_text](https://github.com/matteorusso27/Microfacet-BRDF-Approximation-Employing-Linearly-Transformed-Cosines/blob/main/fit_comparison_r01.png)
Multiscattering BRDF (top) and LTC fit (bottom) for roughness: 0.10. The
angle of incidence increases from left to right

![alt_text](https://github.com/matteorusso27/Microfacet-BRDF-Approximation-Employing-Linearly-Transformed-Cosines/blob/main/fit_comparison_r03.png)
Multiscattering BRDF (top) and LTC fit (bottom) for roughness: 0.30. The
angle of incidence increases from left to right

![alt_text](https://github.com/matteorusso27/Microfacet-BRDF-Approximation-Employing-Linearly-Transformed-Cosines/blob/main/fit_comparison_r04.png)
Multiscattering BRDF (top) and LTC fit (bottom) for roughness: 0.40. The
angle of incidence increases from left to right

![alt_text](https://github.com/matteorusso27/Microfacet-BRDF-Approximation-Employing-Linearly-Transformed-Cosines/blob/main/fit_comparison_r05.png)
Multiscattering BRDF (top) and LTC fit (bottom) for roughness: 0.50. The
angle of incidence increases from left to right

![alt_text](https://github.com/matteorusso27/Microfacet-BRDF-Approximation-Employing-Linearly-Transformed-Cosines/blob/main/fit_comparison_r08.png)
Multiscattering BRDF (top) and LTC fit (bottom) for roughness: 0.80. The
angle of incidence increases from left to right

![alt_text](https://github.com/matteorusso27/Microfacet-BRDF-Approximation-Employing-Linearly-Transformed-Cosines/blob/main/fit_comparison_r10.png)
Multiscattering BRDF (top) and LTC fit (bottom) for roughness: 1.0 The
angle of incidence increases from left to right

To evaluate our model, we based our results on the Furnace test, a sample scene in
which we compute the Rendering Equation [ Kaj86] using a constant environment
map. In order to have no energy loss or gain, the objects in the scene should present
the same color as the background. To run our renderings, we used Yocto GL [FP19 ],
a tiny C++ library for physically-based graphics that provided libraries to support
base scenes, math functions, file management, Path Tracing, and more.
As we can see from the furnace test of Fig. 4.1, the top row represents spheres that
were rendered using the GGX Microfacet BRDF [BW07], the second row represents
the renders using the Multiple Scattering BRDF [ EH16a ] and the last row is our
Linearly Transformed Cosines [EH16b ] approximation of the Energy Conserving
BRDF. Roughness increases from left to right. In the second and third row, since
energy is not lost or gained, the spheres blend in with the environment, resulting in
an energy-conserving model even at high roughness values.

![alt_text](https://github.com/matteorusso27/Microfacet-BRDF-Approximation-Employing-Linearly-Transformed-Cosines/blob/main/basic_furnace_comparison.png)
Furnace test for conductors. Top row: no energy compensation. Middle row:
energy compensation with the use of multiple scattering algorithm. Bottom row: energy
compensation with the use of LTC approximation. Roughness increases from left to
right

![alt_text](https://github.com/matteorusso27/Microfacet-BRDF-Approximation-Employing-Linearly-Transformed-Cosines/blob/main/dragon_ggx.png)
![alt_text](https://github.com/matteorusso27/Microfacet-BRDF-Approximation-Employing-Linearly-Transformed-Cosines/blob/main/dragon_ltc.png)
Comparison between no Energy Compensation model(top) and LTC algo-
rithm(bottom)
