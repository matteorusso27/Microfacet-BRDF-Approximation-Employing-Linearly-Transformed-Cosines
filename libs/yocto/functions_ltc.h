 struct mat33 {
  operator glm::mat3() const {
    return glm::mat3(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]);
  }

  double m[9];
};

#include "tab_ltc.h"
inline mat3 M_GGX(const float theta, const float alpha) {
  int t = std::max<int>(
      0, std::min<int>(dim - 1,
             (int)floorf(theta / (0.5f * 3.14159f) * dim)));
  int a = std::max<int>(
      0, std::min<int>(dim - 1, (int)floorf(sqrtf(alpha) * dim)));

  return tabM[a + t * dim];
}

inline float amplitude_GGX(const float theta, const float alpha) {
  int t = std::max<int>(
      0, std::min<int>(dim - 1,
             (int)floorf(theta / (0.5f * 3.14159f) * dim)));
  int a = std::max<int>(
      0, std::min<int>(dim - 1, (int)floorf(sqrtf(alpha) * dim)));

  return tabAmplitude[a + t * dim];
}

