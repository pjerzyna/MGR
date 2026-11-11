// calibration.hpp: IMU calibration constants and gyro transformation
#pragma once

#include "camera.hpp"

namespace evimo2 {

// Stałe kalibracyjne zestawu sensorów Prophesee (z sekwencji EVIMO2 o wyraźnym ruchu)
// Ekstrynsyka IMU -> Kamera: obrót o 180 deg wokół osi Z
constexpr double R_IMU2CAM_ARR[9] = {
    -1.0,  0.0, 0.0,
     0.0, -1.0, 0.0,
     0.0,  0.0, 1.0
};

// Sprzętowy bias żyroskopu [rad/s] w układzie IMU
constexpr double GYRO_BIAS_ARR[3] = {-0.0295, 0.0154, 0.0180};

// Typowe opóźnienie zegara IMU względem GT [s] (mediana z 19 sekwencji)
constexpr double LAG_CONST_S = -0.0157;

inline Mat3 default_R_imu2cam() {
    Mat3 R;
    for (int i = 0; i < 9; ++i) R.m[i] = R_IMU2CAM_ARR[i];
    return R;
}

inline Vec3 default_gyro_bias() {
    return Vec3{GYRO_BIAS_ARR[0], GYRO_BIAS_ARR[1], GYRO_BIAS_ARR[2]};
}

// Prędkość kątowa kamery [rad/s] z surowego żyroskopu: R * (omega_imu - bias)
inline Vec3 calibrated_gyro(const Vec3& raw_gyro,
                            const Mat3& R = default_R_imu2cam(),
                            const Vec3& bias = default_gyro_bias()) {
    return R * (raw_gyro - bias);
}

} // namespace evimo2
