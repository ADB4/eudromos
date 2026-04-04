#pragma once
// tire.h — Pacejka "Magic Formula" tire model.
//
// The whole point of this file: given how fast the tire is spinning and how
// fast the contact patch is moving over the ground, figure out how much force
// the tire is generating longitudinally (traction/braking) and laterally
// (cornering). The mapping from slip inputs to forces is the Pacejka formula:
//
//   F = D * sin(C * atan(B*x - E*(B*x - atan(B*x))))
//
// B = stiffness, C = shape, D = peak (mu * Fz), E = curvature.
// Separate coefficients for longitudinal (from slip ratio) and lateral
// (from slip angle). See Beckman Parts 2, 3, 6, 7.

#include "math_types.h"
#include <cmath>
#include <algorithm>

struct PacejkaCoeffs {
    double B;  // stiffness
    double C;  // shape
    double D;  // peak scaling (usually 1.0, multiplied by mu*Fz at runtime)
    double E;  // curvature
};

struct TireParams {
    double radius = 0.33;              // [m]
    double rolling_resistance = 0.015;

    // Longitudinal: peaks around 15% slip ratio.
    //   0%→0  5%→74%  10%→96%  15%→100%  20%→100%  30%→99%
    PacejkaCoeffs longitudinal = {10.0, 1.9, 1.0, 0.97};

    // Lateral: peaks around 8° slip angle, falls off past that.
    //   2°→44%  4°→78%  6°→95%  8°→100%  10°→97%  15%→82%  20°→69%
    PacejkaCoeffs lateral = {7.0, 1.9, 1.0, -0.5};

    double mu_peak = 1.0;  // dry tarmac
    double Fz_nominal = 3434;   // [N] nominal load for load sensitivity (~W/4)

    // Relaxation length — how far the tire rolls before Fy reaches steady
    // state. Without this the yaw response is way too twitchy. ~0.3-0.5m
    // for a passenger car tire. It acts as a first-order lag: tau = L / |Vx|
    double relaxation_length = 0.4;  // [m]
};

struct TireState {
    double normal_load  = 0.0;  // Fz [N]
    double slip_ratio   = 0.0;  // [-1, 1]
    double slip_angle   = 0.0;  // [rad]
    double Fx           = 0.0;  // longitudinal force, tire frame [N]
    double Fy           = 0.0;  // lateral force, tire frame [N]
    double omega        = 0.0;  // wheel angular velocity [rad/s]
    double Fy_filtered  = 0.0;  // after relaxation length lag
};

inline double pacejka_formula(double x, const PacejkaCoeffs& c, double Fz, double mu) {
    double D = mu * Fz * c.D;
    double Bx = c.B * x;
    return D * std::sin(c.C * std::atan(Bx - c.E * (Bx - std::atan(Bx))));
}

// Longitudinal slip (Beckman Part 6). Positive = wheel spinning faster than
// ground (wheelspin). Negative = wheel slower than ground (braking).
// The tanh ramp kills the force at very low speed to stop the tire/accel
// feedback loop from chattering. A brush model or Dahl model would be better.
inline double compute_slip_ratio(double omega, double tire_radius, double Vx) {
    double tire_speed = omega * tire_radius;
    constexpr double V_thresh = 2.0;
    double ref = std::max({std::abs(Vx), std::abs(tire_speed), V_thresh});
    double raw = (tire_speed - Vx) / ref;
    double fade = std::tanh(std::max(std::abs(Vx), std::abs(tire_speed)) / V_thresh);
    return raw * fade;
}

// Slip angle (Beckman Part 2). The angle between where the tire points and
// where it's actually going. Same low-speed fade as slip ratio.
inline double compute_slip_angle(double Vy, double Vx) {
    constexpr double V_thresh = 2.0;
    double abs_vx = std::max(std::abs(Vx), V_thresh);
    double raw = std::atan2(Vy, abs_vx);
    double fade = std::tanh(std::abs(Vx) / V_thresh);
    return raw * fade;
}

inline void compute_tire_forces(
    TireState& st, const TireParams& p,
    double Vx_local, double Vy_local, double dt
) {
    st.slip_ratio = std::clamp(compute_slip_ratio(st.omega, p.radius, Vx_local), -1.0, 1.0);
    st.slip_angle = std::clamp(compute_slip_angle(Vy_local, Vx_local), -M_PI/2, M_PI/2);

    double Fz = st.normal_load;
    if (Fz < 1.0) { st.Fx = st.Fy = st.Fy_filtered = 0; return; }

    st.Fx = pacejka_formula(st.slip_ratio, p.longitudinal, Fz, p.mu_peak);

    // Fy opposes slip angle — this is a restoring force. The negation matters:
    // without it, lateral force reinforces sideslip and the car spins instantly.
    double Fy_raw = -pacejka_formula(st.slip_angle, p.lateral, Fz, p.mu_peak);

    // Relaxation length filter (first-order lag on Fy). The tire carcass has to
    // deform before it generates lateral force. tau = relaxation_length / |Vx|.
    // Without this, steering response is unrealistically instant and you get
    // violent yaw oscillations even at 120 Hz.
    {
        double vx = std::max(std::abs(Vx_local), 0.5);
        double tau = std::max(p.relaxation_length / vx, dt);
        double a = 1.0 - std::exp(-dt / tau);
        st.Fy_filtered += (Fy_raw - st.Fy_filtered) * a;
    }
    st.Fy = st.Fy_filtered;

    // Rolling resistance. Small but it matters for top speed.
    // Applied before the friction circle so it's part of the friction budget.
    // SIMPLIFICATION: constant Crr, no temp/pressure/speed dependence.
    st.Fx -= p.rolling_resistance * Fz;

    // Friction circle (Beckman's "traction circle"). Total force can't exceed
    // mu * Fz. If combined Fx+Fy would bust through, scale both back.
    // SIMPLIFICATION: real tires are slightly elliptical, not circular.
    double f_max = p.mu_peak * Fz;
    double f_total = std::sqrt(st.Fx * st.Fx + st.Fy * st.Fy);
    if (f_total > f_max && f_total > 1e-6) {
        double s = f_max / f_total;
        st.Fx *= s;
        st.Fy *= s;
    }
}