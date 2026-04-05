#pragma once
// tire.h — Pacejka "Magic Formula" tire model with combined slip.
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
//
// Load sensitivity: real tires are sublinear — double the normal load gives
// maybe 1.8x the grip, not 2x. The effective mu decreases with load:
//   mu_eff = mu_peak * (Fz_nom / Fz)^load_sensitivity
//
// Combined slip: the old model computed Fx and Fy independently and then
// clamped them with a friction circle. That's wrong — in reality, lateral
// force degrades as soon as you add longitudinal slip, even well below the
// friction limit. You can feel this: trail-brake into a corner and the car
// understeers more even though the tires aren't sliding yet.
//
// The approach here: evaluate each Pacejka channel at its own slip input to
// get the "pure" force (Fx0, Fy0), then scale each down based on how much
// of the friction budget the other channel is consuming. Cosine/sine
// weighting on a normalized slip vector (Beckman Part 7):
//
//   sigma_x = sr / sr_peak             (normalized longitudinal slip)
//   sigma_y = tan(sa) / tan(sa_peak)   (normalized lateral slip)
//   sigma   = sqrt(sigma_x^2 + sigma_y^2)
//
//   Fx = Fx0 * (sigma_x / sigma)
//   Fy = Fy0 * (sigma_y / sigma)
//
// When sigma_y = 0 (pure braking), Fx = Fx0. When sigma_x = 0 (pure
// cornering), Fy = Fy0. When both are present, each gets reduced. The
// friction circle shape emerges naturally from the weighting.

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
    //   0%->0  5%->74%  10%->96%  15%->100%  20%->100%  30%->99%
    PacejkaCoeffs longitudinal = {10.0, 1.9, 1.0, 0.97};

    // Lateral: peaks around 8 deg slip angle, falls off past that.
    //   2->44%  4->78%  6->95%  8->100%  10->97%  15->82%  20->69%
    PacejkaCoeffs lateral = {7.0, 1.9, 1.0, -0.5};

    double mu_peak = 1.0;  // dry tarmac

    // Load sensitivity.
    double load_sensitivity = 0.1;   // [-]
    double Fz_nominal       = 3433;  // [N] ~1400kg/4, overwritten at init

    // Slip values where each Pacejka channel peaks. Used to normalize the
    // combined slip vector so that longitudinal and lateral contribute
    // proportionally to their respective saturation points.
    double sr_peak = 0.15;   // 15% slip ratio
    double sa_peak = 0.14;   // [rad] ~8 deg

    // Relaxation length
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
    double mu_effective = 1.0;  // after load sensitivity (diagnostic)
    double combined_sigma = 0.0; // combined slip magnitude (diagnostic)
};

inline double load_sensitive_mu(double Fz, double mu_peak, double Fz_nom, double sensitivity) {
    if (sensitivity < 1e-6 || Fz_nom < 1.0) return mu_peak;
    double ratio = std::clamp(Fz / Fz_nom, 0.1, 10.0);
    return mu_peak * std::pow(ratio, -sensitivity);
}

inline double pacejka_formula(double x, const PacejkaCoeffs& c, double Fz, double mu) {
    double D = mu * Fz * c.D;
    double Bx = c.B * x;
    return D * std::sin(c.C * std::atan(Bx - c.E * (Bx - std::atan(Bx))));
}

inline double compute_slip_ratio(double omega, double tire_radius, double Vx) {
    double tire_speed = omega * tire_radius;
    constexpr double V_thresh = 2.0;
    double ref = std::max({std::abs(Vx), std::abs(tire_speed), V_thresh});
    double raw = (tire_speed - Vx) / ref;
    double fade = std::tanh(std::max(std::abs(Vx), std::abs(tire_speed)) / V_thresh);
    return raw * fade;
}

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
    if (Fz < 1.0) {
        st.Fx = st.Fy = st.Fy_filtered = 0;
        st.mu_effective = p.mu_peak;
        st.combined_sigma = 0;
        return;
    }

    double mu = load_sensitive_mu(Fz, p.mu_peak, p.Fz_nominal, p.load_sensitivity);
    st.mu_effective = mu;

    // Pure-slip forces — what each channel would produce in isolation.
    double Fx0 = pacejka_formula(st.slip_ratio, p.longitudinal, Fz, mu);
    double Fy0 = -pacejka_formula(st.slip_angle, p.lateral, Fz, mu);

    // Combined slip weighting. Normalize each slip component by its peak
    // value so they contribute proportionally — 50% of peak longitudinal
    // slip "uses up" the same fraction of the friction budget as 50% of
    // peak lateral slip.
    double sigma_x = std::abs(st.slip_ratio) / p.sr_peak;
    double sigma_y = std::abs(std::tan(st.slip_angle)) / std::tan(p.sa_peak);
    double sigma = std::sqrt(sigma_x * sigma_x + sigma_y * sigma_y);
    st.combined_sigma = sigma;

    if (sigma > 1e-6) {
        // Weight each force by its share of the combined slip vector.
        // Pure cornering (sigma_x=0) -> Fy unchanged. Pure braking
        // (sigma_y=0) -> Fx unchanged. Combined -> both reduced.
        double wx = sigma_x / sigma;
        double wy = sigma_y / sigma;
        st.Fx = Fx0 * wx;
        st.Fy = Fy0 * wy;
    } else {
        st.Fx = Fx0;
        st.Fy = Fy0;
    }

    // Friction circle — still needed as a hard backstop. The combined slip
    // weighting handles the smooth interaction, but at extreme slips (both
    // channels past their peaks) the weighting alone can overshoot.
    double f_max = mu * Fz;
    double f_total = std::sqrt(st.Fx * st.Fx + st.Fy * st.Fy);
    if (f_total > f_max && f_total > 1e-6) {
        double s = f_max / f_total;
        st.Fx *= s;
        st.Fy *= s;
    }

    // Relaxation length filter on combined-slip Fy (not pure Fy0). The
    // carcass doesn't know whether the force was reduced by combined slip —
    // it just sees the target force and lags toward it.
    {
        double vx = std::max(std::abs(Vx_local), 0.5);
        double tau = std::max(p.relaxation_length / vx, dt);
        double a = 1.0 - std::exp(-dt / tau);
        double Fy_target = st.Fy;
        st.Fy_filtered += (Fy_target - st.Fy_filtered) * a;
        st.Fy = st.Fy_filtered;
    }

    // Rolling resistance.
    // SIMPLIFICATION: constant Crr, no temp/pressure/speed dependence.
    st.Fx -= p.rolling_resistance * Fz;
}