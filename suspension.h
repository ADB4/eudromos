#pragma once
// suspension.h — per-corner spring + damper + anti-roll bars.
//
// The old compute_normal_loads() was purely algebraic: Fz = static ± m*a*h/t.
// That's the steady-state answer, but it assumes infinitely stiff springs —
// weight transfer is instant, the car feels telepathic. This adds actual
// spring/damper dynamics so the load takes ~100-150ms to settle, with possible
// overshoot if you're underdamped. Same idea as tire relaxation length: the
// physics needs a first-order lag or the response is unrealistically fast.
//
// I tried a full mass-spring-damper ODE first (second-order on deflection)
// and it blew up. The sprung mass is already being integrated in vehicle.step()
// via F=ma — adding a redundant mass term to the suspension double-counts
// inertia. The fix is to model deflection as a first-order filter toward
// equilibrium, with tau = C_damper / K_spring. Same exponential form as the
// tire relaxation filter. It's stable, it's cheap, and it gives you the
// transient delay without fighting the integrator.
//
// ARBs (anti-roll bars) resist the difference in deflection between left and
// right wheels on the same axle. They add roll stiffness without changing ride
// stiffness — the main knob for tuning under/oversteer balance.

#include <cmath>
#include <algorithm>
#include <array>

struct SuspensionCornerParams {
    double spring_rate      = 25000;   // [N/m], ~143 lb/in
    double damping_bump     = 2500;    // [N·s/m] compression
    double damping_rebound  = 3500;    // [N·s/m] extension — stiffer to control oscillation
    double max_compression  = 0.10;    // [m] bump stop
    double max_extension    = 0.12;    // [m] droop stop
    double bump_stop_rate   = 150000;  // [N/m] progressive rubber, last 20mm
    double bump_stop_engage = 0.08;    // [m] 20mm before max
};

struct SuspensionParams {
    std::array<SuspensionCornerParams, 4> corners;  // all identical by default

    // ARBs. Front/rear ratio is the primary balance knob.
    //   more front ARB → front transfers more → front saturates first → understeer
    //   more rear ARB  → rear transfers more  → rear saturates first → oversteer
    double arb_rate_front = 25000;  // [N/m]
    double arb_rate_rear  = 15000;  // [N/m]

    double unsprung_mass = 40;  // [kg] per corner (wheel + hub + brake + knuckle)
};

struct SuspensionCornerState {
    double deflection = 0.0;    // [m] positive = compression
    double velocity   = 0.0;    // [m/s] d(deflection)/dt
    double force      = 0.0;    // [N] normal load at contact patch
    double spring_force = 0.0;  // [N] diagnostic
    double damper_force = 0.0;  // [N] diagnostic
    double arb_force    = 0.0;  // [N] diagnostic
};

struct SuspensionState {
    std::array<SuspensionCornerState, 4> corners;
    double roll_angle  = 0.0;  // [rad] positive = body leans right
    double pitch_angle = 0.0;  // [rad] positive = nose down
};

inline void update_suspension(
    SuspensionState& state,
    const SuspensionParams& params,
    double sprung_mass,
    double ax_body,
    double ay_body,
    double cg_height,
    double cg_to_front,
    double cg_to_rear,
    double track_width,
    double dt
) {
    constexpr double g = 9.81;
    double wheelbase = cg_to_front + cg_to_rear;
    double half_track = track_width / 2.0;

    // Static load per corner, sprung mass only
    double W_sprung = sprung_mass * g;
    double front_frac = cg_to_rear / wheelbase;
    double rear_frac  = cg_to_front / wheelbase;
    double unsprung_weight = params.unsprung_mass * g;

    std::array<double, 4> static_load = {{
        W_sprung * front_frac / 2.0,
        W_sprung * front_frac / 2.0,
        W_sprung * rear_frac  / 2.0,
        W_sprung * rear_frac  / 2.0,
    }};

    // Roll stiffness per axle determines how lateral load transfer splits
    // front/rear. This is the thing that actually controls steady-state balance.
    //
    // The ARB contributes to the roll stiffness that determines the front/rear
    // split of total lateral load transfer. The deflection target is then
    // derived from the target load (which includes the ARB's effect on the
    // split). The force output uses spring + damper only — no separate ARB
    // force term — because the target deflection already encodes the ARB's
    // contribution. Adding an explicit ARB force would double-count it.
    //
    //   K_roll = (K_spring / 2) * (t/2)^2 + K_arb * (t/2)^2
    //
    // The /2 on the spring: in roll, one side compresses while the other
    // extends. The pair acts at half rate across the full track.
    double ht2 = half_track * half_track;
    double K_roll_front = (params.corners[0].spring_rate / 2.0) * ht2
                        + params.arb_rate_front * ht2;
    double K_roll_rear  = (params.corners[2].spring_rate / 2.0) * ht2
                        + params.arb_rate_rear  * ht2;
    double K_roll_total = K_roll_front + K_roll_rear;

    double front_roll_frac = (K_roll_total > 1.0) ? K_roll_front / K_roll_total : 0.5;

    // Lateral load transfer per axle, split by roll stiffness
    double M_roll = sprung_mass * ay_body * cg_height;
    double dFz_lat_front = M_roll * front_roll_frac / track_width;
    double dFz_lat_rear  = M_roll * (1.0 - front_roll_frac) / track_width;

    // Longitudinal load transfer — same as the old algebraic model
    double dFz_long = sprung_mass * ax_body * cg_height / wheelbase;

    // Target load: the old algebraic answer. The suspension filters toward this.
    std::array<double, 4> target_load = {{
        static_load[0] - dFz_long / 2.0 - dFz_lat_front,
        static_load[1] - dFz_long / 2.0 + dFz_lat_front,
        static_load[2] + dFz_long / 2.0 - dFz_lat_rear,
        static_load[3] + dFz_long / 2.0 + dFz_lat_rear,
    }};

    // Target deflection: how far the spring needs to move from static preload.
    // This already includes the ARB's effect via the K_roll front/rear split —
    // no separate ARB force term in the force output.
    std::array<double, 4> target_defl;
    for (int i = 0; i < 4; ++i)
        target_defl[i] = (target_load[i] - static_load[i]) / params.corners[i].spring_rate;

    // Per-corner: first-order filter on deflection, asymmetric damping.
    //
    // tau = C / K. With these defaults:
    //   bump:    2500 / 25000 = 100ms
    //   rebound: 3500 / 25000 = 140ms
    //
    // Compression is faster than rebound — outside wheel loads up quicker
    // than the inside wheel unloads. Real cars do this too.
    for (int i = 0; i < 4; ++i) {
        auto& cs = state.corners[i];
        const auto& cp = params.corners[i];
        double defl_old = cs.deflection;

        double c_damp = (target_defl[i] >= cs.deflection) ? cp.damping_bump : cp.damping_rebound;
        double tau = std::max(c_damp / cp.spring_rate, dt);

        // Same exponential filter as tire relaxation length
        double alpha = 1.0 - std::exp(-dt / tau);
        cs.deflection += (target_defl[i] - cs.deflection) * alpha;

        cs.deflection = std::clamp(cs.deflection, -cp.max_extension, cp.max_compression);

        // Bump stop: stiff progressive spring near full compression
        double bump_force = 0.0;
        if (cs.deflection > cp.bump_stop_engage)
            bump_force = cp.bump_stop_rate * (cs.deflection - cp.bump_stop_engage);

        cs.velocity = (cs.deflection - defl_old) / dt;

        cs.spring_force = cp.spring_rate * cs.deflection + static_load[i] + bump_force;
        cs.damper_force = c_damp * cs.velocity;
        cs.arb_force = 0;  // ARB effect is in K_roll → target_defl, not a separate force

        cs.force = std::max(0.0, cs.spring_force + cs.damper_force + unsprung_weight);
    }

    // Body attitude from deflection geometry (small angle)
    double left_avg  = (state.corners[0].deflection + state.corners[2].deflection) / 2.0;
    double right_avg = (state.corners[1].deflection + state.corners[3].deflection) / 2.0;
    double front_avg = (state.corners[0].deflection + state.corners[1].deflection) / 2.0;
    double rear_avg  = (state.corners[2].deflection + state.corners[3].deflection) / 2.0;

    state.roll_angle  = (right_avg - left_avg) / track_width;
    state.pitch_angle = (front_avg - rear_avg) / wheelbase;
}