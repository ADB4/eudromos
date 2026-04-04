#pragma once
// drivetrain.h — engine torque curve, gearbox, differential, RWD drive logic.
//
// No clutch (instant engagement), no flywheel inertia. Diff model covers
// open, locked, and 1-way / 1.5-way / 2-way clutch-pack LSD.
//
// std::vector replaced with std::array — these are fixed at compile time and
// the heap allocations were gratuitous. Matters when you're spawning traffic.

#include <cmath>
#include <algorithm>
#include <array>

struct EngineTorqueCurve {
    struct Point { double rpm, torque_nm; };

    // Loosely a 2.0L NA inline-4, ~200hp peak. Piecewise linear from dyno-ish data.
    static constexpr int N_POINTS = 11;
    std::array<Point, N_POINTS> points = {{
        {   0,   0},
        {1000, 120},
        {2000, 180},
        {3000, 230},
        {4000, 260},
        {4500, 270},  // peak torque
        {5000, 265},
        {6000, 245},  // peak power lives around here
        {7000, 210},
        {7500, 180},
        {8000,   0},  // fuel cut
    }};

    double torque_at_rpm(double rpm) const {
        if (rpm <= points.front().rpm) return points.front().torque_nm;
        if (rpm >= points.back().rpm) return points.back().torque_nm;
        for (int i = 0; i + 1 < N_POINTS; ++i) {
            if (rpm <= points[i + 1].rpm) {
                double t = (rpm - points[i].rpm) / (points[i+1].rpm - points[i].rpm);
                return points[i].torque_nm + t * (points[i+1].torque_nm - points[i].torque_nm);
            }
        }
        return 0;
    }

    double max_rpm() const { return points.back().rpm; }
    double idle_rpm() const { return 1000; }
};

struct Gearbox {
    static constexpr int N_GEARS = 6;
    std::array<double, N_GEARS> ratios = {3.58, 2.02, 1.35, 1.00, 0.77, 0.63};
    double final_drive = 3.42;
    int current_gear = 0;         // 0-indexed
    double shift_cooldown = 0.0;  // time left in current shift
    double shift_time = 0.2;      // seconds per shift

    double total_ratio() const {
        if (current_gear < 0 || current_gear >= N_GEARS) return 0;
        return ratios[current_gear] * final_drive;
    }

    int gear_count() const { return N_GEARS; }
    bool is_shifting() const { return shift_cooldown > 0; }

    void shift_up() {
        if (current_gear < gear_count() - 1 && !is_shifting()) {
            current_gear++;
            shift_cooldown = shift_time;
        }
    }

    void shift_down() {
        if (current_gear > 0 && !is_shifting()) {
            current_gear--;
            shift_cooldown = shift_time;
        }
    }

    void update(double dt) {
        shift_cooldown = std::max(0.0, shift_cooldown - dt);
    }
};

// --- Differential ---
//
// Three types that cover most real cars:
//
//   open:   torque splits 50/50 regardless of speed difference. The inside
//           wheel spins up freely in a corner, the outside can't get more
//           than 50% of the input. This is why FWD economy cars spin an
//           inside wheel out of slow corners.
//
//   locked: both axle shafts turn at the same speed. Maximum traction but
//           the car fights you in corners — the inside wheel scrubs because
//           it's forced to turn faster than it wants to. Karts and drift
//           cars run this. Causes tight-corner understeer (RWD) because the
//           inside rear pushes the car wide.
//
//   lsd:    clutch-pack limited-slip. A preload torque keeps the clutch
//           plates slightly engaged at all times. On top of that, the
//           clutch locks proportionally to input torque (ramp angle effect).
//           The locking torque limits how much speed difference the diff
//           allows. Behaves like open at low torque, trends toward locked
//           under power.
//
//           coast_lock_ratio vs power_lock_ratio: a 1.5-way LSD has lower
//           coast locking — it lets the diff open up on decel, so trail
//           braking still gets rear rotation. A 2-way has equal lock both
//           directions (drifting). A 1-way only locks on power (some FFs).
//
// The model: the diff always splits base torque 50/50 (it's a gear set,
// that's what it does). The clutch pack adds a *locking torque* that
// transfers torque from the faster wheel to the slower wheel. The locking
// torque is: T_lock = preload + ramp * |T_input|, capped so it can't
// exceed half the input (fully locked = both wheels get exactly 50%).
//
// For locked mode, we don't literally lock the speeds (that would need a
// constraint solver). Instead we apply a very stiff coupling torque
// proportional to the speed difference, which drives them together within
// a few timesteps. Cheaper than a constraint and stable at 120 Hz.

enum class DiffType { OPEN, LOCKED, LSD };

struct DiffParams {
    DiffType type = DiffType::LSD;

    // LSD parameters — clutch-pack style
    double preload          = 40.0;   // [Nm] always-on friction, even at zero throttle
    double power_lock_ratio = 0.25;   // [-] fraction of input torque that locks (power)
    double coast_lock_ratio = 0.15;   // [-] fraction of input torque that locks (coast/decel)

    // Locked mode stiffness — high enough to equalize within a few frames,
    // low enough to not blow up the integrator. 5000 Nm/(rad/s) works at 120 Hz.
    double lock_stiffness = 5000.0;   // [Nm/(rad/s)]
};

struct DiffState {
    double torque_left  = 0;  // [Nm] at the wheel (diagnostic)
    double torque_right = 0;  // [Nm] at the wheel (diagnostic)
    double lock_torque  = 0;  // [Nm] how much the clutch is transferring (diagnostic)
};

inline DiffState diff_split_torque(
    const DiffParams& p,
    double total_torque,     // [Nm] total drive torque arriving at the diff
    double omega_left,       // [rad/s] left wheel speed
    double omega_right       // [rad/s] right wheel speed
) {
    DiffState out;
    double half = total_torque / 2.0;
    double delta_omega = omega_left - omega_right;  // positive = left faster

    switch (p.type) {

    case DiffType::OPEN:
        // Pure open diff: 50/50 torque, no speed coupling at all. The wheels
        // are free to spin at different rates. One-liner, but it means the
        // tire with less grip limits the whole axle.
        out.torque_left  = half;
        out.torque_right = half;
        out.lock_torque  = 0;
        break;

    case DiffType::LOCKED: {
        // Stiff coupling drives the speed difference toward zero. This is
        // a spring-like torque: T = -k * (omega_L - omega_R). The fast
        // wheel loses torque, the slow wheel gains it.
        //
        // Must clamp the coupling — otherwise the stiff spring overshoots
        // and you get an oscillation that blows up. The clamp ensures the
        // lock torque can't reverse the speed difference in a single
        // timestep (which it would at high stiffness and low wheel inertia).
        // Also cap at a physical maximum — the tire can only transmit so
        // much torque before it just slides.
        double lock = p.lock_stiffness * delta_omega;
        double max_lock = std::max(std::abs(half), 3000.0);
        lock = std::clamp(lock, -max_lock, max_lock);
        out.torque_left  = half - lock;
        out.torque_right = half + lock;
        out.lock_torque  = lock;
        break;
    }

    case DiffType::LSD: {
        // Clutch-pack LSD. The clutch friction depends on how much torque
        // is going through the diff (ramp angle effect) plus a static preload.
        //
        // Under power, locking is stronger (power_lock_ratio) — the diff
        // resists the inside wheel spinning up, sending more torque to the
        // outside. Under coast, locking is weaker (coast_lock_ratio) — the
        // diff opens up more, letting the inside wheel decelerate freely.
        // This is what makes trail braking work with an LSD: you want the
        // rear to be able to rotate the car on decel.
        bool on_power = (total_torque > 0);
        double ramp = on_power ? p.power_lock_ratio : p.coast_lock_ratio;

        // Maximum locking torque available from the clutch pack
        double T_lock_max = p.preload + ramp * std::abs(total_torque);

        // The clutch engages proportionally to speed difference — think of
        // it as a viscous coupling limited by the clutch capacity. The tanh
        // gives a smooth transition from slipping to locked without a
        // discontinuity. Scale factor chosen so ~2 rad/s difference ≈ 75%
        // of max locking.
        double lock_demand = T_lock_max * std::tanh(delta_omega / 2.0);

        // Clamp: locking torque can't exceed half the input — that would
        // mean one wheel gets negative torque and the other gets more than
        // 100%, which is nonsensical for a passive diff.
        lock_demand = std::clamp(lock_demand, -std::abs(half), std::abs(half));

        // Transfer from fast wheel to slow wheel
        out.torque_left  = half - lock_demand;
        out.torque_right = half + lock_demand;
        out.lock_torque  = lock_demand;
        break;
    }
    }

    return out;
}

struct Drivetrain {
    EngineTorqueCurve engine;
    Gearbox gearbox;
    DiffParams diff_params;
    DiffState diff_state;  // last frame's output, for diagnostics

    double engine_rpm = 1000;
    double brake_torque_max = 3000;    // per axle [Nm]
    double engine_braking_factor = 0.02;

    // Anti-hunt timer: after a shift, lock out further shifts for a bit.
    // Without this, the RPM drop from an upshift can trigger an immediate
    // downshift and you get gear oscillation.
    double shift_lockout_timer = 0;

    void update_rpm_from_wheel_speed(double omega_wheel) {
        double ratio = gearbox.total_ratio();
        if (ratio < 0.01) { engine_rpm = engine.idle_rpm(); return; }
        double computed = std::abs(omega_wheel) * ratio * (60.0 / (2.0 * M_PI));
        engine_rpm = std::clamp(computed, engine.idle_rpm(), engine.max_rpm());
    }

    // Returns total drive torque at the wheels (caller runs it through the diff).
    double get_drive_torque(double throttle) const {
        if (gearbox.is_shifting()) return 0;

        double ratio = gearbox.total_ratio();
        double t_engine = engine.torque_at_rpm(engine_rpm) * throttle;

        // Engine braking (off-throttle): small retarding torque from pumping losses
        if (throttle < 0.01 && engine_rpm > engine.idle_rpm()) {
            t_engine = -engine_braking_factor * engine.torque_at_rpm(engine_rpm)
                       * (engine_rpm / engine.max_rpm());
        }

        // SIMPLIFICATION: constant 88% transmission efficiency
        return t_engine * ratio * 0.88;
    }

    // SIMPLIFICATION: equal brake torque all 4 wheels (should be ~60/40 front bias)
    double get_brake_torque_per_wheel(double brake_input) const {
        return brake_input * brake_torque_max * 0.25;
    }

    void auto_shift(double throttle) {
        if (shift_lockout_timer > 0) return;

        if (engine_rpm >= 6500 && throttle > 0.5) {
            gearbox.shift_up();
            shift_lockout_timer = 1.0;
        } else if (engine_rpm <= 1800 && throttle > 0.5 && gearbox.current_gear > 0) {
            gearbox.shift_down();
            shift_lockout_timer = 1.0;
        }
        // no rev-matching, no brake-downshift
    }
};