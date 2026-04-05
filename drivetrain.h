#pragma once
// drivetrain.h — engine, clutch, gearbox, differential, RWD.
//
// The engine is a separate rotational DOF (omega_engine) coupled to the
// wheels through a friction-disk clutch and gearbox. omega_engine is always
// integrated independently — the clutch transfers torque between the engine
// and the gearbox input shaft, but never locks them rigidly.
//
// Clutch model: uses the same tanh-based smooth coupling as the LSD diff.
// The clutch has a torque capacity (T_clutch_max) set by spring pressure.
// When the speed difference between engine and gearbox input is small, the
// clutch transmits nearly all engine torque. When the difference is large
// (launch, shift re-engagement), it slips — the transmitted torque is
// limited by T_clutch_max. This avoids a discrete lock/slip state machine
// and is unconditionally stable at 120 Hz.
//
// Launch behavior: the clutch engages progressively based on a clutch_input
// ramp (0 = disengaged, 1 = fully clamped). At standstill with throttle
// applied, the engine spins up freely until the clutch begins to bite, then
// torque flows to the wheels proportional to the engagement and speed diff.
//
// Shifts: clutch opens (clutch_input = 0), gear changes, clutch re-engages
// with a fast ramp (~125ms). Ignition is cut during upshifts to prevent the
// engine from free-revving to the limiter.

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
    int current_gear = 0;
    double shift_cooldown = 0.0;
    double shift_time = 0.2;

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

enum class DiffType { OPEN, LOCKED, LSD };

struct DiffParams {
    DiffType type = DiffType::LSD;
    double preload          = 40.0;
    double power_lock_ratio = 0.25;
    double coast_lock_ratio = 0.15;
    double lock_stiffness = 5000.0;
};

struct DiffState {
    double torque_left  = 0;
    double torque_right = 0;
    double lock_torque  = 0;
};

inline DiffState diff_split_torque(
    const DiffParams& p,
    double total_torque,
    double omega_left,
    double omega_right
) {
    DiffState out;
    double half = total_torque / 2.0;
    double delta_omega = omega_left - omega_right;

    switch (p.type) {

    case DiffType::OPEN:
        out.torque_left  = half;
        out.torque_right = half;
        out.lock_torque  = 0;
        break;

    case DiffType::LOCKED: {
        double lock = p.lock_stiffness * delta_omega;
        double max_lock = std::max(std::abs(half), 3000.0);
        lock = std::clamp(lock, -max_lock, max_lock);
        out.torque_left  = half - lock;
        out.torque_right = half + lock;
        out.lock_torque  = lock;
        break;
    }

    case DiffType::LSD: {
        bool on_power = (total_torque > 0);
        double ramp = on_power ? p.power_lock_ratio : p.coast_lock_ratio;
        double T_lock_max = p.preload + ramp * std::abs(total_torque);
        double lock_demand = T_lock_max * std::tanh(delta_omega / 2.0);
        lock_demand = std::clamp(lock_demand, -std::abs(half), std::abs(half));
        out.torque_left  = half - lock_demand;
        out.torque_right = half + lock_demand;
        out.lock_torque  = lock_demand;
        break;
    }
    }

    return out;
}

// --- Clutch ---
//
// Friction-disk clutch between engine and gearbox input shaft. Same tanh
// pattern as the LSD diff: smooth transition, no state machine for lock/slip,
// unconditionally stable at 120 Hz.
//
//   omega_gearbox_input = omega_wheel_avg × gear_ratio
//   delta_omega = omega_engine - omega_gearbox_input
//   T_clutch = clutch_input × T_max × tanh(delta_omega / omega_scale)
//
// The clutch torque acts on both sides: it decelerates the engine and
// accelerates the wheels. The engine sees:
//   I_engine × d(omega_engine)/dt = T_engine - T_clutch
// The wheels see: T_clutch × ratio × efficiency as drive torque.

struct ClutchParams {
    // Maximum torque when fully clamped. 1.3-1.5× peak engine torque.
    // Peak is 270 Nm → 400 Nm gives good margin without excess slip.
    double torque_capacity = 400.0;  // [Nm]

    // Speed difference scale for tanh. At omega_scale, the clutch transmits
    // ~76% of capacity. 5 rad/s ≈ 48 RPM — quite stiff, realistic for a
    // properly-sized friction disk. Most slip happens in the first ~100 RPM.
    double omega_scale = 5.0;  // [rad/s]

    // Auto-clutch engagement rate for launch [1/s]. 2.0 = 0→1 in 0.5s.
    double engage_rate = 2.0;

    // Disengage rate for shifts. Fast — clutch opens in ~50ms.
    double disengage_rate = 20.0;  // [1/s]

    // Re-engage rate after shift completion. ~125ms to full engagement.
    double shift_engage_rate = 8.0;  // [1/s]
};

struct ClutchState {
    double clutch_input = 1.0;        // [0,1] engagement level
    double torque_transmitted = 0.0;  // [Nm] actual torque this frame
    double slip_speed = 0.0;          // [rad/s] omega_engine - omega_gearbox_input
};

// --- Drivetrain ---

struct Drivetrain {
    EngineTorqueCurve engine;
    Gearbox gearbox;
    DiffParams diff_params;
    DiffState diff_state;
    ClutchParams clutch_params;
    ClutchState clutch_state;

    double engine_rpm = 1000;
    double brake_torque_max = 3000;
    double engine_braking_factor = 0.02;
    double brake_bias_front = 0.6;
    double shift_lockout_timer = 0;

    // Engine inertia
    double I_engine = 0.18;           // [kg·m²]
    double omega_engine = 104.72;     // [rad/s] = 1000 rpm
    double transmission_efficiency = 0.88;

    // Shift state machine: NONE → DISENGAGE → SHIFTING → ENGAGE → NONE
    enum class ShiftPhase { NONE, DISENGAGE, SHIFTING, ENGAGE };
    ShiftPhase shift_phase = ShiftPhase::NONE;
    bool ignition_cut = false;
    bool shift_is_upshift = false;

    // ── RPM ──

    void update_rpm_from_omega_engine() {
        engine_rpm = std::clamp(
            std::abs(omega_engine) * (60.0 / (2.0 * M_PI)),
            engine.idle_rpm(),
            engine.max_rpm()
        );
    }

    // ── Engine torque ──
    //
    // Includes ignition cut (upshifts) and idle governor (anti-stall).

    double get_engine_torque(double throttle) const {
        if (ignition_cut) return 0;

        double t_engine = engine.torque_at_rpm(engine_rpm) * throttle;

        // Engine braking
        if (throttle < 0.01 && engine_rpm > engine.idle_rpm()) {
            t_engine = -engine_braking_factor * engine.torque_at_rpm(engine_rpm)
                       * (engine_rpm / engine.max_rpm());
        }

        // Idle governor: prevent stalling. Linear ramp from 0 Nm at 1100 RPM
        // to 30 Nm at idle (1000 RPM). Gentle enough not to fight the clutch
        // during normal operation.
        if (engine_rpm < 1100) {
            double idle_authority = 30.0 * (1.0 - (engine_rpm - engine.idle_rpm()) / 100.0);
            idle_authority = std::clamp(idle_authority, 0.0, 30.0);
            t_engine = std::max(t_engine, idle_authority);
        }

        return t_engine;
    }

    // ── Clutch torque ──

    double compute_clutch_torque(double omega_gearbox_input) {
        double delta_omega = omega_engine - omega_gearbox_input;
        clutch_state.slip_speed = delta_omega;

        double T_max = clutch_state.clutch_input * clutch_params.torque_capacity;
        if (T_max < 0.1) {
            clutch_state.torque_transmitted = 0;
            return 0;
        }

        double T_clutch = T_max * std::tanh(delta_omega / clutch_params.omega_scale);
        clutch_state.torque_transmitted = T_clutch;
        return T_clutch;
    }

    // ── Drive torque at wheels ──

    double get_drive_torque_from_clutch() const {
        double ratio = gearbox.total_ratio();
        return clutch_state.torque_transmitted * ratio * transmission_efficiency;
    }

    // Legacy interface for energy auditor
    double get_drive_torque(double /*throttle*/) const {
        return get_drive_torque_from_clutch();
    }

    // ── Brakes ──

    double get_brake_torque_per_wheel(double brake_input, bool is_front) const {
        double total = brake_input * brake_torque_max;
        double axle_frac = is_front ? brake_bias_front : (1.0 - brake_bias_front);
        return total * axle_frac / 2.0;
    }

    // ── Auto-clutch for launch ──

    void update_auto_clutch_launch(double throttle, double omega_wheel_avg, double dt) {
        if (shift_phase != ShiftPhase::NONE) return;

        bool at_standstill = std::abs(omega_wheel_avg) < 2.0;

        if (at_standstill && throttle > 0.05 && clutch_state.clutch_input < 1.0) {
            double rate = clutch_params.engage_rate * (0.5 + throttle * 0.5);
            clutch_state.clutch_input += rate * dt;
            clutch_state.clutch_input = std::min(clutch_state.clutch_input, 1.0);
        } else if (at_standstill && throttle < 0.01) {
            clutch_state.clutch_input = 1.0;
        }

        if (!at_standstill && shift_phase == ShiftPhase::NONE) {
            clutch_state.clutch_input = 1.0;
        }
    }

    // ── Shift logic ──

    void auto_shift(double throttle) {
        if (shift_lockout_timer > 0) return;
        if (shift_phase != ShiftPhase::NONE) return;

        if (engine_rpm >= 6500 && throttle > 0.5) {
            shift_phase = ShiftPhase::DISENGAGE;
            shift_is_upshift = true;
            shift_lockout_timer = 1.0;
        } else if (engine_rpm <= 1800 && throttle > 0.5 && gearbox.current_gear > 0) {
            shift_phase = ShiftPhase::DISENGAGE;
            shift_is_upshift = false;
            shift_lockout_timer = 1.0;
        }
    }

    void update_shift_state(double dt) {
        switch (shift_phase) {

        case ShiftPhase::NONE:
            ignition_cut = false;
            break;

        case ShiftPhase::DISENGAGE:
            clutch_state.clutch_input -= clutch_params.disengage_rate * dt;
            if (clutch_state.clutch_input <= 0) {
                clutch_state.clutch_input = 0;
                if (shift_is_upshift) {
                    gearbox.shift_up();
                    ignition_cut = true;
                } else {
                    gearbox.shift_down();
                    ignition_cut = false;
                }
                shift_phase = ShiftPhase::SHIFTING;
            }
            break;

        case ShiftPhase::SHIFTING:
            clutch_state.clutch_input = 0;
            if (!gearbox.is_shifting()) {
                shift_phase = ShiftPhase::ENGAGE;
                ignition_cut = false;
            }
            break;

        case ShiftPhase::ENGAGE:
            clutch_state.clutch_input += clutch_params.shift_engage_rate * dt;
            if (clutch_state.clutch_input >= 1.0) {
                clutch_state.clutch_input = 1.0;
                shift_phase = ShiftPhase::NONE;
            }
            break;
        }
    }

    // ── Main update: integrate omega_engine ──
    //
    // The engine is always a free rotational body:
    //   I_engine × d(omega_engine)/dt = T_engine - T_clutch
    //
    // When the clutch is fully engaged and speeds are matched, T_clutch ≈
    // T_engine and omega_engine barely changes — the system behaves like
    // the old reflected-inertia model. When the clutch is slipping or open,
    // the engine evolves independently.

    void update_engine(double throttle, double omega_wheel_avg, double dt) {
        double ratio = gearbox.total_ratio();
        double omega_gearbox_input = (ratio > 0.01)
            ? std::abs(omega_wheel_avg) * ratio
            : 0.0;

        double T_clutch = compute_clutch_torque(omega_gearbox_input);
        double T_engine = get_engine_torque(throttle);

        double alpha_engine = (T_engine - T_clutch) / I_engine;
        omega_engine += alpha_engine * dt;

        double omega_idle = engine.idle_rpm() * (2.0 * M_PI / 60.0);
        double omega_max  = engine.max_rpm()  * (2.0 * M_PI / 60.0);
        omega_engine = std::clamp(omega_engine, omega_idle, omega_max);

        update_rpm_from_omega_engine();
    }
};