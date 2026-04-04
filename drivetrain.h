#pragma once
// drivetrain.h — engine torque curve, gearbox, RWD drive logic.
//
// No clutch (instant engagement), no diff (50/50 split), no flywheel inertia.
// Just enough to get plausible torque to the rear wheels through a 6-speed box.

#include <cmath>
#include <algorithm>
#include <vector>

struct EngineTorqueCurve {
    struct Point { double rpm, torque_nm; };
    std::vector<Point> points;

    // Loosely a 2.0L NA inline-4, ~200hp peak. Piecewise linear from dyno-ish data.
    EngineTorqueCurve() {
        points = {
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
        };
    }

    double torque_at_rpm(double rpm) const {
        if (points.empty()) return 0;
        if (rpm <= points.front().rpm) return points.front().torque_nm;
        if (rpm >= points.back().rpm) return points.back().torque_nm;
        for (size_t i = 0; i + 1 < points.size(); ++i) {
            if (rpm <= points[i + 1].rpm) {
                double t = (rpm - points[i].rpm) / (points[i+1].rpm - points[i].rpm);
                return points[i].torque_nm + t * (points[i+1].torque_nm - points[i].torque_nm);
            }
        }
        return 0;
    }

    double max_rpm() const { return points.empty() ? 0 : points.back().rpm; }
    double idle_rpm() const { return 1000; }
};

struct Gearbox {
    std::vector<double> ratios;
    double final_drive = 3.42;
    int current_gear = 0;         // 0-indexed
    double shift_cooldown = 0.0;  // time left in current shift
    double shift_time = 0.2;      // seconds per shift

    Gearbox() {
        ratios = {3.58, 2.02, 1.35, 1.00, 0.77, 0.63}; // typical sport sedan
    }

    double total_ratio() const {
        if (current_gear < 0 || current_gear >= (int)ratios.size()) return 0;
        return ratios[current_gear] * final_drive;
    }

    int gear_count() const { return (int)ratios.size(); }
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

struct Drivetrain {
    EngineTorqueCurve engine;
    Gearbox gearbox;

    double throttle = 0;
    double brake_input = 0;
    double brake_torque_max = 3000;    // per axle [Nm]
    double engine_rpm = 1000;
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

    // Returns total drive torque at the wheels (caller splits between L/R).
    double get_drive_torque() const {
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
    double get_brake_torque_per_wheel() const {
        return brake_input * brake_torque_max * 0.25;
    }

    void auto_shift() {
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
