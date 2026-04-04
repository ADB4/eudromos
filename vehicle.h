#pragma once
// vehicle.h — rigid body + 4 wheels + weight transfer.
//
// Coordinate frame: X forward, Y left, Z up (SAE-ish but Z-up).
// Position and velocity live in world frame. Tire forces are computed in
// body/tire-local frames, accumulated in body frame, then rotated to world
// for integration. This avoids having to deal with Coriolis/centripetal
// coupling terms in the integrator — learned that the hard way.
//
// Beckman references: Part 1 (weight transfer), Part 2 (slip angle),
// Part 3 (cornering force), Part 5 (centripetal accel), Parts 6-7 (Pacejka).

#include "math_types.h"
#include "tire.h"
#include "drivetrain.h"
#include <array>
#include <cmath>

enum WheelID { FL = 0, FR = 1, RL = 2, RR = 3 };

struct VehicleParams {
    double mass = 1400;         // [kg], ~3100 lbs
    double wheelbase = 2.6;     // [m]
    double track_width = 1.55;  // [m]

    // CG sits 1.17m behind the front axle → ~55/45 front weight bias (FR car)
    double cg_to_front = 1.17;
    double cg_to_rear = 1.43;   // cg_to_front + cg_to_rear = wheelbase
    double cg_height = 0.50;    // [m] above ground. SIMPLIFICATION: no ride height change.

    // Only Izz matters on flat ground. SIMPLIFICATION: no roll/pitch DOF.
    double Izz = 2800;  // yaw moment of inertia [kg·m²]

    // Aero drag. SIMPLIFICATION: no downforce.
    // At ~56 m/s (200 kph): F_drag ≈ 0.5 * 1.225 * 0.35 * 2.2 * 56² ≈ 1500 N
    double drag_coefficient = 0.35;
    double frontal_area = 2.2;
    double air_density = 1.225;

    // SIMPLIFICATION: no suspension (instantaneous weight transfer), same
    // tires all around.
    TireParams tire_params;
};

struct VehicleState {
    Vec3 position = {0, 0, 0};
    Vec3 velocity = {0, 0, 0};
    double yaw = 0;
    double yaw_rate = 0;
    double steer_angle = 0;    // front wheels [rad]
    std::array<TireState, 4> tires;
    Drivetrain drivetrain;

    double speed_mps() const { return velocity.length(); }
    double speed_kph() const { return speed_mps() * 3.6; }
    double lateral_accel_g = 0;

    // Body-frame velocity, recomputed each tick
    double Vx_body = 0;
    double Vy_body = 0;
};

struct Vehicle {
    VehicleParams params;
    VehicleState state;

    // Weight transfer per Beckman Part 1.
    //
    // Static distribution:
    //   F_front = W * (cg_to_rear / wheelbase)     — further CG from axle = more load
    //   F_rear  = W * (cg_to_front / wheelbase)
    //
    // Dynamic transfer:
    //   longitudinal: dFz = m * ax * h / L          — accel loads rear, braking loads front
    //   lateral:      dFz = m * ay * h / track       — cornering loads outside wheels
    //
    // SIMPLIFICATION: lateral transfer is 50/50 front/rear (no roll stiffness model).
    void compute_normal_loads(double ax_body, double ay_body) {
        constexpr double g = 9.81;
        double W = params.mass * g;
        double L = params.wheelbase;
        double tw = params.track_width;
        double h = params.cg_height;

        double F_front = W * (params.cg_to_rear / L);
        double F_rear  = W * (params.cg_to_front / L);
        double dFz_long = params.mass * ax_body * h / L;
        double dFz_lat  = params.mass * ay_body * h / tw;

        state.tires[FL].normal_load = std::max(0.0, (F_front - dFz_long) / 2 - dFz_lat / 2);
        state.tires[FR].normal_load = std::max(0.0, (F_front - dFz_long) / 2 + dFz_lat / 2);
        state.tires[RL].normal_load = std::max(0.0, (F_rear  + dFz_long) / 2 - dFz_lat / 2);
        state.tires[RR].normal_load = std::max(0.0, (F_rear  + dFz_long) / 2 + dFz_lat / 2);
    }

    // Contact patch velocity for each wheel, in the tire's local frame.
    //
    // Each patch velocity = body CG velocity + yaw_rate × offset, then for
    // the front wheels we rotate into the steered tire frame.
    //
    // This is where individual tire slip angles come from — the yaw rate
    // contribution means each wheel sees a different velocity even though
    // they're all attached to the same body. (Beckman Part 2)
    struct WheelVel { double Vx, Vy; };

    std::array<WheelVel, 4> compute_wheel_velocities() const {
        std::array<WheelVel, 4> wv;
        double ht = params.track_width / 2;

        // Wheel positions relative to CG (x fwd, y left)
        struct Pos { double x, y; };
        std::array<Pos, 4> wp = {{
            { params.cg_to_front,  ht},  // FL
            { params.cg_to_front, -ht},  // FR
            {-params.cg_to_rear,   ht},  // RL
            {-params.cg_to_rear,  -ht},  // RR
        }};

        for (int i = 0; i < 4; ++i) {
            double vxp = state.Vx_body - state.yaw_rate * wp[i].y;
            double vyp = state.Vy_body + state.yaw_rate * wp[i].x;

            if (i == FL || i == FR) {
                // rotate into steered frame
                double cd = std::cos(state.steer_angle);
                double sd = std::sin(state.steer_angle);
                wv[i].Vx =  vxp * cd + vyp * sd;
                wv[i].Vy = -vxp * sd + vyp * cd;
            } else {
                wv[i].Vx = vxp;
                wv[i].Vy = vyp;
            }
        }
        return wv;
    }

    // Main simulation step. Semi-implicit Euler at 120 Hz.
    // (Symplectic Euler would be more correct terminology. First-order but
    // energy-conserving, good enough on flat ground without stiff springs.)
    void step(double dt) {
        constexpr double g = 9.81;

        // World → body frame velocity
        Mat3 R = Mat3::from_yaw(state.yaw);
        Mat3 Rt = R.transposed();
        Vec3 vb = Rt * state.velocity;
        state.Vx_body = vb.x;
        state.Vy_body = vb.y;

        // Normal loads. For lateral transfer, use centripetal accel (V * yaw_rate)
        // which is more stable than re-using previous-frame Fy sums.
        double ay_centripetal = state.Vx_body * state.yaw_rate;
        double ax_est = 0;
        for (int i = 0; i < 4; ++i) ax_est += state.tires[i].Fx;
        ax_est /= params.mass;
        compute_normal_loads(ax_est, ay_centripetal);

        // Drivetrain
        double omega_rear = (state.tires[RL].omega + state.tires[RR].omega) / 2;
        state.drivetrain.update_rpm_from_wheel_speed(omega_rear);
        state.drivetrain.auto_shift();
        state.drivetrain.gearbox.update(dt);
        if (state.drivetrain.shift_lockout_timer > 0)
            state.drivetrain.shift_lockout_timer -= dt;

        double drive_torque = state.drivetrain.get_drive_torque();
        double brake_per_wheel = state.drivetrain.get_brake_torque_per_wheel();

        // Tire forces
        auto wvel = compute_wheel_velocities();
        for (int i = 0; i < 4; ++i)
            compute_tire_forces(state.tires[i], params.tire_params, wvel[i].Vx, wvel[i].Vy, dt);

        // Wheel spin dynamics: I_w * domega/dt = drive - brake - Fx*R
        // SIMPLIFICATION: lumped inertia, no separate brake rotor/hub
        constexpr double I_wheel = 1.5;  // [kg·m²]
        double tire_R = params.tire_params.radius;

        for (int i = 0; i < 4; ++i) {
            double torque = 0;

            if (i == RL || i == RR)
                torque += drive_torque / 2;  // RWD, 50/50 split, no diff

            if (state.tires[i].omega > 0.01)       torque -= brake_per_wheel;
            else if (state.tires[i].omega < -0.01)  torque += brake_per_wheel;

            torque -= state.tires[i].Fx * tire_R;  // ground reaction

            state.tires[i].omega += (torque / I_wheel) * dt;

            // Crude wheel-lock prevention (not really ABS, just a clamp)
            if (state.drivetrain.brake_input > 0 && state.tires[i].omega < 0 && state.Vx_body > 0.5)
                state.tires[i].omega = 0;
        }

        // Accumulate forces in body frame, compute yaw torque.
        // Front tire forces get rotated back out of the steered frame.
        double ht = params.track_width / 2;
        Vec3 F_body = {0, 0, 0};
        double yaw_torque = 0;

        for (int i = 0; i < 4; ++i) {
            double fx, fy;

            if (i == FL || i == FR) {
                double cd = std::cos(state.steer_angle);
                double sd = std::sin(state.steer_angle);
                fx = state.tires[i].Fx * cd - state.tires[i].Fy * sd;
                fy = state.tires[i].Fx * sd + state.tires[i].Fy * cd;
            } else {
                fx = state.tires[i].Fx;
                fy = state.tires[i].Fy;
            }

            F_body.x += fx;
            F_body.y += fy;

            // yaw moment = r × F (z component of 2D cross product)
            double rx = (i <= FR) ? params.cg_to_front : -params.cg_to_rear;
            double ry = (i == FL || i == RL) ? ht : -ht;
            yaw_torque += rx * fy - ry * fx;
        }

        // Aero drag in body frame
        double spd = state.speed_mps();
        if (spd > 0.1) {
            double Fd = 0.5 * params.air_density * params.drag_coefficient
                      * params.frontal_area * spd * spd;
            F_body.x -= Fd * (state.Vx_body / spd);
            F_body.y -= Fd * (state.Vy_body / spd);
        }

        // Integrate in world frame. Body-frame force → world, then F=ma.
        Vec3 F_world = R * F_body;
        Vec3 accel = F_world / params.mass;

        state.lateral_accel_g = (F_body.y / params.mass) / g;

        state.velocity += accel * dt;
        state.position += state.velocity * dt;

        state.yaw_rate += (yaw_torque / params.Izz) * dt;
        state.yaw += state.yaw_rate * dt;

        // keep yaw in [-pi, pi]
        while (state.yaw > M_PI)  state.yaw -= 2 * M_PI;
        while (state.yaw < -M_PI) state.yaw += 2 * M_PI;
    }
};
