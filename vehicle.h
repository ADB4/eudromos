#pragma once
// vehicle.h — rigid body + 4 wheels + suspension + weight transfer.
//
// Coordinate frame: X forward, Y left, Z up (SAE-ish but Z-up).
// Position and velocity live in world frame. Tire forces are computed in
// body/tire-local frames, accumulated in body frame, then rotated to world
// for integration. This avoids having to deal with Coriolis/centripetal
// coupling terms in the integrator — learned that the hard way.
//
// The physics pipeline is split in two:
//   compute_forces() — yours. Tire model, suspension, drivetrain. Produces
//                      a body-frame force + yaw torque. Pure output, no
//                      position/velocity mutation.
//   integrate()      — trivial. Semi-implicit Euler on velocity/position/yaw.
//                      In a game engine this gets replaced by the engine's
//                      physics server — you just feed it the force/torque
//                      from compute_forces().
//   step()           — calls both. For the standalone sim.
//
// Godot note: coordinate system is Z-up here, Godot is Y-up. Isolate the
// conversion at the boundary (one rotation), don't let it leak into the
// physics. The Vec3/Mat3 types are internal to the force computation;
// the engine only sees the VehicleInput/VehicleForces structs.
//
// Beckman references: Part 1 (weight transfer), Part 2 (slip angle),
// Part 3 (cornering force), Part 5 (centripetal accel), Parts 6-7 (Pacejka).

#include "math_types.h"
#include "tire.h"
#include "drivetrain.h"
#include "suspension.h"
#include <array>
#include <cmath>

enum WheelID { FL = 0, FR = 1, RL = 2, RR = 3 };

// --- Input: what the game feeds us each frame ---

struct VehicleInput {
    double throttle   = 0;    // [0, 1]
    double brake      = 0;    // [0, 1]
    double steer_angle = 0;   // front wheels [rad], positive = left
};

// --- Output: what we hand back to the engine ---

struct VehicleForces {
    // Body-frame net force (X fwd, Y left). The engine rotates this to world
    // frame using the body's current orientation, then applies it.
    double Fx_body = 0;
    double Fy_body = 0;

    // Yaw torque about the body's Z axis [N·m].
    double yaw_torque = 0;

    // Lateral accel in g (diagnostic, used for HUD / telemetry)
    double lateral_accel_g = 0;
};

// --- Params: everything that describes the car but doesn't change at runtime ---

struct VehicleParams {
    double mass = 1400;         // [kg], ~3100 lbs
    double wheelbase = 2.6;     // [m]
    double track_width = 1.55;  // [m]

    // CG sits 1.17m behind the front axle → ~55/45 front weight bias (FR car)
    double cg_to_front = 1.17;
    double cg_to_rear = 1.43;   // cg_to_front + cg_to_rear = wheelbase
    double cg_height = 0.50;    // [m] above ground

    // SIMPLIFICATION: no roll/pitch DOF in integrator — just yaw on flat ground.
    double Izz = 2800;  // yaw moment of inertia [kg·m²]

    // Aero drag. SIMPLIFICATION: no downforce.
    // At ~56 m/s (200 kph): F_drag ≈ 0.5 * 1.225 * 0.35 * 2.2 * 56² ≈ 1500 N
    double drag_coefficient = 0.35;
    double frontal_area = 2.2;
    double air_density = 1.225;

    TireParams tire_params;
    SuspensionParams suspension_params;

    double sprung_mass() const {
        return mass - 4.0 * suspension_params.unsprung_mass;
    }
};

// --- State: everything that evolves between frames ---
//
// position/velocity/yaw live here for the standalone sim. When the engine
// owns the rigid body, you still need the rest (tire states, suspension
// deflections, drivetrain internals, wheel spin). The engine feeds you back
// the current velocity/yaw_rate each frame through the body-frame transform.

struct VehicleState {
    // Rigid body — owned by us in standalone, by the engine in-game
    Vec3 position = {0, 0, 0};
    Vec3 velocity = {0, 0, 0};
    double yaw = 0;
    double yaw_rate = 0;

    // Internal physics state — always ours
    std::array<TireState, 4> tires;
    Drivetrain drivetrain;
    SuspensionState suspension;

    // Body-frame velocity, recomputed each tick from the world-frame velocity
    // and current yaw. In-engine, you'd get this from the physics server.
    double Vx_body = 0;
    double Vy_body = 0;

    double speed_mps() const { return velocity.length(); }
    double speed_kph() const { return speed_mps() * 3.6; }
};

struct Vehicle {
    VehicleParams params;
    VehicleState state;

    // Set static normal loads so the car doesn't start with Fz=0 on frame 1.
    // Also sets Fz_nominal for load sensitivity — each corner's reference load
    // is its own static weight. Front corners are heavier on this car (55/45),
    // so they lose proportionally more mu under lateral transfer.
    void init() {
        constexpr double g = 9.81;
        double L = params.wheelbase;
        double W = params.mass * g;
        double F_front = W * (params.cg_to_rear / L) / 2.0;
        double F_rear  = W * (params.cg_to_front / L) / 2.0;
        std::array<double, 4> loads = {{ F_front, F_front, F_rear, F_rear }};
        for (int i = 0; i < 4; ++i) {
            state.tires[i].normal_load = loads[i];
            state.suspension.corners[i].force = loads[i];
        }

        // Fz_nominal = average static corner load. Using the average rather
        // than per-corner values means the front (heavier) tires start with
        // mu slightly below mu_peak and the rear slightly above, which is a
        // mild built-in understeer bias. You could use per-corner Fz_nom
        // instead if you want load sensitivity to be purely about *transfer*.
        params.tire_params.Fz_nominal = W / 4.0;
    }

    // --- Force computation: the part you own forever ---
    //
    // Reads current state + input, writes tire/suspension/drivetrain internals,
    // returns net body-frame force and yaw torque. Does NOT touch position,
    // velocity, yaw, or yaw_rate.

    struct WheelVel { double Vx, Vy; };

    std::array<WheelVel, 4> compute_wheel_velocities(double steer_angle) const {
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
                double cd = std::cos(steer_angle);
                double sd = std::sin(steer_angle);
                wv[i].Vx =  vxp * cd + vyp * sd;
                wv[i].Vy = -vxp * sd + vyp * cd;
            } else {
                wv[i].Vx = vxp;
                wv[i].Vy = vyp;
            }
        }
        return wv;
    }

    VehicleForces compute_forces(const VehicleInput& input, double dt) {
        constexpr double g = 9.81;

        // World → body frame velocity.
        // In-engine, you'd get Vx_body/Vy_body from the physics server's
        // local velocity instead of doing this rotation yourself.
        Mat3 R = Mat3::from_yaw(state.yaw);
        Mat3 Rt = R.transposed();
        Vec3 vb = Rt * state.velocity;
        state.Vx_body = vb.x;
        state.Vy_body = vb.y;

        // Estimate body-frame accelerations for the suspension.
        // Lateral: centripetal (Vx * yaw_rate) is more stable than summing Fy.
        // Longitudinal: sum previous-frame Fx. One frame stale, but it's fine.
        double ay_centripetal = state.Vx_body * state.yaw_rate;
        double ax_est = 0;
        for (int i = 0; i < 4; ++i) ax_est += state.tires[i].Fx;
        ax_est /= params.mass;

        // Suspension
        update_suspension(
            state.suspension, params.suspension_params,
            params.sprung_mass(), ax_est, ay_centripetal,
            params.cg_height, params.cg_to_front, params.cg_to_rear,
            params.track_width, dt
        );
        for (int i = 0; i < 4; ++i)
            state.tires[i].normal_load = state.suspension.corners[i].force;

        // Drivetrain
        double omega_rear = (state.tires[RL].omega + state.tires[RR].omega) / 2;
        state.drivetrain.update_rpm_from_wheel_speed(omega_rear);
        state.drivetrain.auto_shift(input.throttle);
        state.drivetrain.gearbox.update(dt);
        if (state.drivetrain.shift_lockout_timer > 0)
            state.drivetrain.shift_lockout_timer -= dt;

        double drive_torque = state.drivetrain.get_drive_torque(input.throttle);
        double brake_per_wheel = state.drivetrain.get_brake_torque_per_wheel(input.brake);

        // Tire forces
        auto wvel = compute_wheel_velocities(input.steer_angle);
        for (int i = 0; i < 4; ++i)
            compute_tire_forces(state.tires[i], params.tire_params, wvel[i].Vx, wvel[i].Vy, dt);

        // Differential — splits drive torque between rear wheels based on
        // their speed difference and the diff type (open/locked/LSD).
        DiffState diff = diff_split_torque(
            state.drivetrain.diff_params,
            drive_torque,
            state.tires[RL].omega,
            state.tires[RR].omega
        );
        state.drivetrain.diff_state = diff;

        // Wheel spin dynamics: I_w * domega/dt = drive - brake - Fx*R
        // SIMPLIFICATION: lumped inertia, no separate brake rotor/hub
        constexpr double I_wheel = 1.5;  // [kg·m²]
        double tire_R = params.tire_params.radius;

        for (int i = 0; i < 4; ++i) {
            double torque = 0;

            if (i == RL)      torque += diff.torque_left;
            else if (i == RR) torque += diff.torque_right;

            if (state.tires[i].omega > 0.01)       torque -= brake_per_wheel;
            else if (state.tires[i].omega < -0.01)  torque += brake_per_wheel;

            torque -= state.tires[i].Fx * tire_R;  // ground reaction

            state.tires[i].omega += (torque / I_wheel) * dt;

            // Crude wheel-lock prevention (not really ABS, just a clamp)
            if (input.brake > 0 && state.tires[i].omega < 0 && state.Vx_body > 0.5)
                state.tires[i].omega = 0;
        }

        // Accumulate forces in body frame, compute yaw torque.
        // Front tire forces get rotated back out of the steered frame.
        double ht = params.track_width / 2;
        VehicleForces out;
        double yaw_torque = 0;

        for (int i = 0; i < 4; ++i) {
            double fx, fy;

            if (i == FL || i == FR) {
                double cd = std::cos(input.steer_angle);
                double sd = std::sin(input.steer_angle);
                fx = state.tires[i].Fx * cd - state.tires[i].Fy * sd;
                fy = state.tires[i].Fx * sd + state.tires[i].Fy * cd;
            } else {
                fx = state.tires[i].Fx;
                fy = state.tires[i].Fy;
            }

            out.Fx_body += fx;
            out.Fy_body += fy;

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
            out.Fx_body -= Fd * (state.Vx_body / spd);
            out.Fy_body -= Fd * (state.Vy_body / spd);
        }

        out.yaw_torque = yaw_torque;
        out.lateral_accel_g = (out.Fy_body / params.mass) / g;

        return out;
    }

    // --- Integration: the part the engine replaces ---
    //
    // Semi-implicit Euler. Body-frame force → world, then F=ma.
    // In Godot, you'd call compute_forces() in _physics_process(), then do:
    //   var f_world = transform.basis * Vector3(forces.Fx_body, 0, -forces.Fy_body)
    //   apply_central_force(f_world)
    //   apply_torque(Vector3(0, forces.yaw_torque, 0))
    // (with the coordinate swap — Z-up to Y-up)

    void integrate(const VehicleForces& forces, double dt) {
        Mat3 R = Mat3::from_yaw(state.yaw);
        Vec3 F_body = {forces.Fx_body, forces.Fy_body, 0};
        Vec3 F_world = R * F_body;
        Vec3 accel = F_world / params.mass;

        state.velocity += accel * dt;
        state.position += state.velocity * dt;

        state.yaw_rate += (forces.yaw_torque / params.Izz) * dt;
        state.yaw += state.yaw_rate * dt;

        while (state.yaw > M_PI)  state.yaw -= 2 * M_PI;
        while (state.yaw < -M_PI) state.yaw += 2 * M_PI;
    }

    // --- Standalone convenience: compute + integrate in one call ---

    void step(const VehicleInput& input, double dt) {
        VehicleForces forces = compute_forces(input, dt);
        integrate(forces, dt);
    }
};