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

#include "math_types.h"
#include "tire.h"
#include "drivetrain.h"
#include "suspension.h"
#include <array>
#include <cmath>

enum WheelID { FL = 0, FR = 1, RL = 2, RR = 3 };

struct VehicleInput {
    double throttle   = 0;    // [0, 1]
    double brake      = 0;    // [0, 1]
    double steer_angle = 0;   // front wheels [rad], positive = left
};

struct VehicleForces {
    double Fx_body = 0;
    double Fy_body = 0;
    double yaw_torque = 0;
    double lateral_accel_g = 0;
};

struct VehicleParams {
    double mass = 1400;
    double wheelbase = 2.6;
    double track_width = 1.55;
    double cg_to_front = 1.17;
    double cg_to_rear = 1.43;
    double cg_height = 0.50;
    double ackermann_pct = 0.6;
    double Izz = 2800;
    double drag_coefficient = 0.35;
    double frontal_area = 2.2;
    double air_density = 1.225;
    double downforce_ClA = 1.2;
    double aero_balance_front = 0.40;

    TireParams tire_params;
    SuspensionParams suspension_params;

    double sprung_mass() const {
        return mass - 4.0 * suspension_params.unsprung_mass;
    }
};

struct VehicleState {
    Vec3 position = {0, 0, 0};
    Vec3 velocity = {0, 0, 0};
    double yaw = 0;
    double yaw_rate = 0;

    std::array<TireState, 4> tires;
    Drivetrain drivetrain;
    SuspensionState suspension;

    double Vx_body = 0;
    double Vy_body = 0;

    double speed_mps() const { return velocity.length(); }
    double speed_kph() const { return speed_mps() * 3.6; }
};

struct Vehicle {
    VehicleParams params;
    VehicleState state;

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
        params.tire_params.Fz_nominal = W / 4.0;

        // Start with clutch disengaged at standstill so the launch
        // sequence works properly (auto-clutch ramps it in).
        state.drivetrain.clutch_state.clutch_input = 0.0;
    }

    struct WheelVel { double Vx, Vy; };
    struct SteerPair { double left, right; };

    SteerPair compute_steer_angles(double steer_input) const {
        double correction = steer_input * (params.track_width / 2.0)
                          / params.wheelbase * params.ackermann_pct;
        return { steer_input + correction, steer_input - correction };
    }

    std::array<WheelVel, 4> compute_wheel_velocities(const SteerPair& steer) const {
        std::array<WheelVel, 4> wv;
        double ht = params.track_width / 2;

        struct Pos { double x, y; };
        std::array<Pos, 4> wp = {{
            { params.cg_to_front,  ht},
            { params.cg_to_front, -ht},
            {-params.cg_to_rear,   ht},
            {-params.cg_to_rear,  -ht},
        }};

        for (int i = 0; i < 4; ++i) {
            double vxp = state.Vx_body - state.yaw_rate * wp[i].y;
            double vyp = state.Vy_body + state.yaw_rate * wp[i].x;

            if (i == FL || i == FR) {
                double delta = (i == FL) ? steer.left : steer.right;
                double cd = std::cos(delta);
                double sd = std::sin(delta);
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

        Mat3 R = Mat3::from_yaw(state.yaw);
        Mat3 Rt = R.transposed();
        Vec3 vb = Rt * state.velocity;
        state.Vx_body = vb.x;
        state.Vy_body = vb.y;

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

        // Aero downforce
        double V2 = state.Vx_body * state.Vx_body + state.Vy_body * state.Vy_body;
        double F_down = 0.5 * params.air_density * params.downforce_ClA * V2;
        double F_down_front = F_down * params.aero_balance_front;
        double F_down_rear  = F_down * (1.0 - params.aero_balance_front);

        state.tires[FL].normal_load += F_down_front / 2.0;
        state.tires[FR].normal_load += F_down_front / 2.0;
        state.tires[RL].normal_load += F_down_rear  / 2.0;
        state.tires[RR].normal_load += F_down_rear  / 2.0;

        // ── Drivetrain with clutch model ──
        //
        // The engine is always integrated as a free rotational body. The
        // clutch transfers torque between engine and gearbox input shaft
        // based on speed difference and engagement level.
        //
        // Flow:
        //   1. Update shift state machine (clutch engage/disengage ramps)
        //   2. Auto-clutch for launch (ramps engagement at standstill)
        //   3. Auto-shift (triggers new shifts based on RPM)
        //   4. Integrate omega_engine (T_engine - T_clutch) / I_engine
        //   5. Drive torque at wheels = T_clutch × ratio × efficiency
        //   6. Diff splits drive torque between RL and RR
        //   7. Wheel spin: bare I_wheel only (no reflected inertia needed —
        //      the clutch coupling handles everything)

        double omega_rear = (state.tires[RL].omega + state.tires[RR].omega) / 2;

        // Shift state machine: manages clutch ramps during gear changes
        state.drivetrain.update_shift_state(dt);
        state.drivetrain.gearbox.update(dt);
        if (state.drivetrain.shift_lockout_timer > 0)
            state.drivetrain.shift_lockout_timer -= dt;

        // Auto-clutch for launch: ramps engagement when accelerating from stop
        state.drivetrain.update_auto_clutch_launch(input.throttle, omega_rear, dt);

        // Auto-shift: triggers new shift sequences based on RPM
        state.drivetrain.auto_shift(input.throttle);

        // Integrate engine speed. This computes T_clutch internally and
        // stores it in clutch_state.torque_transmitted.
        state.drivetrain.update_engine(input.throttle, omega_rear, dt);

        // Drive torque at wheels comes from the clutch, through the gearbox
        double drive_torque = state.drivetrain.get_drive_torque_from_clutch();

        // Tire forces
        auto steer = compute_steer_angles(input.steer_angle);
        auto wvel = compute_wheel_velocities(steer);
        for (int i = 0; i < 4; ++i)
            compute_tire_forces(state.tires[i], params.tire_params, wvel[i].Vx, wvel[i].Vy, dt);

        // Differential
        DiffState diff = diff_split_torque(
            state.drivetrain.diff_params,
            drive_torque,
            state.tires[RL].omega,
            state.tires[RR].omega
        );
        state.drivetrain.diff_state = diff;

        // Wheel spin dynamics — bare wheel inertia only.
        // No reflected engine inertia needed: the clutch coupling handles
        // the engine-wheel interaction through torque, not through shared
        // inertia. The engine is its own DOF now.
        constexpr double I_wheel = 1.5;  // [kg·m²]
        double tire_R = params.tire_params.radius;

        for (int i = 0; i < 4; ++i) {
            double torque = 0;

            if (i == RL)      torque += diff.torque_left;
            else if (i == RR) torque += diff.torque_right;

            bool is_front = (i == FL || i == FR);
            double brake_torque = state.drivetrain.get_brake_torque_per_wheel(input.brake, is_front);

            if (state.tires[i].omega > 0.01)       torque -= brake_torque;
            else if (state.tires[i].omega < -0.01)  torque += brake_torque;

            torque -= state.tires[i].Fx * tire_R;

            state.tires[i].omega += (torque / I_wheel) * dt;

            if (input.brake > 0 && state.tires[i].omega < 0 && state.Vx_body > 0.5)
                state.tires[i].omega = 0;
        }

        // Force accumulation
        double ht = params.track_width / 2;
        VehicleForces out;
        double yaw_torque = 0;

        for (int i = 0; i < 4; ++i) {
            double fx, fy;

            if (i == FL || i == FR) {
                double delta = (i == FL) ? steer.left : steer.right;
                double cd = std::cos(delta);
                double sd = std::sin(delta);
                fx = state.tires[i].Fx * cd - state.tires[i].Fy * sd;
                fy = state.tires[i].Fx * sd + state.tires[i].Fy * cd;
            } else {
                fx = state.tires[i].Fx;
                fy = state.tires[i].Fy;
            }

            out.Fx_body += fx;
            out.Fy_body += fy;

            double rx = (i <= FR) ? params.cg_to_front : -params.cg_to_rear;
            double ry = (i == FL || i == RL) ? ht : -ht;
            yaw_torque += rx * fy - ry * fx;
        }

        // Aero drag
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

    void step(const VehicleInput& input, double dt) {
        VehicleForces forces = compute_forces(input, dt);
        integrate(forces, dt);
    }
};