// main.cpp — runs the sim, dumps CSV.
//
// g++ -std=c++17 -O2 -o vehicle_sim main.cpp -lm
// ./vehicle_sim [scenario] > output.csv
//
// Scenarios:
//   mixed       (default) original 25s: launch → brake → corner → trail brake
//   skidpad     constant-radius left turn, slowly increasing speed
//   brake_stop  full throttle to ~100 kph, then brake to standstill
//   step_steer  cruise at ~80 kph, then step steer input

#include "vehicle.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstring>

constexpr double DT       = 1.0 / 120.0;  // 120 Hz
constexpr int    LOG_SKIP = 12;            // log every 12th tick → ~10 Hz output

// ── Energy audit ────────────────────────────────────────────────────────
//
// Tracks power flow through the drivetrain every tick and integrates to
// get cumulative energy [J] in each pathway. The balance check:
//
//   energy_in - energy_drag - energy_brake - energy_tire_slip ≈ ΔKE
//
// Any persistent drift means the integrator is creating or destroying energy.
// Reported as energy_balance = (energy_in - losses) - ΔKE. Should stay near 0.
//
// Flywheel KE (0.5 * I_engine * omega_engine²) is included in ke_current.
// The reflected inertia approach means wheel KE uses bare I_wheel and the
// flywheel energy is tracked separately.

struct EnergyAuditor {
    // Cumulative energy through each pathway [J]
    double energy_engine    = 0;  // drive torque × wheel omega (signed; negative = engine braking)
    double energy_drag      = 0;  // aero drag power dissipated
    double energy_brake     = 0;  // brake torque × |wheel omega|
    double energy_body      = 0;  // F_body · V_body + Mz · yaw_rate (what the integrator applies)

    // Initial KE (set on first call)
    double ke_initial = 0;
    double ke_body_initial = 0;
    bool   initialized = false;

    // Per-frame values for CSV
    double ke_current       = 0;
    double balance          = 0;  // integrator error: energy_body - ΔKE
    double tire_loss        = 0;  // derived: engine - drag - brake - body power

    void update(const Vehicle& car, const VehicleInput& input,
                const VehicleForces& forces, double dt) {
        const auto& s = car.state;
        const auto& p = car.params;
        constexpr double I_wheel = 1.5;

        // ── Kinetic energy, split into body and wheels ──
        double ke_body = 0.5 * p.mass * s.velocity.length_sq()
                       + 0.5 * p.Izz * s.yaw_rate * s.yaw_rate;
        double ke_wheels = 0;
        for (int i = 0; i < 4; ++i)
            ke_wheels += 0.5 * I_wheel * s.tires[i].omega * s.tires[i].omega;

        // Flywheel KE. Since we use bare I_wheel in the wheel KE calculation
        // (not I_eff), we need to add the flywheel's rotational energy
        // separately. When the clutch is engaged, omega_engine = omega_wheel
        // × ratio, so this captures the reflected energy correctly. When the
        // clutch is open, the flywheel has its own independent speed.
        double ke_flywheel = 0.5 * s.drivetrain.I_engine
                           * s.drivetrain.omega_engine
                           * s.drivetrain.omega_engine;
        ke_current = ke_body + ke_wheels + ke_flywheel;

        if (!initialized) {
            ke_initial = ke_current;
            ke_body_initial = ke_body;
            initialized = true;
            return;
        }

        // ── Engine power: drive torque at each driven wheel × omega ──
        double drive_torque = s.drivetrain.get_drive_torque(input.throttle);
        for (int i = RL; i <= RR; ++i)
            energy_engine += (drive_torque / 2.0) * s.tires[i].omega * dt;

        // ── Aero drag dissipation: Fd · V ──
        double spd = s.speed_mps();
        if (spd > 0.1) {
            double Fd = 0.5 * p.air_density * p.drag_coefficient
                      * p.frontal_area * spd * spd;
            energy_drag += Fd * spd * dt;
        }

        // ── Brake dissipation: brake torque × |omega| per wheel ──
        for (int i = 0; i < 4; ++i) {
            bool front = (i == FL || i == FR);
            double bt = s.drivetrain.get_brake_torque_per_wheel(input.brake, front);
            energy_brake += bt * std::abs(s.tires[i].omega) * dt;
        }

        // ── Body power: what the integrator applies to the rigid body ──
        // P = F_body · V_body + M_yaw · yaw_rate
        double p_body = forces.Fx_body * s.Vx_body
                      + forces.Fy_body * s.Vy_body
                      + forces.yaw_torque * s.yaw_rate;
        energy_body += p_body * dt;

        // ── Derived quantities ──
        tire_loss = energy_engine - energy_drag - energy_brake - energy_body;

        // Integrator error: compare body power integral to actual body ΔKE.
        // These use the same force and the same velocity — the only difference
        // is whether you integrate P·dt or look at 0.5·m·V². Any gap is
        // numerical error from the semi-implicit Euler scheme.
        double delta_ke_body = ke_body - ke_body_initial;
        balance = energy_body - delta_ke_body;
    }

    void report() const {
        double delta_ke = ke_current - ke_initial;
        std::cerr << "\n── energy audit ──\n"
                  << "  KE initial:      " << ke_initial / 1000 << " kJ\n"
                  << "  KE final:        " << ke_current / 1000 << " kJ\n"
                  << "  ΔKE:             " << delta_ke / 1000 << " kJ\n"
                  << "  engine in:       " << energy_engine / 1000 << " kJ\n"
                  << "  losses:\n"
                  << "    aero drag:     " << energy_drag / 1000 << " kJ\n"
                  << "    brakes:        " << energy_brake / 1000 << " kJ\n"
                  << "    tire/rr:       " << tire_loss / 1000 << " kJ"
                  << "  (slip + rolling resistance)\n"
                  << "  body power:      " << energy_body / 1000 << " kJ"
                  << "  (applied to rigid body)\n"
                  << "  integrator err:  " << balance / 1000 << " kJ"
                  << "  (" << (std::abs(energy_body) > 1
                     ? balance / energy_body * 100 : 0)
                  << "% of body power)\n";
    }
};

// ── CSV logging ─────────────────────────────────────────────────────────

void csv_header() {
    std::cout << "time,phase,"
        "pos_x,pos_y,yaw_deg,"
        "speed_kph,vx_body,vy_body,yaw_rate_dps,"
        "throttle,brake,steer_deg,gear,rpm,"
        "fz_fl,fz_fr,fz_rl,fz_rr,"
        "sr_fl,sr_fr,sr_rl,sr_rr,"
        "sa_fl_deg,sa_fr_deg,sa_rl_deg,sa_rr_deg,"
        "fx_fl,fx_fr,fx_rl,fx_rr,"
        "fy_fl,fy_fr,fy_rl,fy_rr,"
        "lat_accel_g,"
        "susp_defl_fl_mm,susp_defl_fr_mm,susp_defl_rl_mm,susp_defl_rr_mm,"
        "susp_vel_fl,susp_vel_fr,susp_vel_rl,susp_vel_rr,"
        "roll_deg,pitch_deg,"
        // ── new diagnostic columns ──
        "fz_total,understeer_deg,"
        "fc_util_fl,fc_util_fr,fc_util_rl,fc_util_rr,"
        "radius_m,yaw_rate_ref_dps,"
        // ── energy audit columns ──
        "ke_total_kj,energy_engine_kj,energy_drag_kj,energy_brake_kj,"
        "energy_tire_loss_kj,energy_body_kj,integrator_err_kj\n";
}

void csv_row(double t, const std::string& phase, const VehicleInput& in,
             const VehicleState& s, double lat_g,
             const EnergyAuditor& energy,
             double ref_yaw_rate = 0.0) {
    constexpr double R2D = 180.0 / M_PI;
    auto& w = s.tires;
    auto& su = s.suspension;

    // Diagnostic: total Fz (should ≈ m*g = 13734 N)
    double fz_total = w[FL].normal_load + w[FR].normal_load
                    + w[RL].normal_load + w[RR].normal_load;

    // Understeer angle: avg front slip angle − avg rear slip angle
    // Positive = understeer (front saturating first)
    double sa_front_avg = (w[FL].slip_angle + w[FR].slip_angle) / 2.0;
    double sa_rear_avg  = (w[RL].slip_angle + w[RR].slip_angle) / 2.0;
    double understeer_deg = (std::abs(sa_front_avg) - std::abs(sa_rear_avg)) * R2D;

    // Friction circle utilization per corner: sqrt(Fx²+Fy²) / (mu * Fz)
    // Should never exceed 1.0 — if it does, the friction circle clamp leaked
    auto fc_util = [](const TireState& t, double mu) -> double {
        double fz = std::max(t.normal_load, 1.0);
        return std::sqrt(t.Fx * t.Fx + t.Fy * t.Fy) / (mu * fz);
    };
    double mu = 1.0;  // matches tire_params.mu_peak default
    double fc_fl = fc_util(w[FL], mu);
    double fc_fr = fc_util(w[FR], mu);
    double fc_rl = fc_util(w[RL], mu);
    double fc_rr = fc_util(w[RR], mu);

    // Instantaneous radius from speed / yaw_rate (0 if straight)
    double radius = 0.0;
    if (std::abs(s.yaw_rate) > 1e-4)
        radius = std::abs(s.Vx_body / s.yaw_rate);

    std::cout << std::fixed << std::setprecision(4)
        << t << "," << phase << ","
        << s.position.x << "," << s.position.y << "," << s.yaw * R2D << ","
        << s.speed_kph() << "," << s.Vx_body << "," << s.Vy_body << ","
        << s.yaw_rate * R2D << ","
        << in.throttle << "," << in.brake << ","
        << in.steer_angle * R2D << ","
        << s.drivetrain.gearbox.current_gear + 1 << "," << s.drivetrain.engine_rpm << ","
        << w[FL].normal_load << "," << w[FR].normal_load << ","
        << w[RL].normal_load << "," << w[RR].normal_load << ","
        << w[FL].slip_ratio << "," << w[FR].slip_ratio << ","
        << w[RL].slip_ratio << "," << w[RR].slip_ratio << ","
        << w[FL].slip_angle*R2D << "," << w[FR].slip_angle*R2D << ","
        << w[RL].slip_angle*R2D << "," << w[RR].slip_angle*R2D << ","
        << w[FL].Fx << "," << w[FR].Fx << "," << w[RL].Fx << "," << w[RR].Fx << ","
        << w[FL].Fy << "," << w[FR].Fy << "," << w[RL].Fy << "," << w[RR].Fy << ","
        << lat_g << ","
        << su.corners[FL].deflection * 1000 << ","
        << su.corners[FR].deflection * 1000 << ","
        << su.corners[RL].deflection * 1000 << ","
        << su.corners[RR].deflection * 1000 << ","
        << su.corners[FL].velocity * 1000 << ","
        << su.corners[FR].velocity * 1000 << ","
        << su.corners[RL].velocity * 1000 << ","
        << su.corners[RR].velocity * 1000 << ","
        << su.roll_angle * R2D << ","
        << su.pitch_angle * R2D << ","
        // new columns
        << fz_total << "," << understeer_deg << ","
        << fc_fl << "," << fc_fr << "," << fc_rl << "," << fc_rr << ","
        << radius << "," << ref_yaw_rate * R2D << ","
        // energy audit (kJ for human readability)
        << energy.ke_current / 1000 << ","
        << energy.energy_engine / 1000 << ","
        << energy.energy_drag / 1000 << ","
        << energy.energy_brake / 1000 << ","
        << energy.tire_loss / 1000 << ","
        << energy.energy_body / 1000 << ","
        << (energy.energy_body - (energy.ke_current - energy.ke_initial)) / 1000 << "\n";
}

// ── Scenario infrastructure ─────────────────────────────────────────────

struct ScenarioFrame {
    VehicleInput input;
    std::string  phase;
    double       ref_yaw_rate;  // analytical reference (0 if N/A)
};

// ── 1. Original mixed scenario (unchanged) ──────────────────────────────
//
//   0-8s    full throttle from standstill (longitudinal grip, gear shifts)
//   8-10s   brake to cornering speed
//   10-20s  steady left turn at ~1.5° steer, partial throttle (~0.4-0.5g)
//   20-25s  trail brake while unwinding steering

constexpr double MIXED_DURATION = 25.0;

ScenarioFrame scenario_mixed(double t, const VehicleState&) {
    constexpr double D2R = M_PI / 180.0;
    if (t < 8.0) {
        return {{std::min(1.0, t / 1.5), 0, 0}, "accel", 0};
    }
    if (t < 10.0) {
        return {{0, std::min(0.5, (t - 8) / 0.5), 0}, "brake", 0};
    }
    if (t < 20.0) {
        double ramp = std::min(1.0, (t - 10) / 1.5);
        return {{0.35, 0, 1.5 * ramp * D2R}, "corner", 0};
    }
    double brake = std::min(0.25, (t - 20) / 2.5);
    double steer = 1.5 * std::max(0.0, 1.0 - (t - 20) / 3.5);
    return {{0, brake, steer * D2R}, "trail_brake", 0};
}

// ── 2. Constant-radius skidpad ──────────────────────────────────────────
//
// The standard vehicle dynamics characterization test. Drive around a circle
// of fixed radius at slowly increasing speed. At each speed the bicycle model
// predicts a required steer angle:
//
//   delta = L / R + K_us * a_lat
//
// where K_us is the understeer gradient. If K_us > 0 (front-heavy or more
// front roll stiffness), steering increases with speed. If K_us < 0, the car
// oversteers. At the limit, slip angles and lat_g tell you where the tires
// are on the Pacejka curve.
//
// We use a PD controller on yaw rate to hold the target radius rather than
// open-loop steer, because the point is to measure what steer angle the
// *sim* needs — that's the data you compare against the bicycle model.
//
// Strategy:
//   0-5s     accelerate straight to ~40 kph
//   5-8s     ramp in steering to establish the circle
//   8-50s    hold radius, slowly ramp throttle → speed climbs from ~40 to ~90+ kph
//            the controller adjusts steer to hold R = 50m
//   at each logged frame we also output the bicycle model's predicted yaw rate
//   for comparison: omega_ref = V / R

constexpr double SKIDPAD_DURATION = 50.0;
constexpr double SKIDPAD_RADIUS   = 50.0;

ScenarioFrame scenario_skidpad(double t, const VehicleState& s) {
    constexpr double D2R = M_PI / 180.0;

    // Phase 1: straight-line acceleration
    if (t < 5.0) {
        double throttle = std::min(0.6, t / 2.0);
        return {{throttle, 0, 0}, "accel", 0};
    }

    // Target yaw rate for constant radius R at current speed
    double V = std::max(s.Vx_body, 2.0);
    double omega_target = V / SKIDPAD_RADIUS;

    // Bicycle model reference steer angle (no understeer gradient — kinematic)
    // delta_kin = L / R, with L = 2.6m, R = 50m → ~2.98°
    constexpr double L = 2.6;
    double delta_kinematic = L / SKIDPAD_RADIUS;

    // PD controller on yaw rate error → steer correction
    // Kp chosen so that 1 deg/s yaw rate error ≈ 0.3° correction
    // Kd for damping oscillation
    double omega_error = omega_target - s.yaw_rate;
    constexpr double Kp = 18.0;   // [rad steer / (rad/s yaw error)]
    constexpr double Kd = 2.0;

    // We don't have previous error stored, so use a simple proportional +
    // feed-forward. The feed-forward (delta_kinematic) does most of the work.
    double steer = delta_kinematic + Kp * omega_error;

    // Phase 2: ramp in the steering over 3 seconds
    if (t < 8.0) {
        double ramp = (t - 5.0) / 3.0;
        steer *= ramp;
    }

    // Clamp steer to something physical (±10°)
    steer = std::clamp(steer, -10.0 * D2R, 10.0 * D2R);

    // Throttle: ramp from 0.3 to 0.7 over the hold phase to push speed up
    double throttle;
    if (t < 8.0) {
        throttle = 0.3;
    } else {
        double progress = (t - 8.0) / (SKIDPAD_DURATION - 8.0);
        throttle = 0.3 + 0.4 * progress;
    }

    std::string phase = (t < 8.0) ? "entry" : "hold";
    return {{throttle, 0, steer}, phase, omega_target};
}

// ── 3. Straight-line brake to stop ──────────────────────────────────────
//
// Full throttle from standstill to ~100 kph, then hard braking to a stop.
// Known analytical checks:
//
//   - During steady braking, front Fz should increase by m*a*h/L per axle
//     and rear Fz should decrease by the same amount.
//   - Stopping distance from V at deceleration a: d = V² / (2*a).
//     At ~100 kph (27.8 m/s) and ~0.8g (7.85 m/s²): d ≈ 49m.
//   - Total Fz should remain m*g throughout.
//   - All slip ratios should be negative during braking; slip angles ≈ 0.
//   - FL should exactly equal FR, RL should exactly equal RR (symmetry check).

constexpr double BRAKE_STOP_DURATION = 18.0;

ScenarioFrame scenario_brake_stop(double t, const VehicleState& s) {
    // Phase 1: launch (same as mixed, but only to ~100 kph)
    if (t < 7.0) {
        double throttle = std::min(1.0, t / 1.5);
        return {{throttle, 0, 0}, "accel", 0};
    }

    // Phase 2: hard braking to standstill
    // Ramp brake to 0.8 over 0.3s to avoid a discontinuity
    double brake = std::min(0.8, (t - 7.0) / 0.3);

    // Stop the sim once the car is nearly stopped
    std::string phase = (s.speed_kph() < 1.0 && t > 8.0) ? "stopped" : "brake";

    return {{0, brake, 0}, phase, 0};
}

// ── 4. Step steer ───────────────────────────────────────────────────────
//
// Cruise at constant speed (~80 kph), then apply a step steer input.
// The transient yaw rate response should show:
//   - A first-order rise with time constant set by tire relaxation length
//     (tau ≈ relaxation_length / Vx ≈ 0.4 / 22 ≈ 18ms tire lag, plus
//      suspension transient ~100ms, so total rise ~100-200ms)
//   - Possible overshoot if underdamped
//   - Steady state yaw rate matching the bicycle model:
//     omega_ss = V * delta / (L + K_us * V²)
//     For small K_us and delta = 1.5°: omega_ss ≈ V * delta / L
//     = 22.2 * 0.0262 / 2.6 ≈ 0.224 rad/s ≈ 12.8 °/s
//
// The step should be a clean discontinuity (no ramp) so you can measure
// the response time precisely. We do ramp throttle to hold speed.

constexpr double STEP_STEER_DURATION = 20.0;
constexpr double STEP_STEER_CRUISE_SPEED_MPS = 22.2;  // ~80 kph
constexpr double STEP_STEER_ANGLE_DEG = 1.5;

ScenarioFrame scenario_step_steer(double t, const VehicleState& s) {
    constexpr double D2R = M_PI / 180.0;
    constexpr double L = 2.6;
    double delta = STEP_STEER_ANGLE_DEG * D2R;

    // Phase 1: accelerate to cruise speed (0 - 8s)
    if (t < 8.0) {
        double throttle = std::min(1.0, t / 1.5);
        return {{throttle, 0, 0}, "accel", 0};
    }

    // Phase 2: cruise at ~80 kph (8 - 10s). Use a simple P controller on
    // speed to hold it. This gives the car time to settle before the step.
    if (t < 10.0) {
        double speed_error = STEP_STEER_CRUISE_SPEED_MPS - s.Vx_body;
        double throttle = std::clamp(0.3 + speed_error * 0.5, 0.0, 1.0);
        return {{throttle, 0, 0}, "cruise", 0};
    }

    // Phase 3: step steer. Apply full steer angle instantaneously at t=10s.
    // Continue holding speed with throttle controller.
    double speed_error = STEP_STEER_CRUISE_SPEED_MPS - s.Vx_body;
    double throttle = std::clamp(0.3 + speed_error * 0.5, 0.0, 1.0);

    // Bicycle model steady-state reference yaw rate: omega = V * delta / L
    // (ignoring understeer gradient — this is the kinematic prediction)
    double omega_ref = s.Vx_body * delta / L;

    return {{throttle, 0, delta}, "step", omega_ref};
}

// ── Scenario dispatch ───────────────────────────────────────────────────

enum ScenarioType { MIXED, SKIDPAD, BRAKE_STOP, STEP_STEER };

struct ScenarioConfig {
    ScenarioType type;
    double duration;
    const char* name;
};

ScenarioConfig parse_scenario(int argc, char* argv[]) {
    if (argc < 2 || std::strcmp(argv[1], "mixed") == 0)
        return {MIXED,      MIXED_DURATION,      "mixed"};
    if (std::strcmp(argv[1], "skidpad") == 0)
        return {SKIDPAD,    SKIDPAD_DURATION,    "skidpad"};
    if (std::strcmp(argv[1], "brake_stop") == 0)
        return {BRAKE_STOP, BRAKE_STOP_DURATION, "brake_stop"};
    if (std::strcmp(argv[1], "step_steer") == 0)
        return {STEP_STEER, STEP_STEER_DURATION, "step_steer"};

    std::cerr << "unknown scenario: " << argv[1] << "\n"
              << "available: mixed, skidpad, brake_stop, step_steer\n";
    return {MIXED, MIXED_DURATION, "mixed"};
}

ScenarioFrame dispatch(ScenarioType type, double t, const VehicleState& s) {
    switch (type) {
        case SKIDPAD:    return scenario_skidpad(t, s);
        case BRAKE_STOP: return scenario_brake_stop(t, s);
        case STEP_STEER: return scenario_step_steer(t, s);
        default:         return scenario_mixed(t, s);
    }
}

// ── Diagnostics summary to stderr ───────────────────────────────────────

struct DiagnosticAccumulator {
    double fz_sum_min =  1e9;
    double fz_sum_max = -1e9;
    double fc_util_max = 0;
    int    fc_violations = 0;     // frames where any corner > 1.0
    double max_lr_asymmetry = 0;  // max |FL-FR| when steer=0 and speed>5
    int    frames = 0;

    void update(const VehicleState& s, const VehicleInput& in) {
        auto& w = s.tires;
        double fz_sum = w[FL].normal_load + w[FR].normal_load
                      + w[RL].normal_load + w[RR].normal_load;
        fz_sum_min = std::min(fz_sum_min, fz_sum);
        fz_sum_max = std::max(fz_sum_max, fz_sum);

        double mu = 1.0;
        for (int i = 0; i < 4; ++i) {
            double fz = std::max(w[i].normal_load, 1.0);
            double util = std::sqrt(w[i].Fx * w[i].Fx + w[i].Fy * w[i].Fy) / (mu * fz);
            fc_util_max = std::max(fc_util_max, util);
            if (util > 1.001) fc_violations++;
        }

        // Symmetry check: when driving straight, L/R should be identical
        if (std::abs(in.steer_angle) < 1e-6 && s.Vx_body > 5.0) {
            double asym_front = std::abs(w[FL].normal_load - w[FR].normal_load);
            double asym_rear  = std::abs(w[RL].normal_load - w[RR].normal_load);
            max_lr_asymmetry = std::max({max_lr_asymmetry, asym_front, asym_rear});
        }

        frames++;
    }

    void report(double mass) const {
        double mg = mass * 9.81;
        std::cerr << "\n── diagnostics ──\n"
                  << "  Fz total: [" << fz_sum_min << ", " << fz_sum_max << "] N"
                  << "  (m*g = " << mg << ")\n"
                  << "  Fz deviation: "
                  << (std::max(std::abs(fz_sum_min - mg), std::abs(fz_sum_max - mg)) / mg * 100)
                  << "%\n"
                  << "  friction circle max utilisation: " << fc_util_max
                  << (fc_util_max > 1.001 ? "  *** EXCEEDED ***" : "  ok") << "\n"
                  << "  FC violation frames: " << fc_violations
                  << " / " << (frames * 4) << " corner-frames\n"
                  << "  max L/R asymmetry (straight line): " << max_lr_asymmetry << " N\n";
    }
};

// ── Main ────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
    auto cfg = parse_scenario(argc, argv);

    std::cerr << "eudromos — vehicle physics sim\n"
              << "scenario: " << cfg.name << "\n"
              << DT * 1000 << " ms timestep, " << cfg.duration << " s duration\n"
              << "CSV on stdout\n";

    Vehicle car;
    car.init();
    csv_header();

    DiagnosticAccumulator diag;
    EnergyAuditor energy;
    double lat_g = 0;
    int tick = 0;
    bool stopped = false;

    for (double t = 0; t <= cfg.duration && !stopped; t += DT, ++tick) {
        auto [input, phase, ref_yaw_rate] = dispatch(cfg.type, t, car.state);

        if (tick % LOG_SKIP == 0)
            csv_row(t, phase, input, car.state, lat_g, energy, ref_yaw_rate);

        VehicleForces forces = car.compute_forces(input, DT);
        lat_g = forces.lateral_accel_g;
        car.integrate(forces, DT);

        energy.update(car, input, forces, DT);
        diag.update(car.state, input);

        // Early termination for brake_stop once the car is stopped
        if (cfg.type == BRAKE_STOP && phase == "stopped")
            stopped = true;
    }

    auto& su = car.state.suspension;
    std::cerr << "\nfinal: " << car.state.position
              << "  " << car.state.speed_kph() << " kph"
              << "  heading " << car.state.yaw * 180 / M_PI << "°"
              << "  gear " << car.state.drivetrain.gearbox.current_gear + 1
              << "  roll " << su.roll_angle * 180 / M_PI << "°"
              << "  pitch " << su.pitch_angle * 180 / M_PI << "°"
              << "  ticks " << tick << "\n";

    diag.report(car.params.mass);
    energy.report();
}