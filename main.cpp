// main.cpp — runs the sim, dumps CSV.
//
// g++ -std=c++17 -O2 -o vehicle_sim main.cpp -lm
// ./vehicle_sim > output.csv

#include "vehicle.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

constexpr double DT       = 1.0 / 120.0;  // 120 Hz
constexpr double DURATION = 25.0;
constexpr int    LOG_SKIP = 12;            // log every 12th tick → ~10 Hz output

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
        "roll_deg,pitch_deg\n";
}

void csv_row(double t, const std::string& phase, const VehicleInput& in,
             const VehicleState& s, double lat_g) {
    constexpr double R2D = 180.0 / M_PI;
    auto& w = s.tires;
    auto& su = s.suspension;
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
        << su.pitch_angle * R2D << "\n";
}

// Scripted scenario. Four phases that exercise different parts of the tire model:
//
//   0-8s    full throttle from standstill (longitudinal grip, gear shifts)
//   8-10s   brake to cornering speed
//   10-20s  steady left turn at ~1.5° steer, partial throttle (~0.4-0.5g)
//   20-25s  trail brake while unwinding steering
//
// Steer angle chosen from bicycle model: at 80 kph targeting ~0.6g lateral,
// delta ≈ atan(a_lat * L / V²) ≈ 1.8° → we use 1.5° with a ramp-in.

struct Scenario { VehicleInput input; std::string phase; };

Scenario scenario(double t) {
    constexpr double D2R = M_PI / 180.0;
    if (t < 8.0) {
        // ramp throttle over 1.5s to limit launch wheelspin
        return {{std::min(1.0, t / 1.5), 0, 0}, "accel"};
    }
    if (t < 10.0) {
        return {{0, std::min(0.5, (t - 8) / 0.5), 0}, "brake"};
    }
    if (t < 20.0) {
        double ramp = std::min(1.0, (t - 10) / 1.5);
        return {{0.35, 0, 1.5 * ramp * D2R}, "corner"};
    }
    // trail brake: gentle brakes, gradually straighten steering
    double brake = std::min(0.25, (t - 20) / 2.5);
    double steer = 1.5 * std::max(0.0, 1.0 - (t - 20) / 3.5);
    return {{0, brake, steer * D2R}, "trail_brake"};
}

int main() {
    std::cerr << "eudromos — vehicle physics sim\n"
              << DT * 1000 << " ms timestep, " << DURATION << " s duration\n"
              << "CSV on stdout\n";

    Vehicle car;
    car.init();
    csv_header();

    double lat_g = 0;  // from previous frame's forces
    int tick = 0;
    for (double t = 0; t <= DURATION; t += DT, ++tick) {
        auto [input, phase] = scenario(t);

        if (tick % LOG_SKIP == 0)
            csv_row(t, phase, input, car.state, lat_g);

        VehicleForces forces = car.compute_forces(input, DT);
        lat_g = forces.lateral_accel_g;
        car.integrate(forces, DT);
    }

    auto& su = car.state.suspension;
    std::cerr << "\nfinal: " << car.state.position
              << "  " << car.state.speed_kph() << " kph"
              << "  heading " << car.state.yaw * 180 / M_PI << "°"
              << "  gear " << car.state.drivetrain.gearbox.current_gear + 1
              << "  roll " << su.roll_angle * 180 / M_PI << "°"
              << "  pitch " << su.pitch_angle * 180 / M_PI << "°"
              << "  ticks " << tick << "\n";
}