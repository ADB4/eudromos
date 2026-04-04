# eudromos

Tire-level vehicle physics. C++17, no dependencies, no game engine.

I wanted to actually understand what's happening between rubber and asphalt
before handing it off to someone else's black box. This is heavily informed by
Brian Beckman's *The Physics of Racing*. if you haven't read it, go do that first.

## building

```
g++ -std=c++17 -O2 -o vehicle_sim main.cpp -lm
./vehicle_sim > output.csv
```

Status/diagnostics go to stderr. CSV goes to stdout.

## what's in here

```
math_types.h    Vec3, Mat3 (hand-rolled, no Eigen)
tire.h          Pacejka Magic Formula, slip ratio/angle, friction circle
drivetrain.h    engine torque curve, 6-speed box, RWD
vehicle.h       rigid body + 4 wheels + weight transfer + force accumulation
main.cpp        120 Hz fixed-step loop, scenario script, CSV logger
```

The sim runs a scripted 25-second scenario: launch from standstill, brake into
a corner, hold a steady-state left turn, then trail-brake out of it. Nothing
is interactive yet. The point is to get numbers you can plot and sanity-check.

## does it work?

Roughly, yes. During steady cornering the tires sit at 1–3 degree slip angle, lateral
accel is 0.3–0.5g, outside wheels pick up load, inside wheels unload. Trail
braking shifts weight forward and the rear slip angles creep up (oversteer
tendency). The Pacejka curves peak where they should: ~15% slip ratio
longitudinal, ~8 degree slip angle lateral.

The numbers aren't dyno-validated against a real car. They're *plausible*, the
right things go up when other things go down, and the magnitudes are in the
right ballpark for a 1400 kg sedan on street tires.

## what's not modeled (yet)

Everything is marked with `// SIMPLIFICATION:` in the code, but the big ones:

- No suspension at all. Weight transfer is instantaneous (infinitely stiff springs).
  This matters a lot for transient handling; it's probably the first thing to add.
- No aero downforce, only drag.
- Flat ground, yaw only; no pitch, no roll, no terrain.
- No clutch. Gears engage instantly.
- No diff. 50/50 torque split to the rear wheels.
- No drivetrain inertia (flywheel, driveshaft).
- Equal brake torque all four corners (real cars bias ~60% front).
- Friction circle is circular. Real tires are slightly elliptical.

## bugs I hit along the way

Worth documenting because these are easy to walk into:

1. **Fy sign convention.** Pacejka returns positive force for positive slip angle,
   but lateral tire force *opposes* slip ; it's a restoring force. Miss the negation
   and the car spins instantly. This one cost me a while.

2. **World-frame vs body-frame integration.** If you integrate velocity in the
   rotating body frame, you need the centripetal coupling terms (dVy/dt includes
   −Vx·ω). Skip them and lateral velocity grows without bound during cornering.
   Cleaner approach: accumulate forces in body frame, integrate in world frame.

3. **Tire relaxation length.** Without a first-order lag on Fy
   (τ = relaxation_length / |Vx|, typically ~0.4m), the tire produces full
   lateral force instantaneously and the yaw response overshoots violently.
   This is the difference between "corners" and "spins."

## CSV columns

```
time, phase, pos_x, pos_y, yaw_deg, speed_kph, vx_body, vy_body, yaw_rate_dps,
throttle, brake, steer_deg, gear, rpm, fz_fl, fz_fr, fz_rl, fz_rr,
sr_fl, sr_fr, sr_rl, sr_rr, sa_fl_deg, sa_fr_deg, sa_rl_deg, sa_rr_deg,
fx_fl, fx_fr, fx_rl, fx_rr, fy_fl, fy_fr, fy_rl, fy_rr, lat_accel_g
```

## quick plot

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('output.csv')
fig, ax = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

ax[0].plot(df.time, df.speed_kph)
ax[0].set_ylabel('speed (kph)')

ax[1].plot(df.time, df.sa_fl_deg, label='FL')
ax[1].plot(df.time, df.sa_rl_deg, label='RL')
ax[1].set_ylabel('slip angle (°)')
ax[1].legend()

ax[2].plot(df.time, df.fz_fl, label='FL')
ax[2].plot(df.time, df.fz_rr, label='RR')
ax[2].set_ylabel('Fz (N)')
ax[2].legend()

ax[3].plot(df.time, df.lat_accel_g)
ax[3].set_ylabel('lat accel (g)')
ax[3].set_xlabel('time (s)')

plt.tight_layout()
plt.savefig('sim_output.png', dpi=150)
```