# Accel-Mag-Calibrators

Portable, allocation-free sensor calibration for embedded flight control:
accelerometer, magnetometer and board-mounting alignment.

Three self-contained C++11 classes that turn raw IMU and magnetometer readings
into corrected ones, together with the minimal vector/matrix types they need.
They were extracted from an STM32F7 flight-controller project and deliberately
carry **no HAL, no board headers and no peripheral dependencies**, so the same
code builds unchanged on another MCU or in a host-side test harness.

---

## Contents

| File | What it is |
|---|---|
| `AccelerometerCalibrator.hpp/.cpp` | Accelerometer calibration: **six-position**, fitting bias and per-axis scale |
| `CompassCalibrator.hpp/.cpp` | Magnetometer calibration: hard-iron offset + soft-iron matrix from an ellipsoid fit over a free-hand sweep |
| `LevelCalibrator.hpp/.cpp` | Board-mounting rotation (roll/pitch) against a known-level surface, plus an optional supplied yaw offset |
| `MathTypes.hpp` | Aggregator header — the maths types and nothing else |
| `Vector3f.hpp`, `Matrix3f.hpp`, `Quaternionf.hpp` | Minimal float vector, 3×3 matrix and quaternion |

Comments in the sources refer to `ICM42688P`, `LIS3MDL`, `Imu` and
`Calibrator::correctBoardFrame()`. Those live in the parent flight-controller
project; this repository is the portable algorithm subset and does not depend
on them.

---

## Design constraints

These are the rules the code is written to, and they explain most of its API:

- **C++11 freestanding.** Needs only `<cmath>` and `<cstring>`.
- **No dynamic allocation, no exceptions, no RTTI.**
- **No HAL.** Including a shared `common.hpp` for `Vector3f` alone used to pull
  `stm32f7xx_hal.h` in behind it, which made every one of these classes
  unbuildable elsewhere. `MathTypes.hpp` exists to break that.
- **Caller-owned sample buffers.** The compass procedure needs
  ~3.6 kB of sample storage. That is passed in rather than held as a member
  array so the integrator decides where it lives — a scratch/CCM region, a
  stack buffer in the calibration routine, or memory shared between the two
  (they are mutually exclusive, so one buffer serves both). Two classes each
  owning their own array would put 7.2 kB permanently in `.bss` for procedures
  that run for a few seconds on the ground and never again.
- **Unit-agnostic.** Nothing assumes m/s², Gauss or raw LSB. Every threshold is
  either relative or expressed in the same units as the samples you feed in.
- **Single-threaded, non-blocking.** Feed one sample per control-loop tick;
  every call returns immediately. There are no timeouts and no internal delays.

---

## What each procedure can and cannot remove

An accelerometer reading is `raw = A·g + bias`, and the gain matrix splits as
`A = R·S`. Which part of that a procedure recovers is the single most important
thing to understand before choosing one:

| Error | Six-position | Level |
|---|---|---|
| Bias (offset) | ✅ | — |
| Per-axis scale | ✅ | — |
| `S` — cross-axis sensitivity (the chip's own axes not quite perpendicular) | ❌ | — |
| `R` — board mounted rotated in the airframe (roll/pitch) | ❌ | ✅ |
| Mounting yaw | ❌ | ⚠️ only if supplied |

The reason the accelerometer fit cannot recover `R` is inherent, not an
implementation gap: it only ever observes the **magnitude** of gravity, and
rotating a sphere leaves the same sphere, so the data carries no information
about `R` at all. No fit of that shape can, whatever its form. Measured, a 1.16°
mounting rotation comes out at 1.16°. Recovering it needs an outside reference,
which is what `LevelCalibrator` provides.

That split is deliberate and worth preserving: **the accelerometer correction is
diagonal, so it contains no rotation** and cannot fight the one `LevelCalibrator`
stores. Two rotations measured against different references would each be
correcting the other's reference.

`S` is left uncorrected. Separating it needs an ellipsoid fit over a free-hand
tumble — a 3.6 kB buffer, a minute of the operator's time, and a procedure whose
result depends on how well the airframe was rotated. On a modern MEMS part the
term is small enough that this was judged not to pay for itself; a sphere plot
of corrected samples is the cheap way to confirm that on your own hardware.

- Six-position leaves cross-axis misalignment essentially 1:1 (0.009° with none
  present, 1.15° with 1.15° of it).
- Level, run afterwards, removes the roll/pitch part of `R`. A mounting yaw left
  at 0 leaks back into roll/pitch as the airframe tilts (1° of mounting yaw
  leaves 0.50° yaw and 0.39° roll/pitch error across a ±34° envelope, against
  1.10° uncorrected).

**Recommended order:** accelerometer → level → compass. Level consumes
accelerometer-corrected samples, and its rotation applies to every sensor on the
board.

---

## Quick start

### Accelerometer — six-position

No sample buffer required; the procedure keeps six running averages.

```cpp
AccelerometerCalibrator accel;
accel.beginSixPosition(0.5f);          // motion threshold, sample units

for (int p = 0; p < (int)AccelPosition::NUM_POSITIONS; p++) {
    // Prompt the operator, wait for their "ready", then:
    accel.startPosition((AccelPosition)p);
    while (!accel.isPositionDone((AccelPosition)p)) {
        accel.addSample(imu.ax, imu.ay, imu.az);   // one per control tick
        if (accel.isStalled()) { /* "HOLD THE AIRFRAME STILL" */ }
    }
}

if (accel.calibrate() == AccelCalStatus::SUCCESS) {
    Vector3f bias   = accel.getBias();
    Vector3f scale  = accel.getScale();
    Matrix3f matrix = accel.getMatrix();   // diag(1/scale), for one stored form
    Vector3f corrected = accel.correct(raw);            // magnitude ~1.0
}
```

Each position needs `ACCEL_CAL_SAMPLES_PER_POSITION` (100) **consecutive**
samples inside the motion threshold; a single excursion discards that
position's average and restarts it.

### Compass

```cpp
Vector3f buffer[COMPASS_CAL_MAX_SAMPLES];             // 300 → 3.6 kB, caller-owned
CompassCalibrator mag;

if (!mag.begin(0.48f, buffer, COMPASS_CAL_MAX_SAMPLES)) {  // nominal field magnitude
    // null buffer, or capacity < COMPASS_CAL_MIN_SAMPLES — never started
}

while (!mag.isReadyToCalibrate()) {
    mag.addSample(lis.mx, lis.my, lis.mz);
    if (mag.isStalled()) {
        // "KEEP TURNING THE AIRFRAME — ALL SIDES, INCLUDING INVERTED"
    }
}

if (mag.calibrate() == CalStatus::SUCCESS) {
    Vector3f offset   = mag.getOffset();              // hard iron
    Matrix3f softiron = mag.getSoftIronMatrix();
    Vector3f corrected = mag.correct(raw);            // magnitude ~ nominal field
}
```

No pausing needed — the airframe is turned continuously, and it genuinely has
to be turned **over**, not merely waved about upright (see *Coverage gates*).

### Level

```cpp
LevelCalibrator lvl;
lvl.begin(Vector3f(0.0f, 0.0f, 1.0f),   // reading when level and upright (+Z up)
          9.80665f,                     // 1 g in sample units
          0.5f,                         // motion threshold
          0.0f);                        // known mounting yaw, degrees

while (!lvl.isReadyToCalibrate()) {
    lvl.addSample(accel.correct(raw));  // AFTER accelerometer calibration
}

if (lvl.calibrate() == LevelCalStatus::SUCCESS) {
    Matrix3f board = lvl.getRotation();  // apply to EVERY sensor on the board
    float tilt_deg = lvl.getTiltDeg();
}
```

Rest the airframe on a surface known to be level. `LEVEL_CAL_SAMPLES` (200)
consecutive still samples are averaged; one excursion past the motion threshold
discards the average and restarts it.

> The rotation is a rotation of the **whole board**, so it applies to
> accelerometer, gyroscope and magnetometer alike. Applying it to some sensors
> and not others leaves them disagreeing about which way the airframe points,
> which is worse than not correcting at all.

---

## Units and what `correct()` returns

Nothing in these classes assumes a unit system; you decide it by what you pass
as the nominal magnitude. The *output* scale differs between the classes, and
this catches people out:

| Call | Output |
|---|---|
| `AccelerometerCalibrator::correct()` | **Normalised to g** — magnitude ≈ 1.0 |
| `CompassCalibrator::correct()` | **Sensor field units** — magnitude ≈ the nominal magnitude passed to `begin()` |
| `LevelCalibrator::correct()` | Unchanged units — a pure rotation |

So `9.80665` (m/s²) and `2048` (raw LSB at ±16 g) are both valid `nominal_g`
values, as are `0.48` (Gauss) and `500` (raw LSB) for the field magnitude —
but the accelerometer hands you g either way, while the compass hands back
whatever scale you gave it. `CompassCalibrator::getNominalRadius()` returns the
magnitude its fit was normalised against, which is what you need to fold its
matrix into a single stored correction matrix.

`correct()` is a **pass-through until the fit has succeeded** on all three
classes — it returns its input unchanged rather than applying half-finished
gains.

---

## Coverage gates — why a fit gets rejected

The compass fit solves `x'Mx + 2n'x = 1` in the least-squares sense over a
9×9 normal-equation system, **pre-centred on the sample cloud's centroid**
(that form cannot represent an ellipsoid passing through the origin and is
ill-conditioned near it; pre-centring keeps the fit in the well-behaved regime
whatever the offset is, and the centre is added back at the end). The recovered
matrix is the symmetric, rotation-free square root of `M/K`, so it maps the
ellipsoid onto a sphere without also rotating the sensor frame.

A solution that satisfies every algebraic check can still be nonsense, so the
result passes through layered gates:

| Gate | Compass | Failure status |
|---|---|---|
| Minimum samples | 150 | `FAILED_NOT_ENOUGH_SAMPLES` |
| Filled spherical bins (of 64, ≤5 per bin) | ≥45 | `FAILED_POOR_COVERAGE` |
| Scatter anisotropy | ≥0.5 | `FAILED_POOR_COVERAGE` |
| 9×9 solve / 3×3 inverse | ✅ | `FAILED_SINGULAR_MATRIX` |
| Positive-definite, non-flat ellipsoid | ✅ | `FAILED_DEGENERATE_ELLIPSOID` |
| RMS fit residual | ≤15 % | `FAILED_POOR_FIT` |

`getLastFitResidual()` and `getScatterRatio()` report the last two as numbers
rather than just a pass/fail. That matters on a rejection: a residual of 0.16 is
a sweep that nearly worked and 0.90 is a sensor or an environment carrying no
usable field, and the difference decides what the operator should do next.

Two of those gates deserve explanation.

**Compass binning is done from the centre of the cloud, not the origin.** Hard
iron routinely exceeds Earth's field, and binning about the origin then
collapses every reading into the same cone — coverage saturates around 35 of 64
bins and the minimum can never be met. The centre is only knowable once samples
are in, so the bin table is rebuilt (`rebin()`) whenever the centroid has moved
by more than `COMPASS_CAL_REBIN_MOVE_FRACTION` of the nominal radius. Rebinning
**drops** samples the new binning makes redundant, so the sample count can step
backwards during collection; that is deliberate, and without it all 300 slots
fill with duplicates and coverage stalls one bin short of the minimum for ever.

**Scatter anisotropy is what catches a one-sided sweep.** Centroid-relative
binning weakens the coverage test: an operator who never inverts the airframe
covers barely half the sphere, the centroid sits inside that cap, directions
measured from it over-disperse, and the bin count still reaches the minimum.
Such a fit is exact on the samples collected and materially wrong on the
orientations that were not (measured at 6 % sensor noise: 3 % error in the
recovered centre, 5 % in corrected magnitude on orientations never visited) —
and `fitResidual()` cannot see it, because it only sees the samples that exist.
The smallest/largest eigenvalue ratio of the sample scatter matrix does catch
it, for one reason: **a scatter matrix is invariant to translation**, so hard
iron cannot influence it at all.

Measured values: full sphere 1.00; full sphere with heavy hard *and* soft iron
0.70; full sphere with very uneven dwell 0.81; two-thirds of a sphere 0.65 —
against 0.38 for a hemisphere plus a little, 0.25 for an exact hemisphere (the
analytic value: variance `r²/12` across the cap axis versus `r²/3` along it)
and 0.22 for a hemisphere with iron. The threshold of 0.5 sits in that gap with
room on both sides. It is part of *readiness*, not just a fit-time check, so a
lazy sweep simply keeps collecting instead of being told "ready" and then
rejected.

The six-position and level procedures have their own gates: a dead or stuck
axis (any per-axis sensitivity under a tenth of the largest) and positions
captured the wrong way round (a negative scale, which would otherwise "succeed"
while silently flying inverted) both give `FAILED_BAD_GEOMETRY`; a level
average whose magnitude is more than 15 % from nominal gives
`FAILED_BAD_MAGNITUDE`, and a mounting error beyond `LEVEL_CAL_MAX_TILT_DEG`
(15°) gives `FAILED_EXCESSIVE_TILT` — past that the likely explanation is a
surface that is not level or the wrong axis convention, and silently storing a
large rotation is the one outcome that would be hard to notice afterwards.

---

## Progress and stall reporting

`getProgressPercent()` returns the **worst** of the requirements, not their
average — all of them must be met, so the smallest is the honest answer to "how
far along is this?". For the compass that is samples, bins and scatter; for
six-position it is simply how many of the 600 required samples are in.

None of these procedures has a timeout. If the operator stops making progress —
vibration defeating the motion gate, an airframe never turned over — collection
simply continues for ever with the progress figure stuck. `isStalled()` is what
surfaces that. It reports true once the mode's limit of samples has arrived
without progress reaching a **new high** (a high-water mark, so the dip from a
motion reset or a rebin is not counted as lost progress). It is **advisory**:
collection continues and the flag clears by itself if progress resumes, so the
caller decides whether to abort and what to tell the operator.

| Procedure | Stall limit (samples) | ≈ time at 200 Hz | Longest gap measured in a successful run |
|---|---|---|---|
| Six-position | 10 000 | 50 s | — (a healthy run never gaps) |
| Compass | 6 000 | 30 s | 1 154 |
| Level | 10 000 | 50 s | — |

These are **sample counts**, and nothing in these classes measures time — so
what they mean in seconds is set entirely by how fast you feed them. The column
above assumes 200 Hz. Drive them from a 1 kHz control loop and every figure is
five times shorter: six-position's allowance becomes 10 seconds, which an
operator settling an airframe can exceed. Decimate the feed or rescale the
limits; the same applies to `LEVEL_CAL_SAMPLES` and
`ACCEL_CAL_SAMPLES_PER_POSITION`, which set averaging times the same way.

Six-position's allowance is really for an operator still settling the airframe
after sending READY; a healthy run never gaps at all, because progress climbs
with every accepted sample.

Vibration tolerances, measured:

- Six-position (0.5 default motion threshold): completes reliably up to
  0.12 m/s² RMS per axis, 1 run in 20 at 0.20, never at 0.30. A part at rest is
  far below this (ICM42688P is ~0.006 m/s²); a bench that shakes or props
  turning are not.

---

## Tuning constants

All are `#define`s at the top of the respective header, with the reasoning for
each value recorded beside it.

**`AccelerometerCalibrator.hpp`**

| Constant | Default | Meaning |
|---|---|---|
| `ACCEL_CAL_SAMPLES_PER_POSITION` | 100 | Consecutive still samples per face |
| `ACCEL_CAL_SIXPOS_STALL_LIMIT` | 10 000 | Stall detection |
| `ACCEL_CAL_MIN_RADIUS` | 1e-6 | Divide-by-zero guard on per-axis sensitivity |
| `ACCEL_CAL_STANDARD_GRAVITY` | 9.80665 | Standard gravity, for callers working in m/s² |

**`CompassCalibrator.hpp`**

| Constant | Default | Meaning |
|---|---|---|
| `COMPASS_CAL_MAX_SAMPLES` / `MIN_SAMPLES` | 300 / 150 | Recommended and minimum capacity |
| `COMPASS_CAL_NUM_BINS` / `MIN_BINS` / `MAX_PER_BIN` | 64 / 45 / 5 | Spherical coverage (>70 %) |
| `COMPASS_CAL_REBIN_MOVE_FRACTION` | 0.10 | Centroid movement that triggers a rebin |
| `COMPASS_CAL_MIN_SCATTER_RATIO` | 0.5 | One-sided-sweep rejection. Do not go below ~0.3 |
| `COMPASS_CAL_MAX_FIT_RESIDUAL` | 0.15 | Max RMS deviation from nominal, as a fraction |
| `COMPASS_CAL_STALL_LIMIT` | 6 000 | Stall detection |
| `COMPASS_CAL_MIN_RADIUS` | 1e-6 | Guard on nominal field magnitude |

**`LevelCalibrator.hpp`**

| Constant | Default | Meaning |
|---|---|---|
| `LEVEL_CAL_SAMPLES` | 200 | Consecutive still samples |
| `LEVEL_CAL_MOTION_THRESHOLD` | 0.5 | Motion gate, sample units |
| `LEVEL_CAL_MAX_MAGNITUDE_ERROR` | 0.15 | How far the average may be from 1 g |
| `LEVEL_CAL_MAX_TILT_DEG` | 15.0 | Largest mounting error accepted |
| `LEVEL_CAL_STALL_LIMIT` | 10 000 | Stall detection |

---

## Status and result codes

```
AccelCalStatus   IDLE = 0, IN_PROGRESS = 1, SUCCESS = 2,
                 FAILED_NOT_ENOUGH_SAMPLES = 3, FAILED_BAD_GEOMETRY = 4

AccelSampleResult ACCEPTED, ACCEPTED_POSITION_DONE, REJECTED_MOTION,
                 REJECTED_POSITION_DONE, REJECTED_NOT_STARTED

CalStatus        IDLE, COLLECTING, READY_TO_FIT, SUCCESS,
                 FAILED_NOT_ENOUGH_SAMPLES, FAILED_POOR_COVERAGE,
                 FAILED_SINGULAR_MATRIX, FAILED_DEGENERATE_ELLIPSOID,
                 FAILED_POOR_FIT

SampleResult     ACCEPTED, REJECTED_TOO_CLOSE, REJECTED_BUFFER_FULL,
                 REJECTED_NOT_COLLECTING

LevelCalStatus   IDLE, IN_PROGRESS, SUCCESS, FAILED_NOT_ENOUGH_SAMPLES,
                 FAILED_BAD_MAGNITUDE, FAILED_EXCESSIVE_TILT

LevelSampleResult ACCEPTED, ACCEPTED_DONE, REJECTED_MOTION,
                 REJECTED_DONE, REJECTED_NOT_STARTED
```

`AccelCalStatus` is numbered explicitly because these values are often reported
over a link. They were renumbered when the tumble mode was removed: the four
codes only an ellipsoid fit could produce went with it, and the two that remain
moved down.

`REJECTED_TOO_CLOSE` (compass) means the sample landed in a bin that is already full —
it is the normal signal that this orientation has enough coverage, not an
error. Compass `addSample()` deliberately checks the buffer **before** the
mode, so a caller can always tell "buffer full" from "never started".

---

## Memory footprint

Objects hold no sample storage of their own; the buffer is yours to place.

| Object | `sizeof` |
|---|---|
| `AccelerometerCalibrator` | 512 B |
| `CompassCalibrator` | 184 B |
| `LevelCalibrator` | 100 B |
| 300-sample buffer (`Vector3f[300]`) | 3 600 B |

Measured with `sizeof` on a 64-bit host; a 32-bit target is a few bytes smaller
(the sample-buffer pointer). **The accelerometer needs no buffer at all** — it
keeps six running averages, which is why it costs nothing beyond the object.
The compass is the only procedure that wants the 3.6 kB.

---

## Building

There is no build system here — these are drop-in sources. Add the `.cpp` files
to your project and put the directory on the include path.

```sh
g++ -std=c++11 -Wall -Wextra -c *.cpp        # compiles clean, no warnings
```

For an embedded target the same applies with your cross-compiler; nothing needs
`-fexceptions`, RTTI or a heap. Include `MathTypes.hpp` (not a HAL-bearing
`common.hpp`) anywhere you need `Vector3f`, `Matrix3f` or `Quaternionf` without
dragging in a vendor SDK.

Hardware integration is your side of the line: feed `addSample()` from whatever
your driver produces, poll `getProgressPercent()` / `isStalled()` to drive
operator prompts, call `calibrate()` once `isReadyToCalibrate()` is true, and
persist the resulting bias/scale/matrix yourself.

---

## Numerical notes

**Results are not bit-reproducible across toolchains.** Binning goes through
`atan2f`, so a 1 ULP difference at a bin boundary puts a sample in a different
bin; coverage is then reached after a different number of samples with a
slightly different set, and the fits differ accordingly. Accuracy is unaffected
(medians moved by under a tenth of a degree across 40 airframes), but **do not
diff results against a reference capture taken on a different libm**. This
applies to the compass only — six-position and level do no binning and are
bit-reproducible. A compass sweep sitting exactly on the coverage minimum can
likewise go either way, complete or stalled, depending on the libm; 2 of 20
borderline seeds were decided differently between implementations. A sweep with
real margin is unaffected.

**Everything is `float`.** No `double` appears anywhere, so this runs sensibly
on a single-precision FPU. The eigen-decomposition is a cyclic Jacobi on a 3×3
symmetric matrix; the linear solve is Gaussian elimination with partial pivoting
on the 9×9 normal equations.

**Cost.** The accelerometer and level paths are O(1) per sample and involve no
matrix work at all. The compass per-sample path is O(1) apart from a rebin (one pass over
the stored samples) and the cached scatter ratio (one pass plus a 3×3
eigen-decomposition, once per accepted sample — accepted samples are capped by
the buffer, so this is bounded). The fit itself is a single 9×9 solve, paid once
at the end of a procedure that has already run for seconds.

---

## Known limitations

1. **Mounting yaw cannot be measured.** Turning about gravity does not change
   what the accelerometer reads. `LevelCalibrator` accepts a `yaw_offset_deg`
   you supply from outside — the nominal angle the board is bolted at, or a
   value trimmed until heading reads true — but it cannot derive one.
2. **The accelerometer fit leaves all cross-axis terms untouched.** Its
   correction is `diag(1/scale)`, so every off-diagonal term passes through.
   That is what keeps it rotation-free and therefore compatible with
   `LevelCalibrator`; separating the genuine cross-axis part needs an ellipsoid
   fit, which this library no longer carries. Plot corrected samples as a sphere
   to see whether it matters on your hardware.
3. **No timeouts.** Every procedure will collect for ever if the operator never
   satisfies its gates. `isStalled()` is the mechanism for noticing; acting on
   it is the caller's job.
4. **`getMatrix()` is identity until the fit succeeds**, as `correct()` is a
   pass-through until then. It returns `diag(1/scale)` — the same correction
   `correct()` applies, in matrix form for callers that store one affine
   correction per sensor.
5. **Compass coverage under a one-sided sweep** is caught by the scatter test at
   fit and readiness time, not by the bin count. Lowering
   `COMPASS_CAL_MIN_SCATTER_RATIO` toward 0.25 lets never-inverted sweeps
   through again.
6. **Single-threaded.** No internal locking; do not call into one object from an
   ISR and a task at the same time.

---

## Author

KAVINDU — extracted from an STM32F7 flight-controller project.
