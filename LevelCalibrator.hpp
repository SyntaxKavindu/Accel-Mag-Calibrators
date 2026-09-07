/*
 * LevelCalibrator.hpp
 *
 *  Created on: Sep 5, 2026
 *      Author: KAVINDU
 */
#ifndef CALIBRATORS_LEVELCALIBRATOR_HPP_
#define CALIBRATORS_LEVELCALIBRATOR_HPP_

// Maths types only -- deliberately NOT common.hpp, which drags in the STM32
// HAL. This class is pure algorithm and must stay portable; see MathTypes.hpp.
#include "MathTypes.hpp"

// Board-mounting rotation.
//
// AccelerometerCalibrator removes bias, per-axis scale and -- in tumble mode
// -- cross-axis sensitivity. What neither of its modes can remove is the
// airframe sitting rotated relative to the sensor: both fits only ever observe
// the MAGNITUDE of gravity, and rotating a sphere leaves the same sphere, so
// the data carries no information about that rotation at all. Measured, a
// 1.16 deg mounting rotation comes out of both procedures at 1.16 deg.
//
// Recovering it needs an outside reference, which is what this provides: rest
// the airframe on a surface known to be level, and whatever direction gravity
// is measured in must in truth be straight down. The rotation taking one onto
// the other is the mounting error.
//
// TWO LIMITS, both inherent rather than implementation choices:
//
//  1. Gravity fixes only two of three axes. This corrects roll and pitch. A
//     board rotated about the VERTICAL axis is invisible here, because turning
//     about gravity does not change what the accelerometer reads. Yaw
//     mounting error needs mechanical alignment or a known-heading reference.
//
//  2. The result is a rotation of the whole BOARD, so it applies to every
//     sensor on it -- accelerometer, gyroscope and magnetometer alike.
//     Applying it to some and not others leaves them disagreeing about which
//     way the airframe points, which is worse than not correcting at all.
//     Calibrator::correctBoardFrame() is deliberately a separate call for
//     exactly this reason; see the note there.

#define LEVEL_CAL_SAMPLES              200
#define LEVEL_CAL_STANDARD_GRAVITY     9.80665f

// Max distance a sample may sit from the running mean before the average is
// discarded and restarted, in sample units (m/s^2 on the Imu path). Same role
// and same default as the six-position motion gate.
#define LEVEL_CAL_MOTION_THRESHOLD     0.5f

// The measured magnitude must look like 1 g, as a fraction of the nominal. A
// larger discrepancy means the accelerometer calibration that ran before this
// is wrong, or the airframe is not at rest, and the rotation derived from it
// would be meaningless.
#define LEVEL_CAL_MAX_MAGNITUDE_ERROR  0.15f

// Largest mounting error this will accept. Past this the likely explanation is
// a surface that is not level, or the wrong axis convention, rather than a
// genuinely crooked sensor -- and silently storing a large rotation is the one
// outcome that would be hard to notice afterwards.
#define LEVEL_CAL_MAX_TILT_DEG         15.0f

// Samples that may arrive without the average growing before isStalled()
// reports it. At a 200 Hz control loop, about 50 seconds.
#define LEVEL_CAL_STALL_LIMIT          10000U

enum class LevelCalStatus : uint8_t {
    IDLE = 0,
    IN_PROGRESS,
    SUCCESS,
    FAILED_NOT_ENOUGH_SAMPLES,
    FAILED_BAD_MAGNITUDE,
    FAILED_EXCESSIVE_TILT
};

enum class LevelSampleResult : uint8_t {
    ACCEPTED,
    ACCEPTED_DONE,
    REJECTED_MOTION,
    REJECTED_DONE,
    REJECTED_NOT_STARTED
};

class LevelCalibrator {
public:
    LevelCalibrator();
    void reset();

    // expected_level_reading: the direction the accelerometer reads when the
    //   airframe is level and upright, as a unit vector in body axes. It reads
    //   specific force, so with +Z up that is (0,0,+1) -- the same convention
    //   AccelPosition::Z_UP assumes. Pass (0,0,-1) for a +Z-down frame.
    // nominal_g: magnitude of 1 g in sample units (9.80665 for m/s^2).
    // yaw_offset_deg: mounting rotation about the vertical axis. Gravity
    //   cannot supply this -- turning about gravity does not change what the
    //   accelerometer reads -- so it has to come from outside: the nominal
    //   angle the board is bolted at, or trimmed until heading reads true.
    //   Left at 0 the fit corrects roll and pitch only, and a mounting yaw
    //   survives and leaks back into roll and pitch as the airframe tilts
    //   (measured: 1 deg of mounting rotation leaves 0.50 deg of yaw and
    //   0.39 deg of roll/pitch error across a +-34 deg envelope, against 1.10
    //   deg uncorrected). Supplying it removes the rest.
    void begin(const Vector3f &expected_level_reading = Vector3f(0.0f, 0.0f, 1.0f),
               float nominal_g = LEVEL_CAL_STANDARD_GRAVITY,
               float motion_threshold = LEVEL_CAL_MOTION_THRESHOLD,
               float yaw_offset_deg = 0.0f);

    // Feed the accelerometer reading AFTER its own calibration has been
    // applied -- this measures what is left over once bias and scale are gone.
    LevelSampleResult addSample(const Vector3f &corrected_accel);

    bool isReadyToCalibrate() const;
    float getProgressPercent() const;
    // True once LEVEL_CAL_STALL_LIMIT samples have arrived without the average
    // growing, i.e. the airframe will not sit still. Advisory, and clears by
    // itself if progress resumes.
    bool isStalled() const { return _samples_since_progress >= LEVEL_CAL_STALL_LIMIT; }

    LevelCalStatus calibrate();
    LevelCalStatus getStatus() const { return _status; }

    // Rotation from measured board axes to true airframe axes. Apply to every
    // sensor on the board, not just the accelerometer.
    Matrix3f getRotation() const { return _rotation; }
    // Mounting error the fit found, in degrees; -1 until it has succeeded.
    float getTiltDeg() const;
    Vector3f correct(const Vector3f &v) const;

    // Smallest rotation taking unit vector `from` onto unit vector `to`.
    // Smallest matters: it turns about the axis perpendicular to both and so
    // introduces no rotation about `to` itself, which is what keeps this from
    // inventing a yaw correction it has no evidence for. False if the two are
    // antiparallel, where no smallest rotation exists.
    static bool rotationBetween(const Vector3f &from, const Vector3f &to, Matrix3f &out);
    // Right-handed rotation about +Z, used to fold in the yaw offset.
    static Matrix3f yawRotation(float radians);

private:
    LevelCalStatus _status;
    Vector3f _target;
    float _nominal_g;
    float _yaw_offset_rad;
    float _motion_threshold_sq;

    Vector3f _sum;
    Vector3f _running_mean;
    uint16_t _count;

    uint16_t _best_count;
    uint32_t _samples_since_progress;

    Matrix3f _rotation;
    float _tilt_deg;
};

#endif /* CALIBRATORS_LEVELCALIBRATOR_HPP_ */
