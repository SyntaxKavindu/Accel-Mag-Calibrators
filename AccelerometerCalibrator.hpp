/*
 * AccelerometerCalibrator.hpp
 *
 *  Created on: Aug 27, 2026
 *      Author: KAVINDU
 *
 * Six-position accelerometer calibration: rest the airframe on each of its six
 * faces in turn, average what the sensor reads, and solve for bias and
 * per-axis scale.
 *
 *     corrected = (raw - bias) / scale
 *
 * which puts the result on a unit sphere -- 1.0 per axis, in whatever units
 * the samples arrived in.
 *
 * What it does NOT fit
 * --------------------
 * Cross-axis terms. The correction is diag(1/scale), so anything off-diagonal
 * passes through untouched: a chip whose axes are not quite perpendicular, or
 * a board sitting rotated in the airframe, is left alone here.
 *
 * That is deliberate rather than a limitation being tolerated. Of the two
 * kinds, only one is this class's to fix. Writing the sensor model as
 * raw = A*g + bias and splitting A = R*S:
 *
 *   S, symmetric -- per-axis scale and cross-axis sensitivity. The scale half
 *     is what this fits. The cross-axis half is small on a modern MEMS part
 *     and needs an ellipsoid fit to separate, which costs a tumble procedure,
 *     a 3.6 kB sample buffer and a minute of the operator's time.
 *   R, a pure rotation -- the board's mounting angle. No fit that only sees
 *     the MAGNITUDE of gravity can recover this at all: rotating a sphere
 *     leaves the same sphere. It needs an outside reference, which is exactly
 *     what LevelCalibrator is -- a known-level surface, and the residual
 *     attitude stored.
 *
 * So the split across the two classes is clean: this one owns bias and scale
 * and produces NO rotation, LevelCalibrator owns the rotation. Nothing here
 * can fight what that stores, which is the property worth protecting.
 */
#ifndef CALIBRATORS_ACCELEROMETERCALIBRATOR_HPP_
#define CALIBRATORS_ACCELEROMETERCALIBRATOR_HPP_

// Maths types only -- deliberately NOT common.hpp, which drags in the STM32
// HAL. This class is pure algorithm and must stay portable; see MathTypes.hpp.
#include "MathTypes.hpp"

#define ACCEL_CAL_SAMPLES_PER_POSITION 100

// The procedure does not notice when it stops making progress: samples that
// fail the motion gate are simply discarded, isReadyToCalibrate() never
// becomes true, and the caller feeds samples for ever. This is how many
// samples may arrive without the progress figure reaching a new high before
// isStalled() starts reporting it.
//
// A position needs SAMPLES_PER_POSITION consecutive good samples and progress
// climbs with each one, so a healthy run never gaps at all. The allowance is
// for an operator still settling the airframe after sending READY.
//
// This is a SAMPLE COUNT, and nothing in this class measures time -- so what it
// means in seconds is set entirely by how fast the caller feeds it. It was
// chosen against a 200 Hz feed, where it is about 50 seconds. Feed it from a
// 1 kHz loop and the same number is 10 seconds, which an operator settling an
// airframe can easily exceed. Decimate the feed, or rescale this.
#define ACCEL_CAL_SIXPOS_STALL_LIMIT   10000U

// Smallest usable nominal 1 g magnitude. Purely a divide-by-zero guard -- the
// value is unit agnostic, so 9.80665 (m/s^2, what ICM42688P::readData outputs)
// and 2048 (raw LSB at ACCEL_FS_16G) are both valid.
#define ACCEL_CAL_MIN_RADIUS           1e-6f

// Standard gravity, matching ICM42688P's raw -> m/s^2 conversion.
#define ACCEL_CAL_STANDARD_GRAVITY     9.80665f

enum class AccelPosition : uint8_t {
    X_UP = 0, X_DOWN = 1, Y_UP = 2, Y_DOWN = 3, Z_UP = 4, Z_DOWN = 5,
    NUM_POSITIONS = 6
};

enum class AccelSampleResult : uint8_t {
    ACCEPTED = 0,            // counted toward the current position
    ACCEPTED_POSITION_DONE,  // that position now has its full average
    REJECTED_MOTION,         // outside the motion gate; the average restarted
    REJECTED_POSITION_DONE,  // this position is already captured
    REJECTED_NOT_STARTED     // no position is being recorded yet
};

// Reported over telemetry as "$CAL,ACCL,FAIL,<n>", so the VALUES are part of
// the wire format. Numbered explicitly for that reason.
//
// Renumbered when the tumble procedure was removed: the four codes only an
// ellipsoid fit could produce (POOR_COVERAGE, SINGULAR_MATRIX,
// DEGENERATE_ELLIPSOID, POOR_FIT) went with it, and the two that remain moved
// down. Nothing unreachable is kept here to preserve old numbers -- a code that
// cannot happen is worse than a renumbering that can be looked up.
enum class AccelCalStatus : uint8_t {
    IDLE = 0,
    IN_PROGRESS = 1,
    SUCCESS = 2,
    FAILED_NOT_ENOUGH_SAMPLES = 3,   // calibrate() before all six were captured
    FAILED_BAD_GEOMETRY = 4          // a dead axis, or a face captured inverted
};

class AccelerometerCalibrator {
public:
    AccelerometerCalibrator();
    void reset();

    // motion_threshold: max distance a sample may sit from the running mean of
    // the current position before the average is discarded and restarted. In
    // the same units as the samples -- Imu/ICM42688P deliver m/s^2, so the
    // default is a few times a typical at-rest noise floor in m/s^2, NOT the
    // raw-LSB figure this originally carried.
    //
    // A position needs ACCEL_CAL_SAMPLES_PER_POSITION consecutive samples
    // inside this radius, and any single excursion restarts the count, so the
    // vibration it tolerates is well below the threshold itself. Measured with
    // the 0.5 default: completes reliably up to 0.12 m/s^2 RMS per axis, 1 run
    // in 20 at 0.20, never at 0.30 -- and with no timeout behind it, failing
    // to complete means waiting for ever rather than being told.
    void beginSixPosition(float motion_threshold = 0.5f);
    void startPosition(AccelPosition pos);
    bool isPositionDone(AccelPosition pos) const;
    AccelPosition getCurrentPosition() const { return _current_pos; }
    bool allPositionsComplete() const;

    AccelSampleResult addSample(float x, float y, float z);
    AccelSampleResult addSample(const Vector3f &s) { return addSample(s.x, s.y, s.z); }

    bool isReadyToCalibrate() const;
    float getProgressPercent() const;

    // True once ACCEL_CAL_SIXPOS_STALL_LIMIT samples have arrived without
    // progress reaching a new high -- almost always vibration defeating the
    // motion gate. Advisory: collection continues and the flag clears by
    // itself if progress resumes, so it is for the caller to decide whether to
    // abort and say so.
    bool isStalled() const;

    AccelCalStatus calibrate();
    AccelCalStatus getStatus() const { return _status; }

    Vector3f correct(const Vector3f &raw) const;
    Vector3f getBias() const;
    Vector3f getScale() const;
    Matrix3f getMatrix() const;

private:
    AccelCalStatus _status;

    AccelPosition _current_pos;
    Vector3f _pos_sum[(size_t)AccelPosition::NUM_POSITIONS];
    uint16_t _pos_count[(size_t)AccelPosition::NUM_POSITIONS];
    bool _pos_done[(size_t)AccelPosition::NUM_POSITIONS];
    Vector3f _pos_avg[(size_t)AccelPosition::NUM_POSITIONS];

    bool _sp_collecting;
    Vector3f _sp_running_mean;
    float _sp_motion_threshold_sq;
    Vector3f _bias_sixpos;
    Vector3f _scale_sixpos;

    // Stall detection. See isStalled().
    float _best_progress;
    uint32_t _samples_since_progress;
    void noteProgress();
};

#endif /* CALIBRATORS_ACCELEROMETERCALIBRATOR_HPP_ */
