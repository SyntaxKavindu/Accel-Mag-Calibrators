/*
 * AccelerometerCalibrator.hpp
 *
 *  Created on: Aug 27, 2026
 *      Author: KAVINDU
 */
#ifndef CALIBRATORS_ACCELEROMETERCALIBRATOR_HPP_
#define CALIBRATORS_ACCELEROMETERCALIBRATOR_HPP_

// Maths types only -- deliberately NOT common.hpp, which drags in the STM32
// HAL. This class is pure algorithm and must stay portable; see MathTypes.hpp.
#include "MathTypes.hpp"

#define ACCEL_CAL_SAMPLES_PER_POSITION 100
// Recommended sample-buffer capacity for a tumble -- the caller now supplies
// the storage, so this sizes it rather than declaring an array. See
// beginTumble(). Six-position needs no buffer at all: it keeps six running
// averages, which is why that mode costs nothing here.
#define ACCEL_CAL_TUMBLE_MAX_SAMPLES   300
#define ACCEL_CAL_TUMBLE_MIN_SAMPLES   150
#define ACCEL_CAL_STILLNESS_WINDOW     12
#define ACCEL_CAL_TUMBLE_NUM_BINS      64
#define ACCEL_CAL_TUMBLE_MIN_BINS      45
#define ACCEL_CAL_TUMBLE_MAX_PER_BIN   5

// A fit can satisfy every structural check and still not describe the data:
// readings with no ellipsoidal structure at all still yield a solution that
// inverts cleanly and is positive definite. The only thing that catches that
// is asking how well the result reproduces the samples. This is the largest
// RMS deviation from 1 g, as a fraction of it, that a tumble fit may leave
// behind and still be accepted. Generous on purpose -- it is here to reject
// nonsense, not to grade a good calibration.
#define ACCEL_CAL_MAX_FIT_RESIDUAL     0.15f

// Note that binning goes through atan2f, so which samples a tumble admits is
// only as reproducible as the libm behind it: a 1 ULP difference at a bin
// boundary puts a sample in a different bin, and coverage is then reached
// after a different number of samples with a slightly different set. Measured
// between two implementations, the same airframe reached coverage at 24033
// samples on one and 26647 on the other, and the fits differed accordingly --
// both correct, neither bit-identical. Accuracy is unaffected (the medians
// moved by under a tenth of a degree across 40 airframes); results simply are
// not bit-reproducible across toolchains, so do not diff them against a
// reference capture taken on a different libm.
//
// Neither procedure notices when it stops making progress: samples that fail
// the motion or stillness gate are simply discarded, isReadyToCalibrate()
// never becomes true, and the caller feeds samples for ever. These are how
// many samples may arrive without the progress figure reaching a new high
// before isStalled() starts reporting it.
//
// The two differ by a lot because of how the gates work. Six-position needs
// SAMPLES_PER_POSITION consecutive good samples and progress climbs with each
// one, so a healthy run never gaps at all; the allowance here is for an
// operator still settling the airframe after sending READY. Tumble only
// accepts a sample once a whole window is still AND lands in an unfilled bin,
// so long gaps are normal -- the longest measured in a successful run was
// 10151 samples, hence the much larger figure. At a 200 Hz control loop these
// are roughly 50 and 200 seconds of no progress.
#define ACCEL_CAL_SIXPOS_STALL_LIMIT   10000U
#define ACCEL_CAL_TUMBLE_STALL_LIMIT   40000U

// Smallest usable nominal 1 g magnitude. Purely a divide-by-zero guard -- the
// value is unit agnostic, so 9.80665 (m/s^2, what ICM42688P::readData outputs)
// and 2048 (raw LSB at ACCEL_FS_16G) are both valid.
#define ACCEL_CAL_MIN_RADIUS           1e-6f

// Standard gravity, matching ICM42688P's raw -> m/s^2 conversion.
#define ACCEL_CAL_STANDARD_GRAVITY     9.80665f

enum class AccelCalMode : uint8_t {
    NONE,
    SIX_POSITION,
    TUMBLE
};

enum class AccelPosition : uint8_t {
    X_UP=0, X_DOWN=1, Y_UP=2, Y_DOWN=3, Z_UP=4, Z_DOWN=5, NUM_POSITIONS=6
};

enum class AccelSampleResult : uint8_t {
    ACCEPTED,
    ACCEPTED_POSITION_DONE,
    REJECTED_MOTION,
    REJECTED_TOO_CLOSE,
    REJECTED_POSITION_DONE,
    REJECTED_BUFFER_FULL,
    REJECTED_NOT_STARTED
};

enum class AccelCalStatus : uint8_t {
    IDLE,
    IN_PROGRESS,
    SUCCESS,
    FAILED_NOT_ENOUGH_SAMPLES,
    FAILED_POOR_COVERAGE,
    FAILED_BAD_GEOMETRY,
    FAILED_SINGULAR_MATRIX,
    FAILED_DEGENERATE_ELLIPSOID,
    FAILED_POOR_FIT
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
    //
    // Note also that six-position fits bias and per-axis scale only -- its
    // correction matrix is diag(1/scale), so every off-diagonal term passes
    // through untouched. Measured, the tilt error it leaves tracks the total
    // misalignment essentially 1:1 (0.009 deg with none, 1.15 deg at 1.15 deg
    // of it). getMaxMisalignmentDeg() reports what is being left behind.
    //
    // Which part of that a tumble can recover depends on the kind, because
    // raw = A*g + bias splits as A = R*S:
    //   S, symmetric (per-axis scale and cross-axis sensitivity, i.e. the
    //     chip's own axes not being quite perpendicular) -- a tumble removes
    //     this: 1.15 deg -> 0.04 deg measured.
    //   R, a pure rotation (the board sitting rotated in the airframe) --
    //     NEITHER procedure can remove this, and no amount of tumbling will.
    //     Both fits only ever see the magnitude of gravity, and rotating a
    //     sphere leaves the same sphere, so the data carries no information
    //     about R at all. Measured, 1.16 deg in and 1.16 deg out of both.
    // With both present a tumble removes the S half only (2.03 -> 1.23 deg).
    // Correcting R needs an outside reference -- levelling the airframe on a
    // known-flat surface and storing the residual attitude -- which is a
    // separate step from anything in this class.
    void beginSixPosition(float motion_threshold = 0.5f);
    void startPosition(AccelPosition pos);
    bool isPositionDone(AccelPosition pos) const;
    AccelPosition getCurrentPosition() const { return _current_pos; }
    bool allPositionsComplete() const;
    float getMaxMisalignmentDeg() const;

    // nominal_g: magnitude the sensor reads under 1 g, in sample units
    //            (9.80665 for m/s^2 input, ~2048 for raw LSB at FS_16G).
    // stillness_threshold: max RMS deviation across the stillness window, in
    //            sample units. Again m/s^2 for the Imu path.
    //
    //            This is compared against the deviation of the whole 3-vector,
    //            not of one axis, so the per-axis noise it tolerates is
    //            threshold/sqrt(3) -- the default 0.2 accepts about 0.115
    //            m/s^2 RMS per axis. Measured, collection completes reliably
    //            up to 0.12 and never completes at 0.15 or above, and there is
    //            no timeout behind it: past that point the procedure runs for
    //            ever with progress stuck. A part at rest is far below this
    //            (ICM42688P is ~0.006 m/s^2), but a bench that shakes or props
    //            turning are not. If calibration stops advancing, that is what
    //            it means.
    //
    // sample_buffer / capacity: where the tumble samples are collected. The
    // buffer is NOT owned by this class and must outlive the procedure.
    //
    // It is a caller parameter rather than a member array so that the ~3.6 kB
    // it costs is the integrator's to place. Two of these classes each holding
    // their own array put 7.2 kB permanently in .bss for procedures that run
    // for a few seconds on the ground and never again -- survivable on a part
    // with 320 kB, fatal on one with 20 kB. Passing it in lets the memory be
    // shared with the compass calibrator (they are mutually exclusive), live in
    // a scratch/CCM region, or be a stack buffer in the calibration routine.
    //
    // Returns false -- and does NOT start -- if the buffer is null or smaller
    // than ACCEL_CAL_TUMBLE_MIN_SAMPLES, which could never complete.
    bool beginTumble(float nominal_g, float stillness_threshold,
                     Vector3f *sample_buffer, uint16_t capacity);

    AccelSampleResult addSample(float x, float y, float z);
    AccelSampleResult addSample(const Vector3f &s) { return addSample(s.x, s.y, s.z); }

    bool isReadyToCalibrate() const;
    float getProgressPercent() const;
    // True once the mode's stall limit of samples has arrived without progress
    // reaching a new high -- almost always vibration defeating the motion or
    // stillness gate. Advisory: collection continues and the flag clears by
    // itself if progress resumes, so it is for the caller to decide whether to
    // abort and say so.
    bool isStalled() const;

    AccelCalStatus calibrate();
    AccelCalStatus getStatus() const { return _status; }
    AccelCalMode getMode() const { return _mode; }

    Vector3f correct(const Vector3f &raw) const;
    Vector3f getBias() const;
    // Sample-unit magnitude of 1 g the tumble fit was normalised against.
    // Needed to fold the tumble matrix into a single correction matrix.
    float getNominalRadius() const { return _nominal_radius; }
    Vector3f getScale() const;
    Matrix3f getMatrix() const;

private:
    AccelCalMode _mode;
    AccelCalStatus _status;

    // Six-position state
    Vector3f _pos_sum[(size_t)AccelPosition::NUM_POSITIONS];
    uint16_t _pos_count[(size_t)AccelPosition::NUM_POSITIONS];
    bool _pos_done[(size_t)AccelPosition::NUM_POSITIONS];
    Vector3f _pos_avg[(size_t)AccelPosition::NUM_POSITIONS];
    AccelPosition _current_pos;
    bool _sp_collecting;
    Vector3f _sp_running_mean;
    float _sp_motion_threshold_sq;
    Vector3f _bias_sixpos;
    Vector3f _scale_sixpos;

    AccelSampleResult addSampleSixPosition(const Vector3f &s);
    AccelCalStatus calibrateSixPosition();

    // Tumble state
    // Stall detection, shared by both modes. _best_progress is the highest
    // getProgressPercent() has reached, not its current value: a rejected
    // six-position sample resets that position's count and moves the figure
    // backwards, and that dip is not progress.
    float _best_progress;
    uint32_t _samples_since_progress;
    void noteProgress();

    // Caller-owned; see beginTumble(). Null whenever no tumble is running, so a
    // stale pointer from a finished procedure can never be dereferenced.
    Vector3f *_samples;
    uint16_t _capacity;
    uint16_t _sample_count;
    uint8_t _bin_count[ACCEL_CAL_TUMBLE_NUM_BINS];
    uint8_t _filled_bins;
    float _nominal_radius;

    Vector3f _window[ACCEL_CAL_STILLNESS_WINDOW];
    uint8_t _window_count;
    uint8_t _window_head;
    float _stillness_threshold_sq;
    Vector3f _offset_tumble;
    Matrix3f _matrix_tumble;

    AccelSampleResult addSampleTumble(const Vector3f &s);
    AccelCalStatus calibrateTumble();
    bool isWindowStill(Vector3f &out_mean) const;
    // Mean of the collected samples -- a first-order stand-in for the bias,
    // used to pre-centre the ellipsoid fit. Unlike the compass this is not
    // used for binning: see addSampleTumble().
    Vector3f centroid() const;
    // RMS deviation of the corrected tumble samples from 1 g, relative to it.
    // Computed against the candidate gains, before they are committed.
    float fitResidual(const Vector3f &offset, const Matrix3f &matrix) const;
    static uint8_t getBinIndex(const Vector3f &s);

    static bool solve9x9(float A[9][9], float b[9], float x[9]);
    static bool invert3x3(const Matrix3f &in, Matrix3f &out);
    static void eigenSymmetric3x3(Matrix3f m, float eigval[3], Matrix3f &eigvec);
};

#endif /* CALIBRATORS_ACCELEROMETERCALIBRATOR_HPP_ */
