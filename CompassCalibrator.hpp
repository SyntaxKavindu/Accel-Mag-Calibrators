/*
 * CompassCalibrator.hpp
 *
 *  Created on: Aug 27, 2026
 *      Author: KAVINDU
 */
#ifndef CALIBRATORS_COMPASSCALIBRATOR_HPP_
#define CALIBRATORS_COMPASSCALIBRATOR_HPP_

// Maths types only -- deliberately NOT common.hpp, which drags in the STM32
// HAL. This class is pure algorithm and must stay portable; see MathTypes.hpp.
#include "MathTypes.hpp"

// Recommended sample-buffer capacity -- the caller now supplies the storage, so
// this sizes it rather than declaring an array. See begin().
#define COMPASS_CAL_MAX_SAMPLES    300
#define COMPASS_CAL_MIN_SAMPLES    150
#define COMPASS_CAL_NUM_BINS       64
#define COMPASS_CAL_MIN_BINS       45  // Requires >70% spherical coverage
#define COMPASS_CAL_MAX_PER_BIN    5   // Uniform point distribution cap

// Coverage is judged from the centre of the collected cloud rather than from
// the origin, so that a hard-iron offset larger than the field itself cannot
// collapse every sample into a handful of bins. That centre is only knowable
// once samples are in, so the bin table has to be rebuilt as it settles --
// see rebin(). Rebinning is triggered by the centre having moved rather than
// by a sample count, because it is only worth redoing while the estimate is
// still shifting; this is how far it must move to trigger one, as a fraction
// of the nominal radius so that it carries over to whatever units the sensor
// reports in.
#define COMPASS_CAL_REBIN_MOVE_FRACTION  0.10f

// A fit can satisfy every structural check and still not describe the data:
// readings with no ellipsoidal structure at all (a failing sensor, or an
// interference source that moved during the sweep) still yield a solution
// that inverts cleanly and is positive definite. The only thing that catches
// that is asking how well the result actually reproduces the samples. This is
// the largest RMS deviation from the nominal magnitude, as a fraction of it,
// that a fit may leave behind and still be accepted. Generous on purpose --
// it is here to reject nonsense, not to grade a good calibration.
#define COMPASS_CAL_MAX_FIT_RESIDUAL 0.15f

// Smallest scatter-matrix eigenvalue ratio a sweep may have and still be
// accepted -- the test for "did the operator actually turn the airframe over?".
//
// Measured on synthetic sweeps: a full sphere gives 1.00, a full sphere with
// heavy hard AND soft iron 0.70, a full sphere with very uneven dwell 0.81,
// two thirds of a sphere 0.65 -- against 0.38 for a hemisphere plus a little,
// 0.25 for an exact hemisphere (the analytic value: variance r^2/12 across the
// cap axis versus r^2/3 along it) and 0.22 for a hemisphere with iron. 0.5 sits
// in the gap with room on both sides.
//
// Raising this rejects lazier sweeps; lowering it toward 0.25 lets a
// never-inverted sweep through again. Do not go below ~0.3: that is the region
// where a fit is exact on the samples collected and materially wrong on the
// orientations that were not.
#define COMPASS_CAL_MIN_SCATTER_RATIO 0.5f

// Collection makes no progress at all if the airframe stops being moved, and
// nothing in the procedure notices: isReadyToCalibrate() simply never becomes
// true and the caller feeds samples for ever. This is how many samples may
// arrive without the progress figure reaching a new high before isStalled()
// starts reporting it. Measured over successful runs the longest such gap was
// 1154 samples, so this leaves roughly 5x headroom; at a 200 Hz control loop
// it is about 30 seconds of no progress.
//
// Note that a sweep sitting exactly on the coverage minimum can go either way
// -- complete, or be reported stalled -- depending on the libm behind atan2f,
// because a 1 ULP difference there moves a sample across a bin boundary.
// Measured between two implementations that decided 2 of 20 borderline seeds
// differently. Both outcomes are correct for input that marginal, and neither
// is a silent hang; a sweep with real margin is unaffected.
#define COMPASS_CAL_STALL_LIMIT    6000U

// Smallest usable nominal field magnitude. Only guards against a zero or
// negative radius -- the value itself is sensor-unit agnostic, so a LIS3MDL
// reporting ~0.5 Gauss is just as valid as a raw-LSB magnitude of ~500.
#define COMPASS_CAL_MIN_RADIUS     1e-6f

enum class CalStatus : uint8_t {
    IDLE = 0,
    COLLECTING,
    READY_TO_FIT,
    SUCCESS,
    FAILED_NOT_ENOUGH_SAMPLES,
    FAILED_POOR_COVERAGE,
    FAILED_SINGULAR_MATRIX,
    FAILED_DEGENERATE_ELLIPSOID,
    FAILED_POOR_FIT
};

enum class SampleResult : uint8_t {
    ACCEPTED,
    REJECTED_TOO_CLOSE,
    REJECTED_BUFFER_FULL,
    REJECTED_NOT_COLLECTING
};

class CompassCalibrator {
public:
    CompassCalibrator();

    // sample_buffer / capacity: where the samples are collected. Caller-owned,
    // must outlive the procedure.
    //
    // A parameter rather than a member array so the ~3.6 kB it costs is the
    // integrator's to place. A procedure that runs for a few seconds on the
    // ground and never again has no business holding that permanently in .bss
    // on a part where it is a meaningful fraction of the RAM -- and on a part
    // where it is not, passing a static array costs nothing.
    //
    // Returns false -- and does NOT start -- if the buffer is null or smaller
    // than COMPASS_CAL_MIN_SAMPLES, which could never complete.
    bool begin(float nominal_field_magnitude, Vector3f *sample_buffer,
               uint16_t capacity);
    SampleResult addSample(float x, float y, float z);
    SampleResult addSample(const Vector3f &s) { return addSample(s.x, s.y, s.z); }

    float getProgressPercent() const;
    bool isReadyToCalibrate() const;
    // True once COMPASS_CAL_STALL_LIMIT samples have arrived without progress
    // reaching a new high -- the airframe is being held still, or the field is
    // too disturbed to bin anything new. Advisory: collection continues, and
    // the flag clears by itself if progress resumes. It is for the caller to
    // decide whether to abort and say so.
    bool isStalled() const { return _samples_since_progress >= COMPASS_CAL_STALL_LIMIT; }
    CalStatus calibrate();

    CalStatus getStatus() const { return _status; }
    uint16_t getSampleCount() const { return _sample_count; }
    uint8_t getFilledBinsCount() const { return _filled_bins; }

    // How well the last calibrate() reproduced its own samples: RMS deviation
    // from the nominal magnitude, as a fraction of it. -1 until one has run.
    //
    // This is the number that decides FAILED_POOR_FIT, and without it a
    // rejection is unactionable -- 0.16 is a sweep that nearly worked and 0.90
    // is a sensor or an environment that carries no usable field. The
    // difference matters enormously to whoever has to decide what to do next,
    // and it costs one float to report.
    float getLastFitResidual() const { return _last_residual; }

    // Shape of the sample cloud, 1.0 for a full sphere and 0.25 for a
    // hemisphere. See scatterAnisotropy().
    float getScatterRatio() const { return _scatter_ratio; }

    Vector3f correct(const Vector3f &raw) const;
    Vector3f getOffset() const { return _offset; }
    float getNominalRadius() const { return _nominal_radius; }
    Matrix3f getSoftIronMatrix() const { return _softiron; }

    void reset();

private:
    // Caller-owned; see begin(). Null whenever no procedure is running.
    Vector3f *_samples;
    uint16_t _capacity;
    uint16_t _sample_count;

    // Running sum and count of every sample ever accepted, used only to derive
    // the cloud centre. Deliberately NOT rolled back when rebin() drops a
    // sample: if dropping moved the centre, the centre would reassign the bins,
    // which would drop a different set, and the two would chase each other
    // instead of converging. Keeping this monotone breaks that loop.
    Vector3f _sum;
    uint32_t _accepted_total;
    // Centroid as it stood at the last rebin; rebin() runs again once the
    // current one has moved away from it far enough to change bin assignments.
    Vector3f _rebin_centre;

    // Stall detection. _best_progress is the highest getProgressPercent() has
    // reached, not its current value: a rebin can drop samples and move the
    // figure backwards, and a temporary dip like that is not progress.
    float _best_progress;
    uint32_t _samples_since_progress;
    
    // Spherical binning storage
    uint8_t _bin_count[COMPASS_CAL_NUM_BINS];
    uint8_t _filled_bins;

    // Cached scatterAnisotropy(), refreshed once per accepted sample. Readiness
    // and progress both read it, and both are queried far more often than the
    // cloud actually changes.
    float _scatter_ratio;

    // Residual from the last calibrate(); see getLastFitResidual().
    float _last_residual;

    float _nominal_radius;
    CalStatus _status;

    Vector3f _offset;
    Matrix3f _softiron;

    // Approximate centre of the collected cloud -- the origin until the first
    // sample lands. A first-order stand-in for the hard-iron offset, which is
    // all that the binning and the fit pre-centring need from it.
    // Updates the stall counter from the current progress figure.
    void noteProgress();
    Vector3f centroid() const;
    // RMS deviation of the corrected samples from the nominal magnitude,
    // relative to it. Computed against the candidate gains, before they are
    // committed to _offset/_softiron.
    float fitResidual(const Vector3f &offset, const Matrix3f &softiron) const;
    // Recompute the whole bin table about the current centroid, dropping the
    // samples the new binning makes redundant.
    void rebin();

    // Smallest/largest eigenvalue of the sample scatter matrix about the
    // centroid: 1 for a cloud spread evenly over a sphere, 0.25 for one
    // covering a hemisphere, ->0 for a narrow cap. See the note at the
    // definition for why this is the check that works.
    float scatterAnisotropy() const;

    static uint8_t getBinIndex(const Vector3f &s);

    static bool solve9x9(float A[9][9], float b[9], float x[9]);
    static bool invert3x3(const Matrix3f &in, Matrix3f &out);
    static void eigenSymmetric3x3(Matrix3f m, float eigval[3], Matrix3f &eigvec);
};

#endif /* CALIBRATORS_COMPASSCALIBRATOR_HPP_ */
