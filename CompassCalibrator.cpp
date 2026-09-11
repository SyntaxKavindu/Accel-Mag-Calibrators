/*
 * CompassCalibrator.cpp
 *
 *  Created on: Aug 27, 2026
 *      Author: KAVINDU
 */
#include "CompassCalibrator.hpp"

CompassCalibrator::CompassCalibrator() {
    reset();
}

void CompassCalibrator::reset() {
    _samples = nullptr;
    _capacity = 0;
    _sample_count = 0;
    _sum = Vector3f();
    _accepted_total = 0;
    _rebin_centre = Vector3f();
    _best_progress = -1.0f;
    _samples_since_progress = 0;
    _filled_bins = 0;
    memset(_bin_count, 0, sizeof(_bin_count));
    _scatter_ratio = 0.0f;
    _last_residual = -1.0f;
    _nominal_radius = 500.0f;
    _status = CalStatus::IDLE;
    _offset = Vector3f();
    _softiron = Matrix3f::identity();
}

bool CompassCalibrator::begin(float nominal_field_magnitude,
                              Vector3f *sample_buffer, uint16_t capacity) {
    reset();

    // Too small to ever reach COMPASS_CAL_MIN_SAMPLES: refuse up front rather
    // than collect for ever and fail the fit. Note rebin() can DROP samples, so
    // a capacity merely equal to the minimum is legitimate -- coverage, not the
    // buffer, is what finishes the run.
    if (sample_buffer == nullptr || capacity < COMPASS_CAL_MIN_SAMPLES) {
        return false;
    }
    _samples  = sample_buffer;
    _capacity = capacity;
    // The nominal magnitude is both the normalisation used by the fit and the
    // output scale of correct(). It must be strictly positive; anything else
    // would make the fit divide by ~zero, so fall back to the default rather
    // than producing silently garbage gains.
    if (nominal_field_magnitude > COMPASS_CAL_MIN_RADIUS) {
        _nominal_radius = nominal_field_magnitude;
    }
    _status = CalStatus::COLLECTING;
    return true;
}

Vector3f CompassCalibrator::centroid() const {
    if (_accepted_total == 0) return Vector3f();
    return _sum / (float)_accepted_total;
}

// The bin a sample falls in depends on the centroid, and the centroid moves as
// the cloud fills out, so the early bins were assigned about a centre that was
// still wrong. Re-deriving the whole table from the stored samples is what
// keeps _filled_bins an honest coverage measure rather than an optimistic one.
//
// Samples the new binning makes redundant are dropped, so _sample_count can
// step backwards during collection. That is deliberate. Readings arrive once
// per control loop while the airframe is turned by hand, so most of them
// repeat an orientation already held; leaving the redundant ones in place
// fills all 300 slots with duplicates and coverage stalls one bin short of
// the minimum, for ever (measured: stuck at 44 of 64 at every hard-iron
// level). Dropping them frees the slots for orientations still missing.
//
// The centroid is deliberately NOT recomputed from the surviving samples --
// see _sum. If dropping moved the centre, the centre would reassign the bins,
// which would drop a different set, and the two would chase each other rather
// than converge (measured: coverage cycling between 51% and 98%, never
// finishing).
void CompassCalibrator::rebin() {
    const Vector3f c = centroid();

    memset(_bin_count, 0, sizeof(_bin_count));
    _filled_bins = 0;

    uint16_t keep = 0;
    for (uint16_t i = 0; i < _sample_count; i++) {
        const Vector3f s = _samples[i];
        const uint8_t bin = getBinIndex(s - c);
        if (_bin_count[bin] >= COMPASS_CAL_MAX_PER_BIN) continue;
        if (_bin_count[bin] == 0) _filled_bins++;
        _bin_count[bin]++;
        _samples[keep++] = s;
    }
    _sample_count = keep;
}

// Smallest/largest eigenvalue of the sample scatter (second-moment) matrix
// taken about the centroid.
//
// This is the check that finally separates "the airframe was turned all the way
// over" from "the airframe was waved around upright", and it works for one
// reason: a scatter matrix is invariant to TRANSLATION. Hard iron is a
// translation, so it cannot influence this number at all -- which is precisely
// what defeated every direction-based coverage test tried before it. A sweep
// covering the whole sphere is isotropic and returns ~1; one covering a
// hemisphere is squashed along the cap axis and returns ~0.25 (analytically
// r^2/12 over r^2/3); a narrow cap tends to 0.
//
// Soft iron does affect it -- it scales the cloud along the ellipsoid axes --
// but only by the square of the axis ratio, which for a realistic +/-15% part
// still leaves ~0.7, far above the threshold.
//
// Cost: one pass over the samples plus one 3x3 symmetric eigen-decomposition,
// paid once at the end of a calibration that has already run for seconds.
float CompassCalibrator::scatterAnisotropy() const {
    if (_samples == nullptr || _sample_count < 4) return 0.0f;

    const Vector3f c = centroid();

    Matrix3f M = Matrix3f();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) M.m[i][j] = 0.0f;
    }
    for (uint16_t k = 0; k < _sample_count; k++) {
        const Vector3f d = _samples[k] - c;
        const float v[3] = { d.x, d.y, d.z };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) M.m[i][j] += v[i] * v[j];
        }
    }
    const float inv_n = 1.0f / (float)_sample_count;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) M.m[i][j] *= inv_n;
    }

    float eig[3];
    Matrix3f vec;
    eigenSymmetric3x3(M, eig, vec);

    float mn = eig[0], mx = eig[0];
    for (int i = 1; i < 3; i++) {
        if (eig[i] < mn) mn = eig[i];
        if (eig[i] > mx) mx = eig[i];
    }
    if (!(mx > 0.0f)) return 0.0f;
    if (mn < 0.0f) mn = 0.0f;   // round-off on a near-degenerate cloud
    return mn / mx;
}

// O(1) Spherical Bin Index mapping. The argument is a direction from the
// centre of the cloud, not a raw reading -- see rebin().
uint8_t CompassCalibrator::getBinIndex(const Vector3f &s) {
    float mag = sqrtf(s.x * s.x + s.y * s.y + s.z * s.z);
    if (mag < 1e-6f) return 0;
    
    // Elevation index (8 bands along Z-axis)
    float zn = s.z / mag;
    int elev = (int)((zn + 1.0f) * 4.0f);
    if (elev < 0) elev = 0; else if (elev > 7) elev = 7;

    // Azimuth index (8 angular slices around XY plane)
    float az_rad = atan2f(s.y, s.x);
    int az = (int)((az_rad + 3.1415926535f) * (8.0f / 6.283185307f));
    if (az < 0) az = 0; else if (az > 7) az = 7;

    return (uint8_t)(elev * 8 + az);
}

SampleResult CompassCalibrator::addSample(float x, float y, float z) {
    // Buffer state is checked before the mode check: once the last slot is
    // taken addSample() flips the status to READY_TO_FIT, so testing the mode
    // first would report every later sample as REJECTED_NOT_COLLECTING and the
    // caller could never tell "buffer full" from "never started".
    if (_samples == nullptr || _sample_count >= _capacity) {
        if (_status == CalStatus::COLLECTING) _status = CalStatus::READY_TO_FIT;
        return SampleResult::REJECTED_BUFFER_FULL;
    }
    if (_status != CalStatus::COLLECTING) {
        return SampleResult::REJECTED_NOT_COLLECTING;
    }

    Vector3f s(x, y, z);
    // Binned by direction from the centre of the cloud. Binning about the
    // origin instead made coverage collapse once the hard iron exceeded the
    // field magnitude -- every reading then points into the same cone, so the
    // bins saturate at ~35 of 64 and COMPASS_CAL_MIN_BINS can never be met.
    //
    // KNOWN TRADE-OFF: this also weakens the check on a one-sided sweep. An
    // operator who never inverts the airframe covers barely half the sphere,
    // and the centroid then sits inside that cap rather than at the sphere's
    // centre, so directions measured from it over-disperse and the bin count
    // still reaches the minimum -- where origin-relative binning would have
    // rejected it. Such a fit is exact on noiseless data but extrapolates
    // poorly: measured at 6% sensor noise, 3% error in the centre and 5% in
    // corrected magnitude on the orientations never visited, while looking
    // fine on the ones that were. fitResidual() cannot see this; it only sees
    // the samples that exist. Re-running the coverage test against the fitted
    // centre was tried and does not separate the cases -- a good full-sphere
    // sweep leaves 35 bins occupied there and a 70%-of-sphere sweep leaves 36.
    //
    // RESOLVED at fit time instead, by scatterAnisotropy(): the SHAPE of the
    // cloud gives the sweep away where its bin count does not. See the note
    // there. This binning stays exactly as it is -- it is what makes coverage
    // reachable at all under heavy hard iron -- and the one-sided case is now
    // caught in calibrate() rather than left to operator discipline.
    uint8_t bin = getBinIndex(s - centroid());

    // Reject if bin has reached maximum point density cap
    if (_bin_count[bin] >= COMPASS_CAL_MAX_PER_BIN) {
        noteProgress();
        return SampleResult::REJECTED_TOO_CLOSE;
    }

    if (_bin_count[bin] == 0) {
        _filled_bins++;
    }
    _bin_count[bin]++;

    _samples[_sample_count] = s;
    _sample_count++;
    _sum = _sum + s;
    _accepted_total++;

    // Only rebin while the centre estimate is still moving. Once it settles
    // the bin assignments stop changing, no more samples are dropped, and
    // coverage climbs at full rate.
    const Vector3f c = centroid();
    const Vector3f moved = c - _rebin_centre;
    const float move_sq = moved.x * moved.x + moved.y * moved.y
            + moved.z * moved.z;
    const float move_tol = COMPASS_CAL_REBIN_MOVE_FRACTION * (_nominal_radius);
    if (move_sq > move_tol * move_tol) {
        _rebin_centre = c;
        rebin();
    }

    // Recomputed here -- once per ACCEPTED sample, after any rebin -- rather
    // than inside getProgressPercent()/isReadyToCalibrate(), which the control
    // loop calls on every sample whether it was accepted or not. Accepted
    // samples are capped by the buffer, so this is bounded work.
    _scatter_ratio = scatterAnisotropy();

    if (_sample_count >= _capacity) {
        _status = CalStatus::READY_TO_FIT;
    }

    noteProgress();
    return SampleResult::ACCEPTED;
}

// Called on every sample the caller offers, accepted or not: a run of samples
// that are all rejected is exactly the situation this has to notice.
void CompassCalibrator::noteProgress() {
    const float p = getProgressPercent();
    if (p > _best_progress) {
        _best_progress = p;
        _samples_since_progress = 0;
    } else if (_samples_since_progress < COMPASS_CAL_STALL_LIMIT) {
        _samples_since_progress++;
    }
}

float CompassCalibrator::getProgressPercent() const {
    float sample_frac = (float)_sample_count / (float)COMPASS_CAL_MIN_SAMPLES;
    if (sample_frac > 1.0f) sample_frac = 1.0f;

    float bin_frac = (float)_filled_bins / (float)COMPASS_CAL_MIN_BINS;
    if (bin_frac > 1.0f) bin_frac = 1.0f;

    // The WORST of the three requirements, not their average: all three must be
    // met, so the smallest is the honest answer to "how far along is this?".
    // The scatter term is what makes the figure stall for an operator who keeps
    // waving the airframe about upright, and a stalled figure is exactly what
    // isStalled() reports.
    float scatter_frac = _scatter_ratio / COMPASS_CAL_MIN_SCATTER_RATIO;
    if (scatter_frac > 1.0f) scatter_frac = 1.0f;

    float combined = sample_frac < bin_frac ? sample_frac : bin_frac;
    if (scatter_frac < combined) combined = scatter_frac;
    return combined * 100.0f;
}

// The scatter requirement is part of READINESS, not just a check at fit time.
//
// Failing it only in calibrate() would be correct but hostile: the operator is
// told "ready", the fit is then rejected, and the whole sweep has to start again
// with no indication of what was wrong. Requiring it here means collection
// simply continues until the airframe has actually been turned over -- and the
// existing stall machinery already carries the right advice for a sweep that
// never gets there ("KEEP TURNING THE AIRFRAME - ALL SIDES, INCLUDING
// INVERTED").
//
// calibrate() re-checks it regardless: this class does not assume its caller
// consulted isReadyToCalibrate() first.
bool CompassCalibrator::isReadyToCalibrate() const {
    return (_sample_count >= COMPASS_CAL_MIN_SAMPLES)
        && (_filled_bins >= COMPASS_CAL_MIN_BINS)
        && (_scatter_ratio >= COMPASS_CAL_MIN_SCATTER_RATIO);
}

bool CompassCalibrator::solve9x9(float A[9][9], float b[9], float x[9]) {
    const int N = 9;
    float M[9][10];

    // Singularity has to be judged against the size of the matrix, not against
    // a fixed number. ATA here is a sum over every sample, so its entries grow
    // with the sample count; a flat 1e-9 is far below float precision at that
    // scale and would let a rank-deficient system through with a garbage
    // solution. Scale off the largest magnitude actually present.
    float scale = 0.0f;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            const float v = fabsf(A[i][j]);
            if (v > scale) scale = v;
        }
    }
    if (scale <= 0.0f) return false;
    const float pivot_tol = 1e-7f * scale;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) M[i][j] = A[i][j];
        M[i][N] = b[i];
    }

    for (int col = 0; col < N; col++) {
        int pivot_row = col;
        float pivot_val = fabsf(M[col][col]);
        for (int r = col + 1; r < N; r++) {
            float v = fabsf(M[r][col]);
            if (v > pivot_val) {
                pivot_val = v;
                pivot_row = r;
            }
        }
        if (pivot_val < pivot_tol) return false;
        if (pivot_row != col) {
            for (int c = 0; c <= N; c++) {
                float tmp = M[col][c];
                M[col][c] = M[pivot_row][c];
                M[pivot_row][c] = tmp;
            }
        }

        float inv_pivot = 1.0f / M[col][col];
        for (int c = col; c <= N; c++) M[col][c] *= inv_pivot;

        for (int r = 0; r < N; r++) {
            if (r == col) continue;
            float factor = M[r][col];
            if (factor == 0.0f) continue;
            for (int c = col; c <= N; c++) {
                M[r][c] -= factor * M[col][c];
            }
        }
    }

    for (int i = 0; i < N; i++) x[i] = M[i][N];
    return true;
}

bool CompassCalibrator::invert3x3(const Matrix3f &in, Matrix3f &out) {
    const float a = in.m[0][0], b = in.m[0][1], c = in.m[0][2];
    const float d = in.m[1][0], e = in.m[1][1], f = in.m[1][2];
    const float g = in.m[2][0], h = in.m[2][1], i = in.m[2][2];

    const float A = (e * i - f * h);
    const float B = -(d * i - f * g);
    const float C = (d * h - e * g);

    const float det = a * A + b * B + c * C;
    if (fabsf(det) < 1e-9f) return false;

    const float inv_det = 1.0f / det;

    out.m[0][0] = A * inv_det;
    out.m[0][1] = -(b * i - c * h) * inv_det;
    out.m[0][2] = (b * f - c * e) * inv_det;

    out.m[1][0] = B * inv_det;
    out.m[1][1] = (a * i - c * g) * inv_det;
    out.m[1][2] = -(a * f - c * d) * inv_det;

    out.m[2][0] = C * inv_det;
    out.m[2][1] = -(a * h - b * g) * inv_det;
    out.m[2][2] = (a * e - b * d) * inv_det;

    return true;
}

static inline void jacobiRotate(Matrix3f &m, Matrix3f &v, int p, int q) {
    if (fabsf(m.m[p][q]) < 1e-12f) return;

    float theta = (m.m[q][q] - m.m[p][p]) / (2.0f * m.m[p][q]);
    float t = (theta >= 0.0f) ? 1.0f / (theta + sqrtf(1.0f + theta * theta))
                              : -1.0f / (-theta + sqrtf(1.0f + theta * theta));
    float c = 1.0f / sqrtf(1.0f + t * t);
    float s = t * c;

    float mpp = m.m[p][p], mqq = m.m[q][q], mpq = m.m[p][q];
    m.m[p][p] = mpp - t * mpq;
    m.m[q][q] = mqq + t * mpq;
    m.m[p][q] = 0.0f;
    m.m[q][p] = 0.0f;

    for (int k = 0; k < 3; k++) {
        if (k != p && k != q) {
            float mkp = m.m[k][p], mkq = m.m[k][q];
            m.m[k][p] = m.m[p][k] = c * mkp - s * mkq;
            m.m[k][q] = m.m[q][k] = s * mkp + c * mkq;
        }
    }

    for (int k = 0; k < 3; k++) {
        float vkp = v.m[k][p], vkq = v.m[k][q];
        v.m[k][p] = c * vkp - s * vkq;
        v.m[k][q] = s * vkp + c * vkq;
    }
}

void CompassCalibrator::eigenSymmetric3x3(Matrix3f m, float eigval[3], Matrix3f &eigvec) {
    eigvec = Matrix3f::identity();

    // Convergence measured against the size of the matrix. An absolute
    // threshold would either never be met or be met immediately, depending
    // only on the units the caller happens to be working in.
    float scale = fabsf(m.m[0][0]) + fabsf(m.m[1][1]) + fabsf(m.m[2][2]);
    if (scale < 1e-30f) scale = 1.0f;
    const float off_tol = 1e-7f * scale;

    for (int sweep = 0; sweep < 12; sweep++) {
        float off = fabsf(m.m[0][1]) + fabsf(m.m[0][2]) + fabsf(m.m[1][2]);
        if (off < off_tol) break;

        jacobiRotate(m, eigvec, 0, 1);
        jacobiRotate(m, eigvec, 0, 2);
        jacobiRotate(m, eigvec, 1, 2);
    }
    eigval[0] = m.m[0][0];
    eigval[1] = m.m[1][1];
    eigval[2] = m.m[2][2];
}

float CompassCalibrator::fitResidual(const Vector3f &offset,
        const Matrix3f &softiron) const {
    if (_sample_count == 0) return 0.0f;

    const float R = (_nominal_radius > COMPASS_CAL_MIN_RADIUS) ? _nominal_radius : 1.0f;
    float acc = 0.0f;
    for (uint16_t i = 0; i < _sample_count; i++) {
        const Vector3f c = softiron.mul(_samples[i] - offset);
        const float mag = sqrtf(c.x * c.x + c.y * c.y + c.z * c.z);
        const float e = (mag - R) / R;
        acc += e * e;
    }
    return sqrtf(acc / (float)_sample_count);
}

CalStatus CompassCalibrator::calibrate() {
    if (_sample_count < COMPASS_CAL_MIN_SAMPLES) {
        _status = CalStatus::FAILED_NOT_ENOUGH_SAMPLES;
        return _status;
    }
    if (_filled_bins < COMPASS_CAL_MIN_BINS) {
        _status = CalStatus::FAILED_POOR_COVERAGE;
        return _status;
    }

    // Second, independent coverage test -- see scatterAnisotropy(). The bin
    // count above is measured from the centroid, which a one-sided sweep moves
    // into the middle of its own cap, so that test alone can be satisfied by an
    // airframe that was never turned over. This one looks at the shape of the
    // cloud, which no amount of hard iron can disguise.
    if (scatterAnisotropy() < COMPASS_CAL_MIN_SCATTER_RATIO) {
        // Recomputed rather than read from _scatter_ratio: calibrate() is public
        // and must not depend on addSample() having run since the last change.
        _status = CalStatus::FAILED_POOR_COVERAGE;
        return _status;
    }

    // Guard against a zero/negative radius only. The old guard was
    // `> 1.0f`, which silently clamped R to 1.0 for any sensor whose nominal
    // magnitude is below one in its native units -- the LIS3MDL reports Gauss,
    // so Earth's field is ~0.5 and every fit came out scaled by 1/0.5.
    const float R = (_nominal_radius > COMPASS_CAL_MIN_RADIUS) ? _nominal_radius : 1.0f;
    const float inv_R = 1.0f / R;

    // Fit about the centre of the cloud rather than about the sensor origin.
    // The form being fitted, x'Mx + 2n'x = 1, cannot represent an ellipsoid
    // that passes through the origin, and degrades either side of that: with
    // the hard iron near the field magnitude the fit is ill-conditioned (the
    // recovered centre was ~30% out), and beyond it M comes out negative
    // definite with K < 0 -- a perfectly valid ellipsoid that the sign checks
    // below would throw away. Pre-centring keeps the fit in the well behaved
    // K ~ 1 regime whatever the hard iron is, and the centre is simply added
    // back at the end.
    const Vector3f C = centroid();

    float ATA[9][9];
    float ATb[9];
    memset(ATA, 0, sizeof(ATA));
    memset(ATb, 0, sizeof(ATb));

    for (uint16_t k = 0; k < _sample_count; k++) {
        Vector3f s = (_samples[k] - C) * inv_R;
        float row[9] = {
            s.x * s.x, s.y * s.y, s.z * s.z,
            2.0f * s.x * s.y, 2.0f * s.x * s.z, 2.0f * s.y * s.z,
            2.0f * s.x, 2.0f * s.y, 2.0f * s.z
        };
        for (int i = 0; i < 9; i++) {
            ATb[i] += row[i];
            for (int j = 0; j < 9; j++) {
                ATA[i][j] += row[i] * row[j];
            }
        }
    }

    float p[9];
    if (!solve9x9(ATA, ATb, p)) {
        _status = CalStatus::FAILED_SINGULAR_MATRIX;
        return _status;
    }

    Matrix3f M;
    M.m[0][0] = p[0]; M.m[0][1] = p[3]; M.m[0][2] = p[4];
    M.m[1][0] = p[3]; M.m[1][1] = p[1]; M.m[1][2] = p[5];
    M.m[2][0] = p[4]; M.m[2][1] = p[5]; M.m[2][2] = p[2];

    Vector3f n(p[6], p[7], p[8]);

    Matrix3f Minv;
    if (!invert3x3(M, Minv)) {
        _status = CalStatus::FAILED_SINGULAR_MATRIX;
        return _status;
    }
    Vector3f V = Minv.mul(n) * -1.0f;

    // Completing the square gives (x-V)'M(x-V) = K, so the surface the samples
    // lie on is (x-V)'(M/K)(x-V) = 1 and the matrix to take the square root of
    // is M/K, not M. Only their ratio has to be positive definite: M and K may
    // both be negative and still describe a real ellipsoid. Folding K in here
    // rather than applying 1/sqrt(K) afterwards handles either sign.
    float K = 1.0f - n.dot(V);
    if (fabsf(K) < 1e-6f) {
        _status = CalStatus::FAILED_DEGENERATE_ELLIPSOID;
        return _status;
    }
    Matrix3f MK = M.scaled(1.0f / K);

    float eigval[3];
    Matrix3f eigvec;
    eigenSymmetric3x3(MK, eigval, eigvec);

    // Reject anything that is not a genuine ellipsoid: a non-positive
    // eigenvalue means a hyperboloid or a paraboloid, and an eigenvalue
    // vanishing relative to the largest means a degenerate, effectively flat
    // one. The relative floor is what makes this test independent of units.
    float max_eig = eigval[0];
    for (int i = 1; i < 3; i++) {
        if (eigval[i] > max_eig) max_eig = eigval[i];
    }
    if (max_eig <= 0.0f) {
        _status = CalStatus::FAILED_DEGENERATE_ELLIPSOID;
        return _status;
    }
    const float min_eig = 1e-6f * max_eig;
    for (int i = 0; i < 3; i++) {
        if (eigval[i] < min_eig) {
            _status = CalStatus::FAILED_DEGENERATE_ELLIPSOID;
            return _status;
        }
    }

    // Every element written explicitly. Only the diagonal carries a value, and
    // the square root below is wrong unless the rest are zero, so this must
    // not depend on how common.hpp happens to default-construct a Matrix3f --
    // `Matrix3f D{}` zeroes an aggregate but calls a user-provided default
    // constructor, which may leave the members uninitialised.
    Matrix3f D = Matrix3f();
    D.m[0][0] = 0.0f; D.m[0][1] = 0.0f; D.m[0][2] = 0.0f;
    D.m[1][0] = 0.0f; D.m[1][1] = 0.0f; D.m[1][2] = 0.0f;
    D.m[2][0] = 0.0f; D.m[2][1] = 0.0f; D.m[2][2] = 0.0f;
    D.m[0][0] = sqrtf(eigval[0]);
    D.m[1][1] = sqrtf(eigval[1]);
    D.m[2][2] = sqrtf(eigval[2]);

    // Symmetric (rotation-free) square root of M/K: it maps the ellipsoid onto
    // a sphere without also rotating the sensor frame, which is what keeps the
    // corrected axes aligned with the physical ones.
    const Matrix3f softiron = eigvec.mul(D).mul(eigvec.transposed());
    const Vector3f offset = V * R + C;

    // Last gate, and the only one that looks at how well the answer fits the
    // data rather than at the shape of the algebra. Everything above can pass
    // on readings that carry no field information at all.
    // Recorded whether it passes or fails: on a failure it is the only number
    // that says how badly, and on a success it is a quality score worth seeing.
    _last_residual = fitResidual(offset, softiron);
    if (_last_residual > COMPASS_CAL_MAX_FIT_RESIDUAL) {
        _status = CalStatus::FAILED_POOR_FIT;
        return _status;
    }

    _softiron = softiron;
    _offset = offset;

    _status = CalStatus::SUCCESS;
    return _status;
}

Vector3f CompassCalibrator::correct(const Vector3f &raw) const {
    if (_status != CalStatus::SUCCESS) {
        return raw;
    }
    return _softiron.mul(raw - _offset);
}