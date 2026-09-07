/*
 * AccelerometerCalibrator.cpp
 *
 *  Created on: Aug 27, 2026
 *      Author: KAVINDU
 */

#include "AccelerometerCalibrator.hpp"

AccelerometerCalibrator::AccelerometerCalibrator() {
    reset();
}

void AccelerometerCalibrator::reset() {
    _mode = AccelCalMode::NONE;
    _status = AccelCalStatus::IDLE;

    for (size_t i = 0; i < (size_t)AccelPosition::NUM_POSITIONS; i++) {
        _pos_sum[i] = Vector3f();
        _pos_count[i] = 0;
        _pos_done[i] = false;
        _pos_avg[i] = Vector3f();
    }
    _current_pos = AccelPosition::X_UP;
    _sp_collecting = false;
    _sp_running_mean = Vector3f();
    _sp_motion_threshold_sq = 0.5f * 0.5f;
    _bias_sixpos = Vector3f();
    _scale_sixpos = Vector3f(1.0f, 1.0f, 1.0f);

    _best_progress = -1.0f;
    _samples_since_progress = 0;

    _samples = nullptr;
    _capacity = 0;
    _sample_count = 0;
    _filled_bins = 0;
    memset(_bin_count, 0, sizeof(_bin_count));
    _nominal_radius = ACCEL_CAL_STANDARD_GRAVITY;
    _window_count = 0;
    _window_head = 0;
    _stillness_threshold_sq = 0.2f * 0.2f;
    _offset_tumble = Vector3f();
    _matrix_tumble = Matrix3f::identity();
}

void AccelerometerCalibrator::beginSixPosition(float motion_threshold) {
    reset();
    _mode = AccelCalMode::SIX_POSITION;
    // A zero or negative threshold would reject every sample after the first
    // and stall the calibration forever; keep reset()'s default in that case.
    if (motion_threshold > 0.0f) {
        _sp_motion_threshold_sq = motion_threshold * motion_threshold;
    }
    _status = AccelCalStatus::IN_PROGRESS;
}

bool AccelerometerCalibrator::beginTumble(float nominal_g, float stillness_threshold,
                                          Vector3f *sample_buffer, uint16_t capacity) {
    reset();

    // A buffer too small to reach the minimum sample count would collect
    // happily and then fail the fit for ever, with progress stuck just under
    // 100%. Refuse at the start, where the caller can still report it, rather
    // than at the end where it looks like an operator error.
    if (sample_buffer == nullptr || capacity < ACCEL_CAL_TUMBLE_MIN_SAMPLES) {
        return false;
    }
    _samples  = sample_buffer;
    _capacity = capacity;

    _mode = AccelCalMode::TUMBLE;
    // _nominal_radius normalises the ellipsoid fit and divides in correct(),
    // so it has to be strictly positive. Validating it here means calibrate()
    // and correct() can both use it directly and can never disagree about the
    // scale they applied.
    if (nominal_g > ACCEL_CAL_MIN_RADIUS) {
        _nominal_radius = nominal_g;
    }
    if (stillness_threshold > 0.0f) {
        _stillness_threshold_sq = stillness_threshold * stillness_threshold;
    }
    _status = AccelCalStatus::IN_PROGRESS;
    return true;
}

void AccelerometerCalibrator::startPosition(AccelPosition pos) {
    if (_mode != AccelCalMode::SIX_POSITION) return;
    size_t idx = (size_t)pos;
    _current_pos = pos;
    _pos_sum[idx] = Vector3f();
    _pos_count[idx] = 0;
    _pos_done[idx] = false;
    _sp_collecting = true;
}

AccelSampleResult AccelerometerCalibrator::addSampleSixPosition(const Vector3f &s) {
    size_t idx = (size_t)_current_pos;

    if (_pos_done[idx]) return AccelSampleResult::REJECTED_POSITION_DONE;
    if (!_sp_collecting) return AccelSampleResult::REJECTED_NOT_STARTED;

    if (_pos_count[idx] == 0) {
        _pos_sum[idx] = s;
        _pos_count[idx] = 1;
        _sp_running_mean = s;
    } else {
        Vector3f d = s - _sp_running_mean;
        float dist_sq = d.x * d.x + d.y * d.y + d.z * d.z;
        if (dist_sq > _sp_motion_threshold_sq) {
            // Fix: Reset accumulated position cleanly on motion without setting outlier sample as running mean
            _pos_sum[idx] = Vector3f();
            _pos_count[idx] = 0;
            return AccelSampleResult::REJECTED_MOTION;
        }
        _pos_sum[idx] = _pos_sum[idx] + s;
        _pos_count[idx]++;
        _sp_running_mean = _pos_sum[idx] / (float)_pos_count[idx];
    }

    if (_pos_count[idx] >= ACCEL_CAL_SAMPLES_PER_POSITION) {
        _pos_avg[idx] = _pos_sum[idx] / (float)_pos_count[idx];
        _pos_done[idx] = true;
        _sp_collecting = false;
        return AccelSampleResult::ACCEPTED_POSITION_DONE;
    }
    return AccelSampleResult::ACCEPTED;
}

bool AccelerometerCalibrator::isPositionDone(AccelPosition pos) const {
    return _pos_done[(size_t)pos];
}

bool AccelerometerCalibrator::allPositionsComplete() const {
    for (size_t i = 0; i < (size_t)AccelPosition::NUM_POSITIONS; i++) {
        if (!_pos_done[i]) return false;
    }
    return true;
}

AccelCalStatus AccelerometerCalibrator::calibrateSixPosition() {
    if (!allPositionsComplete()) {
        _status = AccelCalStatus::FAILED_NOT_ENOUGH_SAMPLES;
        return _status;
    }
    const Vector3f &xu = _pos_avg[(size_t)AccelPosition::X_UP];
    const Vector3f &xd = _pos_avg[(size_t)AccelPosition::X_DOWN];
    const Vector3f &yu = _pos_avg[(size_t)AccelPosition::Y_UP];
    const Vector3f &yd = _pos_avg[(size_t)AccelPosition::Y_DOWN];
    const Vector3f &zu = _pos_avg[(size_t)AccelPosition::Z_UP];
    const Vector3f &zd = _pos_avg[(size_t)AccelPosition::Z_DOWN];

    _bias_sixpos.x = (xu.x + xd.x) * 0.5f;
    _bias_sixpos.y = (yu.y + yd.y) * 0.5f;
    _bias_sixpos.z = (zu.z + zd.z) * 0.5f;

    _scale_sixpos.x = (xu.x - xd.x) * 0.5f;
    _scale_sixpos.y = (yu.y - yd.y) * 0.5f;
    _scale_sixpos.z = (zu.z - zd.z) * 0.5f;

    // correct() divides by these, so a dead or stuck axis has to be caught
    // here. The guard is relative rather than a fixed "< 1.0", which was only
    // meaningful for raw-LSB samples and sat exactly on the boundary for
    // samples already expressed in g. A healthy tri-axis accelerometer has
    // per-axis sensitivities within a few percent of one another, so anything
    // under a tenth of the largest is a broken axis in any unit system.
    const float ax = fabsf(_scale_sixpos.x);
    const float ay = fabsf(_scale_sixpos.y);
    const float az = fabsf(_scale_sixpos.z);
    const float max_scale = fmaxf(ax, fmaxf(ay, az));
    if (max_scale < ACCEL_CAL_MIN_RADIUS) {
        _status = AccelCalStatus::FAILED_BAD_GEOMETRY;
        return _status;
    }
    const float min_scale = 0.1f * max_scale;
    if (ax < min_scale || ay < min_scale || az < min_scale) {
        _status = AccelCalStatus::FAILED_BAD_GEOMETRY;
        return _status;
    }

    // An accelerometer reads specific force, so with the +X axis pointing up
    // it reports +1 g on X and (xu.x - xd.x) is positive. A negative scale
    // therefore means X_UP and X_DOWN were captured the wrong way round. The
    // magnitude checks above use fabsf and cannot see this, and the fit would
    // otherwise "succeed" while silently negating that axis -- correcting a
    // true +1 g to -1 g, which flies inverted. Reject it instead.
    if (_scale_sixpos.x <= 0.0f || _scale_sixpos.y <= 0.0f
            || _scale_sixpos.z <= 0.0f) {
        _status = AccelCalStatus::FAILED_BAD_GEOMETRY;
        return _status;
    }

    _status = AccelCalStatus::SUCCESS;
    return _status;
}

float AccelerometerCalibrator::getMaxMisalignmentDeg() const {
    // correct() is a pass-through until the fit has succeeded, so without this
    // status check the function would silently report an angle derived from
    // uncorrected samples instead of the documented "unavailable" -1.
    if (_mode != AccelCalMode::SIX_POSITION || !allPositionsComplete()) return -1.0f;
    if (_status != AccelCalStatus::SUCCESS) return -1.0f;

    struct Entry { AccelPosition pos; int primary_axis; };
    static const Entry entries[6] = {
        {AccelPosition::X_UP, 0}, {AccelPosition::X_DOWN, 0},
        {AccelPosition::Y_UP, 1}, {AccelPosition::Y_DOWN, 1},
        {AccelPosition::Z_UP, 2}, {AccelPosition::Z_DOWN, 2},
    };
    float max_deg = 0.0f;
    for (const Entry &e : entries) {
        const Vector3f &raw_avg = _pos_avg[(size_t)e.pos];
        Vector3f c = correct(raw_avg);
        float primary = (e.primary_axis == 0) ? c.x : (e.primary_axis == 1 ? c.y : c.z);
        float o1 = (e.primary_axis == 0) ? c.y : (e.primary_axis == 1) ? c.x : c.x;
        float o2 = (e.primary_axis == 0) ? c.z : (e.primary_axis == 1) ? c.z : c.y;
        float off_mag = sqrtf(o1 * o1 + o2 * o2);
        if (fabsf(primary) < 1e-6f) continue;
        float deg = atanf(off_mag / fabsf(primary)) * 57.29577951f;
        if (deg > max_deg) max_deg = deg;
    }
    return max_deg;
}

bool AccelerometerCalibrator::isWindowStill(Vector3f &out_mean) const {
    if (_window_count < ACCEL_CAL_STILLNESS_WINDOW) return false;
    Vector3f sum = Vector3f();
    for (uint8_t i = 0; i < ACCEL_CAL_STILLNESS_WINDOW; i++) sum = sum + _window[i];
    Vector3f mean = sum / (float)ACCEL_CAL_STILLNESS_WINDOW;
    float sq_dev_sum = 0.0f;
    for (uint8_t i = 0; i < ACCEL_CAL_STILLNESS_WINDOW; i++) {
        Vector3f d = _window[i] - mean;
        sq_dev_sum += d.x * d.x + d.y * d.y + d.z * d.z;
    }
    float mean_sq_dev = sq_dev_sum / (float)ACCEL_CAL_STILLNESS_WINDOW;
    if (mean_sq_dev > _stillness_threshold_sq) return false;
    out_mean = mean;
    return true;
}

Vector3f AccelerometerCalibrator::centroid() const {
    if (_sample_count == 0) return Vector3f();
    Vector3f sum = Vector3f();
    for (uint16_t i = 0; i < _sample_count; i++) sum = sum + _samples[i];
    return sum / (float)_sample_count;
}

uint8_t AccelerometerCalibrator::getBinIndex(const Vector3f &s) {
    float mag = sqrtf(s.x * s.x + s.y * s.y + s.z * s.z);
    if (mag < 1e-6f) return 0;
    
    float zn = s.z / mag;
    int elev = (int)((zn + 1.0f) * 4.0f);
    if (elev < 0) elev = 0; else if (elev > 7) elev = 7;

    float az_rad = atan2f(s.y, s.x);
    int az = (int)((az_rad + 3.1415926535f) * (8.0f / 6.283185307f));
    if (az < 0) az = 0; else if (az > 7) az = 7;

    return (uint8_t)(elev * 8 + az);
}

AccelSampleResult AccelerometerCalibrator::addSampleTumble(const Vector3f &s) {
    if (_samples == nullptr || _sample_count >= _capacity) {
        return AccelSampleResult::REJECTED_BUFFER_FULL;
    }
    _window[_window_head] = s;
    _window_head = (uint8_t)((_window_head + 1) % ACCEL_CAL_STILLNESS_WINDOW);
    if (_window_count < ACCEL_CAL_STILLNESS_WINDOW) _window_count++;

    Vector3f settled_mean;
    if (!isWindowStill(settled_mean)) return AccelSampleResult::REJECTED_MOTION;

    // Binned about the origin, unlike the compass, which bins about the centre
    // of its cloud. That distinction only matters when the offset is
    // comparable to the radius being measured: hard iron routinely exceeds
    // Earth's field, but an accelerometer bias is always a small fraction of
    // 1 g, and one large enough to distort this coverage measure would mean a
    // dead part. Binning about the centre here was measured and costs ~40%
    // more tumbling for no reachable benefit -- the stillness window makes
    // this input far more repetitive than the compass's, so the rebinning it
    // requires discards much more of the work.
    uint8_t bin = getBinIndex(settled_mean);
    if (_bin_count[bin] >= ACCEL_CAL_TUMBLE_MAX_PER_BIN) {
        return AccelSampleResult::REJECTED_TOO_CLOSE;
    }

    if (_bin_count[bin] == 0) _filled_bins++;
    _bin_count[bin]++;

    _samples[_sample_count] = settled_mean;
    _sample_count++;

    return AccelSampleResult::ACCEPTED;
}

bool AccelerometerCalibrator::solve9x9(float A[9][9], float b[9], float x[9]) {
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

    for (int i = 0; i < N; i++) { for (int j = 0; j < N; j++) M[i][j] = A[i][j]; M[i][N] = b[i]; }
    for (int col = 0; col < N; col++) {
        int pivot_row = col;
        float pivot_val = fabsf(M[col][col]);
        for (int r = col + 1; r < N; r++) {
            float v = fabsf(M[r][col]);
            if (v > pivot_val) { pivot_val = v; pivot_row = r; }
        }
        if (pivot_val < pivot_tol) return false;
        if (pivot_row != col) {
            for (int c = 0; c <= N; c++) { float t=M[col][c]; M[col][c]=M[pivot_row][c]; M[pivot_row][c]=t; }
        }
        float inv_pivot = 1.0f / M[col][col];
        for (int c = col; c <= N; c++) M[col][c] *= inv_pivot;
        for (int r = 0; r < N; r++) {
            if (r == col) continue;
            float factor = M[r][col];
            if (factor == 0.0f) continue;
            for (int c = col; c <= N; c++) M[r][c] -= factor * M[col][c];
        }
    }
    for (int i = 0; i < N; i++) x[i] = M[i][N];
    return true;
}

bool AccelerometerCalibrator::invert3x3(const Matrix3f &in, Matrix3f &out) {
    const float a=in.m[0][0], b=in.m[0][1], c=in.m[0][2];
    const float d=in.m[1][0], e=in.m[1][1], f=in.m[1][2];
    const float g=in.m[2][0], h=in.m[2][1], i=in.m[2][2];
    const float A=(e*i-f*h), B=-(d*i-f*g), C=(d*h-e*g);
    const float det = a*A + b*B + c*C;
    if (fabsf(det) < 1e-9f) return false;
    const float inv_det = 1.0f / det;
    out.m[0][0]=A*inv_det; out.m[0][1]=-(b*i-c*h)*inv_det; out.m[0][2]=(b*f-c*e)*inv_det;
    out.m[1][0]=B*inv_det; out.m[1][1]=(a*i-c*g)*inv_det;  out.m[1][2]=-(a*f-c*d)*inv_det;
    out.m[2][0]=C*inv_det; out.m[2][1]=-(a*h-b*g)*inv_det; out.m[2][2]=(a*e-b*d)*inv_det;
    return true;
}

static inline void jacobiRotateU(Matrix3f &m, Matrix3f &v, int p, int q) {
    if (fabsf(m.m[p][q]) < 1e-12f) return;
    float theta = (m.m[q][q] - m.m[p][p]) / (2.0f * m.m[p][q]);
    float t = (theta >= 0.0f) ? 1.0f/(theta+sqrtf(1.0f+theta*theta)) : -1.0f/(-theta+sqrtf(1.0f+theta*theta));
    float c = 1.0f/sqrtf(1.0f+t*t), s = t*c;
    float mpp=m.m[p][p], mqq=m.m[q][q], mpq=m.m[p][q];
    m.m[p][p]=mpp-t*mpq; m.m[q][q]=mqq+t*mpq; m.m[p][q]=0.0f; m.m[q][p]=0.0f;
    for (int k=0;k<3;k++) if (k!=p && k!=q) {
        float mkp=m.m[k][p], mkq=m.m[k][q];
        m.m[k][p]=m.m[p][k]=c*mkp-s*mkq;
        m.m[k][q]=m.m[q][k]=s*mkp+c*mkq;
    }
    for (int k=0;k<3;k++) {
        float vkp=v.m[k][p], vkq=v.m[k][q];
        v.m[k][p]=c*vkp-s*vkq;
        v.m[k][q]=s*vkp+c*vkq;
    }
}

void AccelerometerCalibrator::eigenSymmetric3x3(Matrix3f m, float eigval[3], Matrix3f &eigvec) {
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
        jacobiRotateU(m, eigvec, 0, 1);
        jacobiRotateU(m, eigvec, 0, 2);
        jacobiRotateU(m, eigvec, 1, 2);
    }
    eigval[0]=m.m[0][0]; eigval[1]=m.m[1][1]; eigval[2]=m.m[2][2];
}

float AccelerometerCalibrator::fitResidual(const Vector3f &offset,
        const Matrix3f &matrix) const {
    if (_sample_count == 0) return 0.0f;

    const float R = _nominal_radius;
    float acc = 0.0f;
    for (uint16_t i = 0; i < _sample_count; i++) {
        const Vector3f c = matrix.mul(_samples[i] - offset);
        const float mag = sqrtf(c.x * c.x + c.y * c.y + c.z * c.z);
        const float e = (mag - R) / R;
        acc += e * e;
    }
    return sqrtf(acc / (float)_sample_count);
}

AccelCalStatus AccelerometerCalibrator::calibrateTumble() {
    if (_sample_count < ACCEL_CAL_TUMBLE_MIN_SAMPLES) {
        _status = AccelCalStatus::FAILED_NOT_ENOUGH_SAMPLES;
        return _status;
    }
    if (_filled_bins < ACCEL_CAL_TUMBLE_MIN_BINS) {
        _status = AccelCalStatus::FAILED_POOR_COVERAGE;
        return _status;
    }

    // beginTumble() already guarantees _nominal_radius > 0. The previous
    // `> 1.0f` guard silently clamped R to 1.0 for any sensor whose 1 g
    // reading is below one in its native units, which desynchronised the fit
    // from the division correct() performs and scaled the output wrongly.
    const float R = _nominal_radius;
    const float inv_R = 1.0f / R;

    // Fit about the centre of the cloud rather than about the sensor origin.
    // The form being fitted, x'Mx + 2n'x = 1, cannot represent an ellipsoid
    // that passes through the origin and is ill-conditioned near that case,
    // which for this sensor means a bias approaching 1 g. Pre-centring keeps
    // the fit in the well behaved K ~ 1 regime whatever the bias is, and the
    // centre is simply added back at the end.
    const Vector3f C = centroid();

    float ATA[9][9]; float ATb[9];
    memset(ATA, 0, sizeof(ATA)); memset(ATb, 0, sizeof(ATb));

    for (uint16_t k = 0; k < _sample_count; k++) {
        Vector3f s = (_samples[k] - C) * inv_R;
        float row[9] = { s.x*s.x, s.y*s.y, s.z*s.z, 2*s.x*s.y, 2*s.x*s.z, 2*s.y*s.z, 2*s.x, 2*s.y, 2*s.z };
        for (int i = 0; i < 9; i++) {
            ATb[i] += row[i];
            for (int j = 0; j < 9; j++) ATA[i][j] += row[i]*row[j];
        }
    }

    float p[9];
    if (!solve9x9(ATA, ATb, p)) { _status = AccelCalStatus::FAILED_SINGULAR_MATRIX; return _status; }

    Matrix3f M;
    M.m[0][0]=p[0]; M.m[0][1]=p[3]; M.m[0][2]=p[4];
    M.m[1][0]=p[3]; M.m[1][1]=p[1]; M.m[1][2]=p[5];
    M.m[2][0]=p[4]; M.m[2][1]=p[5]; M.m[2][2]=p[2];
    Vector3f n(p[6], p[7], p[8]);

    Matrix3f Minv;
    if (!invert3x3(M, Minv)) { _status = AccelCalStatus::FAILED_SINGULAR_MATRIX; return _status; }
    Vector3f V = Minv.mul(n) * -1.0f;

    // Completing the square gives (x-V)'M(x-V) = K, so the surface the samples
    // lie on is (x-V)'(M/K)(x-V) = 1 and the matrix to take the square root of
    // is M/K, not M. Only their ratio has to be positive definite: M and K may
    // both be negative and still describe a real ellipsoid. Folding K in here
    // rather than applying 1/sqrt(K) afterwards handles either sign.
    float K = 1.0f - (n.x*V.x + n.y*V.y + n.z*V.z);
    if (fabsf(K) < 1e-6f) { _status = AccelCalStatus::FAILED_DEGENERATE_ELLIPSOID; return _status; }
    Matrix3f MK = M.scaled(1.0f / K);

    float eigval[3]; Matrix3f eigvec;
    eigenSymmetric3x3(MK, eigval, eigvec);

    // Reject anything that is not a genuine ellipsoid: a non-positive
    // eigenvalue means a hyperboloid or a paraboloid, and an eigenvalue
    // vanishing relative to the largest means a degenerate, effectively flat
    // one. The relative floor is what makes this test independent of units.
    float max_eig = eigval[0];
    for (int i = 1; i < 3; i++) { if (eigval[i] > max_eig) max_eig = eigval[i]; }
    if (max_eig <= 0.0f) { _status = AccelCalStatus::FAILED_DEGENERATE_ELLIPSOID; return _status; }
    const float min_eig = 1e-6f * max_eig;
    for (int i = 0; i < 3; i++) {
        if (eigval[i] < min_eig) { _status = AccelCalStatus::FAILED_DEGENERATE_ELLIPSOID; return _status; }
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
    D.m[0][0]=sqrtf(eigval[0]); D.m[1][1]=sqrtf(eigval[1]); D.m[2][2]=sqrtf(eigval[2]);

    // Symmetric (rotation-free) square root of M/K: it maps the ellipsoid onto
    // a sphere without also rotating the sensor frame, which is what keeps the
    // corrected axes aligned with the physical ones.
    const Matrix3f matrix = eigvec.mul(D).mul(eigvec.transposed());
    const Vector3f offset = V * R + C;

    // Last gate, and the only one that looks at how well the answer fits the
    // data rather than at the shape of the algebra.
    if (fitResidual(offset, matrix) > ACCEL_CAL_MAX_FIT_RESIDUAL) {
        _status = AccelCalStatus::FAILED_POOR_FIT;
        return _status;
    }

    _matrix_tumble = matrix;
    _offset_tumble = offset;

    _status = AccelCalStatus::SUCCESS;
    return _status;
}

AccelSampleResult AccelerometerCalibrator::addSample(float x, float y, float z) {
    Vector3f s(x, y, z);
    AccelSampleResult r;
    if (_mode == AccelCalMode::SIX_POSITION)   r = addSampleSixPosition(s);
    else if (_mode == AccelCalMode::TUMBLE)    r = addSampleTumble(s);
    else return AccelSampleResult::REJECTED_NOT_STARTED;

    // Counted on every sample the caller offers, accepted or not: a long run
    // of rejected samples is exactly what this has to notice.
    noteProgress();
    return r;
}

void AccelerometerCalibrator::noteProgress() {
    const float p = getProgressPercent();
    if (p > _best_progress) {
        _best_progress = p;
        _samples_since_progress = 0;
    } else if (_samples_since_progress < 0xFFFFFFFFU) {
        _samples_since_progress++;
    }
}

bool AccelerometerCalibrator::isStalled() const {
    if (_mode == AccelCalMode::SIX_POSITION) {
        return _samples_since_progress >= ACCEL_CAL_SIXPOS_STALL_LIMIT;
    }
    if (_mode == AccelCalMode::TUMBLE) {
        return _samples_since_progress >= ACCEL_CAL_TUMBLE_STALL_LIMIT;
    }
    return false;
}

bool AccelerometerCalibrator::isReadyToCalibrate() const {
    if (_mode == AccelCalMode::SIX_POSITION) return allPositionsComplete();
    if (_mode == AccelCalMode::TUMBLE) {
        return (_sample_count >= ACCEL_CAL_TUMBLE_MIN_SAMPLES) && (_filled_bins >= ACCEL_CAL_TUMBLE_MIN_BINS);
    }
    return false;
}

float AccelerometerCalibrator::getProgressPercent() const {
    if (_mode == AccelCalMode::SIX_POSITION) {
        uint32_t total = 0;
        for (size_t i = 0; i < (size_t)AccelPosition::NUM_POSITIONS; i++) {
            total += _pos_done[i] ? ACCEL_CAL_SAMPLES_PER_POSITION : _pos_count[i];
        }
        float pct = 100.0f * (float)total / (float)(ACCEL_CAL_SAMPLES_PER_POSITION * (int)AccelPosition::NUM_POSITIONS);
        return pct > 100.0f ? 100.0f : pct;
    }
    if (_mode == AccelCalMode::TUMBLE) {
        float sample_frac = (float)_sample_count / (float)ACCEL_CAL_TUMBLE_MIN_SAMPLES;
        if (sample_frac > 1.0f) sample_frac = 1.0f;
        float bin_frac = (float)_filled_bins / (float)ACCEL_CAL_TUMBLE_MIN_BINS;
        if (bin_frac > 1.0f) bin_frac = 1.0f;
        float combined = sample_frac < bin_frac ? sample_frac : bin_frac;
        return combined * 100.0f;
    }
    return 0.0f;
}

AccelCalStatus AccelerometerCalibrator::calibrate() {
    if (_mode == AccelCalMode::SIX_POSITION) return calibrateSixPosition();
    if (_mode == AccelCalMode::TUMBLE) return calibrateTumble();
    _status = AccelCalStatus::FAILED_NOT_ENOUGH_SAMPLES;
    return _status;
}

Vector3f AccelerometerCalibrator::correct(const Vector3f &raw) const {
    if (_status != AccelCalStatus::SUCCESS) return raw;
    if (_mode == AccelCalMode::SIX_POSITION) {
        return Vector3f(
            (raw.x - _bias_sixpos.x) / _scale_sixpos.x,
            (raw.y - _bias_sixpos.y) / _scale_sixpos.y,
            (raw.z - _bias_sixpos.z) / _scale_sixpos.z);
    }

    if (_mode == AccelCalMode::TUMBLE) {
        Vector3f corrected_raw_units = _matrix_tumble.mul(raw - _offset_tumble);
        return corrected_raw_units / _nominal_radius;
    }
    return raw;
}

Vector3f AccelerometerCalibrator::getBias() const {
    if (_mode == AccelCalMode::SIX_POSITION) return _bias_sixpos;
    if (_mode == AccelCalMode::TUMBLE) return _offset_tumble;
    return Vector3f();
}

Vector3f AccelerometerCalibrator::getScale() const { return _scale_sixpos; }
Matrix3f AccelerometerCalibrator::getMatrix() const { return _matrix_tumble; }
