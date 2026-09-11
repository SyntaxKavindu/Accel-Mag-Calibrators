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
}

void AccelerometerCalibrator::beginSixPosition(float motion_threshold) {
    reset();
    // A zero or negative threshold would reject every sample after the first
    // and stall the calibration forever; keep reset()'s default in that case.
    if (motion_threshold > 0.0f) {
        _sp_motion_threshold_sq = motion_threshold * motion_threshold;
    }
    _status = AccelCalStatus::IN_PROGRESS;
}

void AccelerometerCalibrator::startPosition(AccelPosition pos) {
    if (_status != AccelCalStatus::IN_PROGRESS) return;
    size_t idx = (size_t)pos;
    _current_pos = pos;
    _pos_sum[idx] = Vector3f();
    _pos_count[idx] = 0;
    _pos_done[idx] = false;
    _sp_collecting = true;
}

AccelSampleResult AccelerometerCalibrator::addSample(float x, float y, float z) {
    const Vector3f s(x, y, z);
    const size_t idx = (size_t)_current_pos;

    AccelSampleResult r = AccelSampleResult::ACCEPTED;
    if (_pos_done[idx]) {
        r = AccelSampleResult::REJECTED_POSITION_DONE;
    } else if (!_sp_collecting) {
        r = AccelSampleResult::REJECTED_NOT_STARTED;
    } else if (_pos_count[idx] == 0) {
        _pos_sum[idx] = s;
        _pos_count[idx] = 1;
        _sp_running_mean = s;
    } else {
        const Vector3f d = s - _sp_running_mean;
        const float dist_sq = d.x * d.x + d.y * d.y + d.z * d.z;
        if (dist_sq > _sp_motion_threshold_sq) {
            // Discard the whole average rather than folding the outlier in.
            // The running mean is deliberately NOT reset to the outlier: doing
            // that would make the next sample look still relative to a value
            // taken mid-movement.
            _pos_sum[idx] = Vector3f();
            _pos_count[idx] = 0;
            r = AccelSampleResult::REJECTED_MOTION;
        } else {
            _pos_sum[idx] = _pos_sum[idx] + s;
            _pos_count[idx]++;
            _sp_running_mean = _pos_sum[idx] / (float)_pos_count[idx];

            if (_pos_count[idx] >= ACCEL_CAL_SAMPLES_PER_POSITION) {
                _pos_avg[idx] = _pos_sum[idx] / (float)_pos_count[idx];
                _pos_done[idx] = true;
                _sp_collecting = false;
                r = AccelSampleResult::ACCEPTED_POSITION_DONE;
            }
        }
    }

    // Counted on every sample the caller offers, accepted or not: a long run
    // of rejected samples is exactly what this has to notice.
    noteProgress();
    return r;
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
    return _samples_since_progress >= ACCEL_CAL_SIXPOS_STALL_LIMIT;
}

bool AccelerometerCalibrator::isReadyToCalibrate() const {
    return allPositionsComplete();
}

float AccelerometerCalibrator::getProgressPercent() const {
    uint32_t total = 0;
    for (size_t i = 0; i < (size_t)AccelPosition::NUM_POSITIONS; i++) {
        total += _pos_done[i] ? ACCEL_CAL_SAMPLES_PER_POSITION : _pos_count[i];
    }
    const float pct = 100.0f * (float)total
            / (float)(ACCEL_CAL_SAMPLES_PER_POSITION * (int)AccelPosition::NUM_POSITIONS);
    return pct > 100.0f ? 100.0f : pct;
}

AccelCalStatus AccelerometerCalibrator::calibrate() {
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

    // With the axis up and down, gravity enters that axis at +1 g and -1 g, so
    // the midpoint of the pair is the bias and half the difference is the
    // sensitivity. The other two components of each pair carry the cross-axis
    // terms; they are deliberately not used -- see the note in the header on
    // why the rotation half of the sensor model belongs to LevelCalibrator.
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

Vector3f AccelerometerCalibrator::correct(const Vector3f &raw) const {
    if (_status != AccelCalStatus::SUCCESS) return raw;
    return Vector3f(
        (raw.x - _bias_sixpos.x) / _scale_sixpos.x,
        (raw.y - _bias_sixpos.y) / _scale_sixpos.y,
        (raw.z - _bias_sixpos.z) / _scale_sixpos.z);
}

Vector3f AccelerometerCalibrator::getBias() const { return _bias_sixpos; }
Vector3f AccelerometerCalibrator::getScale() const { return _scale_sixpos; }

// The same correction correct() applies, as one matrix, so the facade can hold
// a single affine form for both sensors. Diagonal by construction: this fit
// produces no rotation and no cross-axis term. Identity until the fit has
// succeeded, so an unfinished calibration cannot scale anything.
Matrix3f AccelerometerCalibrator::getMatrix() const {
    Matrix3f m = Matrix3f::identity();
    if (_status != AccelCalStatus::SUCCESS) return m;
    m.m[0][0] = 1.0f / _scale_sixpos.x;
    m.m[1][1] = 1.0f / _scale_sixpos.y;
    m.m[2][2] = 1.0f / _scale_sixpos.z;
    return m;
}
