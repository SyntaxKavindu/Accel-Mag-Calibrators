/*
 * LevelCalibrator.cpp
 *
 *  Created on: Sep 5, 2026
 *      Author: KAVINDU
 */

#include "LevelCalibrator.hpp"

LevelCalibrator::LevelCalibrator() {
    reset();
}

void LevelCalibrator::reset() {
    _status = LevelCalStatus::IDLE;
    _target = Vector3f(0.0f, 0.0f, 1.0f);
    _nominal_g = LEVEL_CAL_STANDARD_GRAVITY;
    _yaw_offset_rad = 0.0f;
    _motion_threshold_sq = LEVEL_CAL_MOTION_THRESHOLD * LEVEL_CAL_MOTION_THRESHOLD;

    _sum = Vector3f();
    _running_mean = Vector3f();
    _count = 0;

    _best_count = 0;
    _samples_since_progress = 0;

    _rotation = Matrix3f::identity();
    _tilt_deg = -1.0f;
}

void LevelCalibrator::begin(const Vector3f &expected_level_reading,
        float nominal_g, float motion_threshold, float yaw_offset_deg) {
    reset();
    _yaw_offset_rad = yaw_offset_deg * 0.01745329252f;

    // Normalised here so callers may pass any non-zero vector, and so that
    // calibrate() can assume a unit target. A degenerate one falls back to the
    // default rather than producing a rotation derived from a zero axis.
    const float len = expected_level_reading.length();
    if (len > 1e-6f) {
        _target = expected_level_reading / len;
    }
    if (nominal_g > 1e-6f) {
        _nominal_g = nominal_g;
    }
    if (motion_threshold > 0.0f) {
        _motion_threshold_sq = motion_threshold * motion_threshold;
    }
    _status = LevelCalStatus::IN_PROGRESS;
}

LevelSampleResult LevelCalibrator::addSample(const Vector3f &corrected_accel) {
    if (_status == LevelCalStatus::SUCCESS) return LevelSampleResult::REJECTED_DONE;
    if (_status != LevelCalStatus::IN_PROGRESS) {
        return LevelSampleResult::REJECTED_NOT_STARTED;
    }
    if (_count >= LEVEL_CAL_SAMPLES) return LevelSampleResult::REJECTED_DONE;

    LevelSampleResult result;
    if (_count == 0) {
        _sum = corrected_accel;
        _count = 1;
        _running_mean = corrected_accel;
        result = LevelSampleResult::ACCEPTED;
    } else {
        const Vector3f d = corrected_accel - _running_mean;
        const float dist_sq = d.x * d.x + d.y * d.y + d.z * d.z;
        if (dist_sq > _motion_threshold_sq) {
            // Same rule as the six-position gate: one excursion discards the
            // whole average rather than letting a disturbed sample into it.
            // The airframe has to be genuinely at rest for the whole run.
            _sum = Vector3f();
            _count = 0;
            result = LevelSampleResult::REJECTED_MOTION;
        } else {
            _sum = _sum + corrected_accel;
            _count++;
            _running_mean = _sum / (float)_count;
            result = (_count >= LEVEL_CAL_SAMPLES)
                    ? LevelSampleResult::ACCEPTED_DONE
                    : LevelSampleResult::ACCEPTED;
        }
    }

    // Progress is the high-water mark of the average, not its current size: a
    // motion reset drops the count to zero and that dip is not a loss of
    // progress, it is the gate doing its job. Counted for rejected samples too,
    // since an airframe that will not settle rejects every one of them.
    if (_count > _best_count) {
        _best_count = _count;
        _samples_since_progress = 0;
    } else if (_samples_since_progress < LEVEL_CAL_STALL_LIMIT) {
        _samples_since_progress++;
    }
    return result;
}

bool LevelCalibrator::isReadyToCalibrate() const {
    return _count >= LEVEL_CAL_SAMPLES;
}

float LevelCalibrator::getProgressPercent() const {
    float pct = 100.0f * (float)_count / (float)LEVEL_CAL_SAMPLES;
    return pct > 100.0f ? 100.0f : pct;
}

float LevelCalibrator::getTiltDeg() const {
    return (_status == LevelCalStatus::SUCCESS) ? _tilt_deg : -1.0f;
}

bool LevelCalibrator::rotationBetween(const Vector3f &from, const Vector3f &to,
        Matrix3f &out) {
    // Rodrigues, built from the cross and dot products directly. With v the
    // cross product, s = |v| and c the dot product:
    //     R = I + [v]x + [v]x^2 * (1-c)/s^2
    // and (1-c)/s^2 reduces to 1/(1+c), which stays well behaved right down to
    // the near-parallel case where s itself vanishes.
    const Vector3f v(from.y * to.z - from.z * to.y,
                     from.z * to.x - from.x * to.z,
                     from.x * to.y - from.y * to.x);
    const float c = from.dot(to);

    // Antiparallel: every rotation taking `from` to `to` is a half turn about
    // some axis perpendicular to both, and nothing picks between them. Refuse
    // rather than choose arbitrarily -- in this application it means the board
    // is mounted upside down, which is not something to paper over.
    if (c <= -1.0f + 1e-6f) {
        out = Matrix3f::identity();
        return false;
    }

    const float k = 1.0f / (1.0f + c);

    out.m[0][0] = c + v.x * v.x * k;
    out.m[0][1] = -v.z + v.x * v.y * k;
    out.m[0][2] = v.y + v.x * v.z * k;

    out.m[1][0] = v.z + v.y * v.x * k;
    out.m[1][1] = c + v.y * v.y * k;
    out.m[1][2] = -v.x + v.y * v.z * k;

    out.m[2][0] = -v.y + v.z * v.x * k;
    out.m[2][1] = v.x + v.z * v.y * k;
    out.m[2][2] = c + v.z * v.z * k;
    return true;
}

Matrix3f LevelCalibrator::yawRotation(float radians) {
    const float c = cosf(radians);
    const float s = sinf(radians);
    Matrix3f r = Matrix3f::identity();
    r.m[0][0] = c;    r.m[0][1] = -s;   r.m[0][2] = 0.0f;
    r.m[1][0] = s;    r.m[1][1] = c;    r.m[1][2] = 0.0f;
    r.m[2][0] = 0.0f; r.m[2][1] = 0.0f; r.m[2][2] = 1.0f;
    return r;
}

LevelCalStatus LevelCalibrator::calibrate() {
    if (_count < LEVEL_CAL_SAMPLES) {
        _status = LevelCalStatus::FAILED_NOT_ENOUGH_SAMPLES;
        return _status;
    }

    const Vector3f mean = _sum / (float)_count;
    const float mag = mean.length();

    // A reading that is not 1 g means the airframe was not at rest, or the
    // accelerometer calibration this runs on top of is wrong. Either way the
    // direction is not trustworthy enough to derive a mounting rotation from.
    if (mag < 1e-6f
            || fabsf(mag - _nominal_g) > LEVEL_CAL_MAX_MAGNITUDE_ERROR * _nominal_g) {
        _status = LevelCalStatus::FAILED_BAD_MAGNITUDE;
        return _status;
    }

    const Vector3f measured = mean / mag;

    float cos_tilt = measured.dot(_target);
    if (cos_tilt > 1.0f) cos_tilt = 1.0f;
    if (cos_tilt < -1.0f) cos_tilt = -1.0f;
    const float tilt_deg = acosf(cos_tilt) * 57.29577951f;

    if (tilt_deg > LEVEL_CAL_MAX_TILT_DEG) {
        _status = LevelCalStatus::FAILED_EXCESSIVE_TILT;
        return _status;
    }

    Matrix3f r;
    if (!rotationBetween(measured, _target, r)) {
        _status = LevelCalStatus::FAILED_EXCESSIVE_TILT;
        return _status;
    }

    // Yaw is applied after the levelling rotation, so it turns about the
    // corrected vertical rather than the crooked measured one.
    _rotation = (_yaw_offset_rad != 0.0f)
            ? yawRotation(_yaw_offset_rad).mul(r)
            : r;
    _tilt_deg = tilt_deg;
    _status = LevelCalStatus::SUCCESS;
    return _status;
}

Vector3f LevelCalibrator::correct(const Vector3f &v) const {
    if (_status != LevelCalStatus::SUCCESS) {
        return v;
    }
    return _rotation.mul(v);
}
