#include "shaping_clipper.h"
#include <algorithm>
#include <complex>
#include <cmath>

// Enable BIN_GAIN_DEBUG to output white noise shaped by the bin gain.
// It can be viewed with a spectrogram.
#define BIN_GAIN_DEBUG 0

// Run calculate_mask_curve again after applying bin_gain to calculate bin_level_in
// instead of scaling the mask_curve bin by bin. This should be more accurate.
#define REAL_BIN_LEVEL_IN 1

shaping_clipper::shaping_clipper(int sample_rate, int fft_size, float clip_level, int max_oversample) {
    this->sample_rate = sample_rate;
    this->size = fft_size;
    max_oversample = 1 << (int)log2(std::max(1, max_oversample));
    this->max_oversample = max_oversample;
    this->oversample = 1;
    this->clip_level = clip_level;
    this->iterations = 6;
    this->adaptive_distortion_strength = 1.0;
    this->overlap = fft_size / 4;
    this->pffft = pffft_new_setup(fft_size, PFFFT_REAL);
    this->pffft_oversampled = NULL;

    // The psy masking calculation is O(n^2),
    // so skip it for frequencies not covered by base sampling rantes (i.e. 44k)
    if (sample_rate <= 50000) {
        this->num_psy_bins = fft_size / 2;
    } else if (sample_rate <= 100000) {
        this->num_psy_bins = fft_size / 4;
    } else {
        this->num_psy_bins = fft_size / 8;
    }

    this->window.resize(fft_size * max_oversample);
    this->inv_window.resize(fft_size * max_oversample);
    generate_hann_window();

    this->in_frame.resize(fft_size);
    this->out_dist_frame.resize(fft_size);
    this->margin_curve.resize(fft_size / 2 + 1);
    this->bin_gain.resize(fft_size / 2 + 1);
    this->bin_hold.resize(fft_size / 2 + 1);
    for (int i = 0; i < fft_size / 2 + 1; i++) {
        this->bin_gain[i] = 1.0;
        this->bin_hold[i] = 0;
    }
    // normally I use __builtin_ctz for this but the intrinsic is different for
    // different compilers.
    // 2 spread_table entries per octave.
    int spread_table_rows = log2(this->num_psy_bins) * 2;
    this->spread_table.resize(spread_table_rows * this->num_psy_bins);
    this->spread_table_range.resize(spread_table_rows);
    this->spread_table_index.resize(this->num_psy_bins);

    // default curve
    int points[][2] = { {0,15}, {125,15}, {250,15}, {500,16}, {1000,16}, {2000,16}, {4000,15}, {8000,15}, {16000,15}, {20000,14} };
    int num_points = 10;
    set_margin_curve(points, num_points);
    generate_spread_table();
    set_compress_speed(300, 300);
}

shaping_clipper::~shaping_clipper() {
    pffft_destroy_setup(this->pffft);
    if (this->pffft_oversampled) {
        pffft_destroy_setup(this->pffft_oversampled);
    }
}

int shaping_clipper::get_feed_size() {
    return this->overlap;
}

void shaping_clipper::set_clip_level(float clip_level) {
    this->clip_level = clip_level;
}

void shaping_clipper::set_iterations(int iterations) {
    this->iterations = iterations;
}

void shaping_clipper::set_adaptive_distortion_strength(float strength) {
    this->adaptive_distortion_strength = strength;
}

void shaping_clipper::set_oversample(int oversample) {
    oversample = std::max(1, std::min(this->max_oversample, oversample));
    oversample = 1 << (int)log2(oversample);
    if (oversample != this->oversample) {
        if (this->pffft_oversampled) {
            pffft_destroy_setup(this->pffft_oversampled);
            this->pffft_oversampled = NULL;
        }
        if (oversample > 1) {
            this->pffft_oversampled = pffft_new_setup(this->size * oversample, PFFFT_REAL);
        }
    }
    this->oversample = oversample;
}

void shaping_clipper::set_compress_speed(float attack_db_per_sec, float release_db_per_sec) {
    float attack_db_per_frame = attack_db_per_sec * this->overlap / this->sample_rate;
    attack_speed = pow(10, -attack_db_per_frame / 20);
    float release_db_per_frame = release_db_per_sec * this->overlap / this->sample_rate;
    release_speed = pow(10, release_db_per_frame / 20);
}

void shaping_clipper::feed(const float* in_samples, float* out_samples, bool diff_only, float* total_margin_shift) {
    // shift in/out buffers
    for (int i = 0; i < this->size - this->overlap; i++) {
        this->in_frame[i] = this->in_frame[i + this->overlap];
        this->out_dist_frame[i] = this->out_dist_frame[i + this->overlap];
    }
    for (int i = 0; i < this->overlap; i++) {
        this->in_frame[i + this->size - this->overlap] = in_samples[i];
        this->out_dist_frame[i + this->size - this->overlap] = 0;
    }

    float peak;
    float* windowed_frame = (float*)alloca(sizeof(float) * this->size * this->oversample);
    float* clipping_delta = (float*)alloca(sizeof(float) * this->size * this->oversample);
    float* spectrum_buf = (float*)alloca(sizeof(float) * this->size * this->oversample);
    float* spectrum_for_atten = (float*)alloca(sizeof(float) * this->size);
    float* mask_curve = (float*)alloca(sizeof(float) * (this->size / 2 + 1));
    float* mask_curve2 = (float*)alloca(sizeof(float) * (this->size / 2 + 1));
    float* bin_level_in = (float*)alloca(sizeof(float) * (this->size / 2 + 1));

    apply_window(this->in_frame.data(), windowed_frame);

    pffft_transform_ordered(this->pffft, windowed_frame, spectrum_buf, NULL, PFFFT_FORWARD);

    calculate_mask_curve(spectrum_buf, mask_curve);

#if REAL_BIN_LEVEL_IN
    // borrow the spectrum_for_atten buffer.
    for (int i = 0; i < this->num_psy_bins; i++) {
        spectrum_for_atten[i * 2] = spectrum_buf[i * 2] * bin_gain[i];
        spectrum_for_atten[i * 2 + 1] = spectrum_buf[i * 2 + 1] * bin_gain[i];
    }
    calculate_mask_curve(spectrum_for_atten, bin_level_in);
#else
    // Reuse the spreading from the mask curve for the bin_level_in.
    // This is not completely accurate because it's not the same as applying the gain then doing the spreading.
    for (int i = 0; i < this->num_psy_bins; i++) {
        bin_level_in[i] = mask_curve[i] * bin_gain[i];
    }
#endif

    mask_curve2[0] = abs(spectrum_buf[0]) * (1.0 - bin_gain[0]);
    for (int i = 1; i < this->num_psy_bins; i++) {
        float real = spectrum_buf[i * 2];
        float imag = spectrum_buf[i * 2 + 1];
        // although the negative frequencies are omitted because they are redundant,
        // the magnitude of the positive frequencies are not doubled.
        // Multiply the magnitude by 2 to simulate adding up the + and - frequencies.
        mask_curve2[i] = abs(std::complex<float>(real, imag)) * 2 * (1.0 - bin_gain[i]);
    }

    // high sample rates: clear mask_curve2 for ultrasonic frequencies
    for (int i = this->num_psy_bins; i < this->size / 2 + 1; i++) {
        mask_curve2[i] = 0;
    }


    int clipping_samples = this->size * this->oversample;
    int window_stride = this->max_oversample / this->oversample;
    PFFFT_Setup* clipping_pffft = this->pffft;
    if (this->oversample > 1) {
        clipping_pffft = this->pffft_oversampled;

        // use IFFT to oversample the windowed frame
        spectrum_buf[this->size] = spectrum_buf[1] / 2;
        spectrum_buf[1] = 0;
        for (int i = this->size + 1; i < this->size * this->oversample; i++) {
            spectrum_buf[i] = 0;
        }
        // adjust scaling for oversampling
        for (int i = 0; i <= this->size; i++) {
            spectrum_buf[i] *= this->oversample;
        }
        for (int i = 0; i < this->size / 2 + 1; i++) {
            mask_curve[i] *= this->oversample;
            mask_curve2[i] *= this->oversample;
        }
        pffft_transform_ordered(this->pffft_oversampled, spectrum_buf, windowed_frame, NULL, PFFFT_BACKWARD);
        for (int i = 0; i < clipping_samples; i++) {
            windowed_frame[i] /= clipping_samples;
        }
    }

    // apply bin gain through clipping delta
    spectrum_buf[0] *= bin_gain[0] - 1.0;
    spectrum_buf[1] = 0;
    for (int i = 1; i < this->num_psy_bins; i++) {
        spectrum_buf[i * 2] *= bin_gain[i] - 1.0;
        spectrum_buf[i * 2 + 1] *= bin_gain[i] - 1.0;
    }
    // if oversampling is active, zero out bins above audible range. bin_gain doesn't apply there.
    for (int i = this->num_psy_bins * 2; i < clipping_samples; i++) {
        spectrum_buf[i] = 0;
    }
    memcpy(spectrum_for_atten, spectrum_buf, sizeof(float) * this->size);
    pffft_transform_ordered(clipping_pffft, spectrum_buf, clipping_delta, NULL, PFFFT_BACKWARD);
    for (int i = 0; i < clipping_samples; i++) {
        clipping_delta[i] /= clipping_samples;
    }

    float orig_peak = 0;
    for (int sample_idx = 0, window_idx = 0; sample_idx < clipping_samples; sample_idx++, window_idx += window_stride) {
        orig_peak = std::max<float>(orig_peak, std::abs(windowed_frame[sample_idx] * inv_window[window_idx]));
    }
    orig_peak /= this->clip_level;
    peak = orig_peak;

    if (total_margin_shift) {
        *total_margin_shift = 1.0;
    }

    // repeat clipping-filtering process a few times to control both the peaks and the spectrum
    for (int i = 0; i < this->iterations; i++) {
        // The last 1/3 of rounds have boosted delta to help reach the peak target faster
        float delta_boost = 1.0;
        if (i >= this->iterations - this->iterations / 3) {
            // boosting the delta when largs peaks are still present is dangerous
            if (peak < 2.0) {
                delta_boost = 2.0;
            }
        }

        for (int i = 0; i < this->size / 2 + 1; i++) {
            mask_curve2[i] = std::max(mask_curve[i], mask_curve2[i]);
        }

        clip_to_window(windowed_frame, clipping_delta, delta_boost);

        pffft_transform_ordered(clipping_pffft, clipping_delta, spectrum_buf, NULL, PFFFT_FORWARD);

        if (this->oversample > 1) {
            // Zero out all frequency bins above the base nyquist rate.
            // limit_clip_spectrum doesn't handle these bins.
            spectrum_buf[1] = 0;
            for (int i = this->size + 1; i < this->size * this->oversample; i++) {
                spectrum_buf[i] = 0;
            }
        }

        limit_clip_spectrum(spectrum_buf, mask_curve2);

        // For bins whose allowed clipping distortion is less than the effect of attenuation hold,
        // forbid all clipping distortion by setting their values back to the attenuation hold values.
        // For bins that are allowed to have more distortion, do distortion control on the clipping
        // distortion on top of the attenuation. Distortion popping in and out can happen otherwise.

        // All these special cases for bin 0 and N are ugly.
        // I should allocate one more slot in the spectrum_buf and move bin N to the end.
        // Also, normalize the scaling for bin 0 and N vs the rest so no x2 is needed for middle bins.
        if (mask_curve[0] < mask_curve2[0]) {
            spectrum_buf[0] = spectrum_for_atten[0];
        } else {
            float atten_real = spectrum_for_atten[0];
            float atten_mag = abs(atten_real);
            float remaining_mag = mask_curve2[0] - atten_mag;
            float clip_real = spectrum_buf[0] - spectrum_for_atten[0];
            float clip_mag = abs(clip_real);
            if (clip_mag > remaining_mag) {
                float scale = remaining_mag / clip_mag;
                spectrum_buf[0] = atten_real + clip_real * scale;
            }
        }
        for (int i = 1; i < this->num_psy_bins; i++) {
            if (mask_curve[i] < mask_curve2[i]) {
                spectrum_buf[i * 2] = spectrum_for_atten[i * 2];
                spectrum_buf[i * 2 + 1] = spectrum_for_atten[i * 2 + 1];
            } else {
                // FIXME: Calculating the atten_mag every time is inefficient.
                // We had it in mask_curve2 before the clipping-filtering loop.
                float atten_real = spectrum_for_atten[i * 2];
                float atten_imag = spectrum_for_atten[i * 2 + 1];
                float atten_mag = abs(std::complex<float>(atten_real, atten_imag)) * 2;
                float remaining_mag = mask_curve2[i] - atten_mag;
                float clip_real = spectrum_buf[i * 2] - atten_real;
                float clip_imag = spectrum_buf[i * 2 + 1] - atten_imag;
                float clip_mag = abs(std::complex<float>(clip_real, clip_imag)) * 2;
                if (clip_mag > remaining_mag) {
                    float scale = remaining_mag / clip_mag;
                    spectrum_buf[i * 2] = atten_real + clip_real * scale;
                    spectrum_buf[i * 2 + 1] = atten_imag + clip_imag * scale;
                }
            }
        }

        pffft_transform_ordered(clipping_pffft, spectrum_buf, clipping_delta, NULL, PFFFT_BACKWARD);
        // see pffft.h
        for (int i = 0; i < clipping_samples; i++) {
            clipping_delta[i] /= clipping_samples;
        }

        peak = 0;
        for (int sample_idx = 0, window_idx = 0; sample_idx < clipping_samples; sample_idx++, window_idx += window_stride) {
            peak = std::max<float>(peak, std::abs((windowed_frame[sample_idx] + clipping_delta[sample_idx]) * inv_window[window_idx]));
        }
        peak /= this->clip_level;

        // Automatically adjust mask_curve as necessary to reach peak target
        float mask_curve_shift = 1.122; // 1.122 is 1dB
        if (orig_peak > 1.0 && peak > 1.0) {
            float diff_achieved = orig_peak - peak;
            if (i + 1 < this->iterations - this->iterations / 3 && diff_achieved > 0) {
                float diff_needed = orig_peak - 1.0;
                float diff_ratio = diff_needed / diff_achieved;
                // If a good amount of peak reduction was already achieved,
                // don't shift the mask_curve by the full peak value
                // On the other hand, if only a little peak reduction was achieved,
                // don't shift the mask_curve by the enormous diff_ratio.
                diff_ratio = std::min<float>(diff_ratio, peak);
                mask_curve_shift = std::max<float>(mask_curve_shift, diff_ratio);
            } else {
                // If the peak got higher than the input or we are in the last 1/3 rounds,
                // go back to the heavy-handed peak heuristic.
                mask_curve_shift = std::max<float>(mask_curve_shift, peak);
            }
        }

        mask_curve_shift = 1.0 + (mask_curve_shift - 1.0) * this->adaptive_distortion_strength;

        if (total_margin_shift && peak > 1.01 && i < this->iterations - 1) {
            *total_margin_shift *= mask_curve_shift;
        }

        // Be less strict in the next iteration.
        // This helps with peak control.
        for (int i = 0; i < this->size / 2 + 1; i++) {
            mask_curve[i] *= mask_curve_shift;
        }
    }

    if (this->oversample > 1) {
        // Downsample back to original rate.
        // We already zeroed out all frequency bins above the base nyquist rate so we can simply drop samples here.
        for (int i = 0; i < this->size; i++) {
            clipping_delta[i] = clipping_delta[i * this->oversample];
            windowed_frame[i] = windowed_frame[i * this->oversample];
        }
    }

#if BIN_GAIN_DEBUG
    // output white noise scaled by bin_gain
    float* debug_temp = (float*)alloca(sizeof(float) * this->size);
    spectrum_buf[0] = bin_gain[0] * (float)rand() * 256 / RAND_MAX;
    for (int i = 0; i < this->num_psy_bins; i++) {
        float gain = bin_gain[i];
        spectrum_buf[i * 2] = gain * (float)rand() * 256 / RAND_MAX;
        spectrum_buf[i * 2 + 1] = gain * (float)rand() * 256 / RAND_MAX;
    }
    pffft_transform_ordered(this->pffft, spectrum_buf, debug_temp, NULL, PFFFT_BACKWARD);
    apply_window(debug_temp, this->out_dist_frame.data(), true);
    diff_only = true;
#else

    // do overlap & add
    apply_window(clipping_delta, this->out_dist_frame.data(), true);

#endif

    for (int i = 0; i < this->overlap; i++) {
        // 4 times overlap with squared hanning window results in 1.5 time increase in amplitude
        out_samples[i] = this->out_dist_frame[i] / 1.5;
    }

    if (this->oversample > 1) {
        // When oversampling is active, the clipping delta can be non-zero even if all the samples in the frame are well below the threshold.
        // This is likely due to the aliasing introduced when windowing near-nyquist frequencies.
        // For purity, detect and remove these unwanted changes to the input signal.
        // The effect on peak control is < 0.1% or -60dB. Oversampling is not exact anyway.
        float max_out_diff = 0;
        for (int i = 0; i < this->overlap; i++) {
            max_out_diff = std::max(max_out_diff, std::abs(out_samples[i]));
        }
        if (max_out_diff < this->clip_level * 0.001) {
            for (int i = 0; i < this->overlap; i++) {
                out_samples[i] = 0;
            }
        }
    }

    if (!diff_only) {
        for (int i = 0; i < this->overlap; i++) {
            out_samples[i] += this->in_frame[i];
        }
    }

    // update bin gain
    for (int i = 0; i < this->size; i++) {
        windowed_frame[i] += clipping_delta[i];
    }
    pffft_transform_ordered(this->pffft, windowed_frame, spectrum_buf, NULL, PFFFT_FORWARD);
    // reuse the spreading of the mask curve to calculate final bin level
    calculate_mask_curve(spectrum_buf, mask_curve);
    update_bin_gain(bin_level_in, mask_curve);
}

void shaping_clipper::generate_hann_window() {
    double pi = acos(-1);
    int window_size = this->size * this->max_oversample;
    for (int i = 0; i < window_size; i++) {
        float value = 0.5 * (1 - cos(2 * pi * i / window_size));
        this->window[i] = value;
        // 1/window to calculate unwindowed peak.
        this->inv_window[i] = value > 0.1 ? 1.0 / value : 0;
    }
}

void shaping_clipper::generate_spread_table() {
    // Calculate tent-shape function in log-log scale.

    // As an optimization, only consider bins close to "bin"
    // This reduces the number of multiplications needed in calculate_mask_curve
    // The masking contribution at faraway bins is negligeable

    // Another optimization to save memory and speed up the calculation of the
    // spread table is to calculate and store only 2 spread functions per
    // octave, and reuse the same spread function for multiple bins.
    int table_index = 0;
    int bin = 0;
    int increment = 1;
    while (bin < this->num_psy_bins) {
        float freq = bin * this->sample_rate / this->size;
        float sum = 0;
        int base_idx = table_index * this->num_psy_bins;
        int start_bin = bin * 3 / 4;
        int end_bin = std::min(this->num_psy_bins, ((bin + 1) * 4 + 2) / 3);

        for (int j = start_bin; j < end_bin; j++) {
            // add 0.5 so i=0 doesn't get log(0)
            float rel_idx_log = std::abs(log((j + 0.5) / (bin + 0.5)));
            float value;
            if (j >= bin) {
                // mask up
                value = exp(-rel_idx_log * 60);
            } else {
                // mask down
                value = exp(-rel_idx_log * 100);
            }
            // the spreading function is centred in the row
            sum += value;
            this->spread_table[base_idx + this->num_psy_bins / 2 + j - bin] = value;
        }
        // now normalize it
        for (int j = start_bin; j < end_bin; j++) {
            this->spread_table[base_idx + this->num_psy_bins / 2 + j - bin] /= sum;
        }

        this->spread_table_range[table_index] = std::make_pair(start_bin - bin, end_bin - bin);

        int next_bin;
        if (bin <= 1) {
            next_bin = bin + 1;
        } else {
            if ((bin & (bin - 1)) == 0) {
                // power of 2
                increment = bin / 2;
            }
            next_bin = bin + increment;
        }

        // set bins between "bin" and "next_bin" to use this table_index
        for (int i = bin; i < next_bin; i++) {
            this->spread_table_index[i] = table_index;
        }

        bin = next_bin;
        table_index++;
    }
}

void shaping_clipper::set_margin_curve(int points[][2], int num_points) {
    this->margin_curve[0] = points[0][1];

    int j = 0;
    for (int i = 0; i < num_points - 1; i++) {
        while (j < this->size / 2 + 1 && j * this->sample_rate / this->size < points[i + 1][0]) {
            // linearly interpolate between points
            float binHz = j * this->sample_rate / this->size;
            margin_curve[j] = points[i][1] + (binHz - points[i][0]) * (points[i + 1][1] - points[i][1]) / (points[i + 1][0] - points[i][0]);
            j++;
        }
    }
    // handle bins after the last point
    while (j < this->size / 2 + 1) {
        margin_curve[j] = points[num_points - 1][1];
        j++;
    }

    // convert margin curve to linear amplitude scale
    for (j = 0; j < this->size / 2 + 1; j++) {
        margin_curve[j] = pow(10, margin_curve[j] / 20);
    }
}

void shaping_clipper::apply_window(const float* in_frame, float* out_frame, const bool add_to_out_frame) {
    const float* window = this->window.data();
    int total_samples = this->size;
    int window_stride = this->max_oversample;
    for (int i = 0; i < total_samples; i++) {
        if (add_to_out_frame) {
            out_frame[i] += in_frame[i] * *window;
        } else {
            out_frame[i] = in_frame[i] * *window;
        }
        window += window_stride;
    }
}

void shaping_clipper::clip_to_window(const float* windowed_frame, float* clipping_delta, float delta_boost) {
    const float* window = this->window.data();
    int total_samples = this->size * this->oversample;
    int window_stride = this->max_oversample / this->oversample;
    for (int i = 0; i < total_samples; i++) {
        float limit = this->clip_level * *window;
        float effective_value = windowed_frame[i] + clipping_delta[i];
        if (effective_value > limit) {
            clipping_delta[i] += (limit - effective_value) * delta_boost;
        } else if (effective_value < -limit) {
            clipping_delta[i] += (-limit - effective_value) * delta_boost;
        }
        window += window_stride;
    }
}

void shaping_clipper::calculate_mask_curve(const float* spectrum, float* mask_curve) {
    for (int i = 0; i < this->size / 2 + 1; i++) {
        mask_curve[i] = 0;
    }
    for (int i = 0; i < this->num_psy_bins; i++) {
        float magnitude;
        if (i == 0) {
            magnitude = std::abs(spectrum[0]);
        } else if (i == this->size / 2) {
            magnitude = std::abs(spectrum[1]);
        } else {
            // although the negative frequencies are omitted because they are redundant,
            // the magnitude of the positive frequencies are not doubled.
            // Multiply the magnitude by 2 to simulate adding up the + and - frequencies.
            magnitude = abs(std::complex<float>(spectrum[2 * i], spectrum[2 * i + 1])) * 2;
        }
        float energy = magnitude * magnitude;
        int table_idx = this->spread_table_index[i];
        std::pair<int, int> range = this->spread_table_range[table_idx];
        int base_idx = table_idx * this->num_psy_bins;
        int start_bin = std::max(0, i + range.first);
        int end_bin = std::min(this->num_psy_bins, i + range.second);
        for (int j = start_bin; j < end_bin; j++) {
            mask_curve[j] += this->spread_table[base_idx + this->num_psy_bins / 2 + j - i] * energy;
        }
    }
    // sqrt back to magnitude
    for (int i = 0; i < this->num_psy_bins; i++) {
        mask_curve[i] = sqrt(mask_curve[i]);
    }

    // for ultrasonic frequencies, skip the O(n^2) spread calculation and just copy the magnitude
    for (int i = this->num_psy_bins; i < this->size / 2 + 1; i++) {
        float magnitude;
        if (i == this->size / 2) {
            magnitude = std::abs(spectrum[1]);
        } else {
            // although the negative frequencies are omitted because they are redundant,
            // the magnitude of the positive frequencies are not doubled.
            // Multiply the magnitude by 2 to simulate adding up the + and - frequencies.
            magnitude = abs(std::complex<float>(spectrum[2 * i], spectrum[2 * i + 1])) * 2;
        }
        mask_curve[i] = magnitude;
    }

    for (int i = 0; i < this->size / 2 + 1; i++) {
        mask_curve[i] = mask_curve[i] / this->margin_curve[i];
    }
}

void shaping_clipper::limit_clip_spectrum(float* clip_spectrum, const float* mask_curve) {
    // bin 0
    float relative_distortion_level = std::abs(clip_spectrum[0]) / mask_curve[0];
    if (relative_distortion_level > 1.0) {
        clip_spectrum[0] /= relative_distortion_level;
    }
    // bin 1..N/2-1
    // When oversampling is on, the base nyquist bin is handled in this loop
    int i = 1;
    int bins_below_nyquist = std::min(this->size * this->oversample / 2, this->size / 2 + 1);
    for (; i < bins_below_nyquist; i++) {
        float real = clip_spectrum[i * 2];
        float imag = clip_spectrum[i * 2 + 1];
        // although the negative frequencies are omitted because they are redundant,
        // the magnitude of the positive frequencies are not doubled.
        // Multiply the magnitude by 2 to simulate adding up the + and - frequencies.
        relative_distortion_level = abs(std::complex<float>(real, imag)) * 2 / mask_curve[i];
        if (relative_distortion_level > 1.0) {
            clip_spectrum[i * 2] /= relative_distortion_level;
            clip_spectrum[i * 2 + 1] /= relative_distortion_level;
        }
    }
    // bin N/2
    // When oversampling is off, the base nyquist bins needs to be handled here.
    if (i == this->size / 2) {
        relative_distortion_level = std::abs(clip_spectrum[1]) / mask_curve[this->size / 2];
        if (relative_distortion_level > 1.0) {
            clip_spectrum[1] /= relative_distortion_level;
        }
    }
}

void shaping_clipper::update_bin_gain(const float* bin_level_in, const float* bin_level_out) {
    float min_bin_gain = 1.0;
    for (int i = 0; i < this->num_psy_bins; i++) {
        float bin_atten = 1.0;
        if (bin_level_in[i]) {
            bin_atten = std::min(1.0f, bin_level_out[i] / bin_level_in[i]);
        }
        // Due to the shortcut used to calculate bin_level_in by scaling mask_curve,
        // there can be small differences between mask_curve and bin_level_in even when there is no attenuation,
        // so use 0.99 here to ignore the small differences.
        if (bin_atten < 0.99) {
            bin_atten = std::max(bin_atten, attack_speed);
            bin_gain[i] *= bin_atten;
            bin_hold[i] = 2;
        } else {
            if (bin_hold[i]) {
                bin_hold[i]--;
            } else {
                bin_gain[i] *= release_speed;
            }
        }
        bin_gain[i] = std::min<float>(bin_gain[i], 1.0);
    }
    // Limit gain difference between bins to avoid extreme EQ changes.
    for (int i = 0; i < this->num_psy_bins; i++) {
        bin_gain[i] = std::min(bin_gain[i], min_bin_gain * 4);
    }
}

void shaping_clipper_stereo_coupler::sync_bin_gain() {
    float max_diff = 1.122;
    for (int i = 0; i < m_clipper_l->num_psy_bins; i++) {
        float &gain_l = m_clipper_l->bin_gain[i];
        float &gain_r = m_clipper_r->bin_gain[i];
        if (gain_l < gain_r) {
            gain_r = std::min(gain_l * max_diff, gain_r);
        } else {
            gain_l = std::min(gain_r * max_diff, gain_l);
        }
    }
}