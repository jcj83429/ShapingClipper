#include "ShapingClipper.h"
#include <algorithm>
#include <complex>
#include <cmath>

ShapingClipper::ShapingClipper(int sampleRate, int fftSize, float clipLevel) {
    this->sampleFreq = sampleRate;
    this->size = fftSize;
    this->clipLevel = clipLevel;
    this->iterations = 6;
    this->adaptiveDistortionStrength = 1.0;
    this->overlap = fftSize / 4;
    this->pffft = pffft_new_setup(fftSize, PFFFT_REAL);

    // The psy masking calculation is O(n^2),
    // so skip it for frequencies not covered by base sampling rantes (i.e. 44k)
    if (sampleRate <= 50000) {
        this->numPsyBins = fftSize / 2;
    } else if (sampleRate <= 100000) {
        this->numPsyBins = fftSize / 4;
    } else {
        this->numPsyBins = fftSize / 8;
    }

    this->window.resize(fftSize);
    this->invWindow.resize(fftSize);
    generateHannWindow();

    this->inFrame.resize(fftSize);
    this->outDistFrame.resize(fftSize);
    this->marginCurve.resize(fftSize / 2 + 1);
    // normally I use __builtin_ctz for this but the intrinsic is different for
    // different compilers.
    // 2 spreadTable entries per octave.
    int spreadTableRows = log2(this->numPsyBins) * 2;
    this->spreadTable.resize(spreadTableRows * this->numPsyBins);
    this->spreadTableRange.resize(spreadTableRows);
    this->spreadTableIndex.resize(this->numPsyBins);

    // default curve
    int points[][2] = { {0,14}, {125,14}, {250,16}, {500,18}, {1000,20}, {2000,20}, {4000,20}, {8000,15}, {16000,5}, {20000,-10} };
    int numPoints = 10;
    setMarginCurve(points, numPoints);
    generateSpreadTable();
}

ShapingClipper::~ShapingClipper() {
    pffft_destroy_setup(this->pffft);
}

int ShapingClipper::getFeedSize() {
    return this->overlap;
}

void ShapingClipper::setClipLevel(float clipLevel) {
    this->clipLevel = clipLevel;
}

void ShapingClipper::setIterations(int iterations) {
    this->iterations = iterations;
}

void ShapingClipper::setAdaptiveDistortionStrength(float strength) {
    this->adaptiveDistortionStrength = strength;
}

void ShapingClipper::feed(const float* inSamples, float* outSamples, bool diffOnly, float* totalMarginShift) {
    // shift in/out buffers
    for (int i = 0; i < this->size - this->overlap; i++) {
        this->inFrame[i] = this->inFrame[i + this->overlap];
        this->outDistFrame[i] = this->outDistFrame[i + this->overlap];
    }
    for (int i = 0; i < this->overlap; i++) {
        this->inFrame[i + this->size - this->overlap] = inSamples[i];
        this->outDistFrame[i + this->size - this->overlap] = 0;
    }

    float peak;
    float* windowedFrame = (float*)alloca(sizeof(float) * this->size);
    float* clippingDelta = (float*)alloca(sizeof(float) * this->size);
    float* spectrumBuf = (float*)alloca(sizeof(float) * this->size);
    float* maskCurve = (float*)alloca(sizeof(float) * this->size / 2 + 1);

    applyWindow(this->inFrame.data(), windowedFrame);
    pffft_transform_ordered(this->pffft, windowedFrame, spectrumBuf, NULL, PFFFT_FORWARD);
    calculateMaskCurve(spectrumBuf, maskCurve);

    // It would be easier to calculate the peak from the unwindowed input.
    // This is just for consistency with the clipped peak calculateion
    // because the invWindow zeros out samples on the edge of the window.
    float origPeak = 0;
    for (int i = 0; i < this->size; i++) {
        origPeak = std::max<float>(origPeak, std::abs(windowedFrame[i] * invWindow[i]));
    }
    origPeak /= this->clipLevel;

    // clear clippingDelta
    for (int i = 0; i < this->size; i++) {
        clippingDelta[i] = 0;
    }

    if (totalMarginShift) {
        *totalMarginShift = 1.0;
    }

    // repeat clipping-filtering process a few times to control both the peaks and the spectrum
    for (int i = 0; i < this->iterations; i++) {
        // The last 1/3 of rounds have boosted delta to help reach the peak target faster
        float deltaBoost = i >= this->iterations - this->iterations / 3 ? 2.0 : 1.0;
        clipToWindow(windowedFrame, clippingDelta, deltaBoost);

        pffft_transform_ordered(this->pffft, clippingDelta, spectrumBuf, NULL, PFFFT_FORWARD);

        limitClipSpectrum(spectrumBuf, maskCurve);

        pffft_transform_ordered(this->pffft, spectrumBuf, clippingDelta, NULL, PFFFT_BACKWARD);
        // see pffft.h
        for (int i = 0; i < this->size; i++) {
            clippingDelta[i] /= this->size;
        }

        peak = 0;
        for (int i = 0; i < this->size; i++)
            peak = std::max<float>(peak, std::abs((windowedFrame[i] + clippingDelta[i]) * invWindow[i]));
        peak /= this->clipLevel;

        // Automatically adjust maskCurve as necessary to reach peak target
        float maskCurveShift = 1.122; // 1.122 is 1dB
        if (origPeak > 1.0 && peak > 1.0) {
            float diffAchieved = origPeak - peak;
            if (i + 1 < this->iterations - this->iterations / 3 && diffAchieved > 0) {
                float diffNeeded = origPeak - 1.0;
                float diffRatio = diffNeeded / diffAchieved;
                // If a good amount of peak reduction was already achieved,
                // don't shift the maskCurve by the full peak value
                // On the other hand, if only a little peak reduction was achieved,
                // don't shift the maskCurve by the enormous diffRatio.
                diffRatio = std::min<float>(diffRatio, peak);
                maskCurveShift = std::max<float>(maskCurveShift, diffRatio);
            } else {
                // If the peak got higher than the input or we are in the last 1/3 rounds,
                // go back to the heavy-handed peak heuristic.
                maskCurveShift = std::max<float>(maskCurveShift, peak);
            }
        }

        maskCurveShift = 1.0 + (maskCurveShift - 1.0) * this->adaptiveDistortionStrength;

        if (totalMarginShift && peak > 1.01 && i < this->iterations - 1) {
            *totalMarginShift *= maskCurveShift;
        }

        // Be less strict in the next iteration.
        // This helps with peak control.
        for (int i = 0; i < this->size / 2 + 1; i++) {
            maskCurve[i] *= maskCurveShift;
        }
    }

    // do overlap & add
    applyWindow(clippingDelta, this->outDistFrame.data(), true);

    for (int i = 0; i < this->overlap; i++) {
        outSamples[i] = this->outDistFrame[i] / 1.5;
        // 4 times overlap with squared hanning window results in 1.5 time increase in amplitude
        if (!diffOnly) {
            outSamples[i] += this->inFrame[i];
        }
    }
}

void ShapingClipper::generateHannWindow() {
    double pi = acos(-1);
    for (int i = 0; i < this->size; i++) {
        float value = 0.5 * (1 - cos(2 * pi * i / this->size));
        this->window[i] = value;
        // 1/window to calculate unwindowed peak.
        this->invWindow[i] = value > 0.1 ? 1.0 / value : 0;
    }
}

void ShapingClipper::generateSpreadTable() {
    // Calculate tent-shape function in log-log scale.

    // As an optimization, only consider bins close to "bin"
    // This reduces the number of multiplications needed in calculateMaskCurve
    // The masking contribution at faraway bins is negligeable

    // Another optimization to save memory and speed up the calculation of the
    // spread table is to calculate and store only 2 spread functions per
    // octave, and reuse the same spread function for multiple bins.
    int tableIndex = 0;
    int bin = 0;
    int increment = 1;
    while (bin < this->numPsyBins) {
        float sum = 0;
        int baseIdx = tableIndex * this->numPsyBins;
        int startBin = bin * 3 / 4;
        int endBin = std::min(this->numPsyBins, ((bin + 1) * 4 + 2) / 3);

        for (int j = startBin; j < endBin; j++) {
            // add 0.5 so i=0 doesn't get log(0)
            float relIdxLog = std::abs(log((j + 0.5) / (bin + 0.5)));
            float value;
            if (j >= bin) {
                // mask up
                value = exp(-relIdxLog * 40);
            } else {
                // mask down
                value = exp(-relIdxLog * 80);
            }
            // the spreading function is centred in the row
            sum += value;
            this->spreadTable[baseIdx + this->numPsyBins / 2 + j - bin] = value;
        }
        // now normalize it
        for (int j = startBin; j < endBin; j++) {
            this->spreadTable[baseIdx + this->numPsyBins / 2 + j - bin] /= sum;
        }

        this->spreadTableRange[tableIndex] = std::make_pair(startBin - bin, endBin - bin);

        int nextBin;
        if (bin <= 1) {
            nextBin = bin + 1;
        } else {
            if ((bin & (bin - 1)) == 0) {
                // power of 2
                increment = bin / 2;
            }
            nextBin = bin + increment;
        }

        // set bins between "bin" and "nextBin" to use this tableIndex
        for (int i = bin; i < nextBin; i++) {
            this->spreadTableIndex[i] = tableIndex;
        }

        bin = nextBin;
        tableIndex++;
    }
}

void ShapingClipper::setMarginCurve(int points[][2], int numPoints) {
    this->marginCurve[0] = points[0][1];

    int j = 0;
    for (int i = 0; i < numPoints - 1; i++) {
        while (j < this->size / 2 + 1 && j * this->sampleFreq / this->size < points[i + 1][0]) {
            // linearly interpolate between points
            marginCurve[j] = points[i][1] + (j - points[i][0] * this->size / this->sampleFreq) * (points[i + 1][1] - points[i][1]) / ((points[i + 1][0] - points[i][0]) * this->size / this->sampleFreq);
            j++;
        }
    }
    // handle bins after the last point
    while (j < this->size / 2 + 1) {
        marginCurve[j] = points[numPoints - 1][1];
        j++;
    }

    // convert margin curve to linear amplitude scale
    for (j = 0; j < this->size / 2 + 1; j++) {
        marginCurve[j] = pow(10, marginCurve[j] / 20);
    }
}

void ShapingClipper::applyWindow(const float* inFrame, float* outFrame, const bool addToOutFrame) {
    const float* window = this->window.data();
    for (int i = 0; i < this->size; i++) {
        if (addToOutFrame) {
            outFrame[i] += inFrame[i] * window[i];
        } else {
            outFrame[i] = inFrame[i] * window[i];
        }
    }
}

void ShapingClipper::clipToWindow(const float* windowedFrame, float* clippingDelta, float deltaBoost) {
    const float* window = this->window.data();
    for (int i = 0; i < this->size; i++) {
        float limit = this->clipLevel * window[i];
        float effectiveValue = windowedFrame[i] + clippingDelta[i];
        if (effectiveValue > limit) {
            clippingDelta[i] += (limit - effectiveValue) * deltaBoost;
        } else if (effectiveValue < -limit) {
            clippingDelta[i] += (-limit - effectiveValue) * deltaBoost;
        }
    }
}

void ShapingClipper::calculateMaskCurve(const float* spectrum, float* maskCurve) {
    for (int i = 0; i < this->size / 2 + 1; i++) {
        maskCurve[i] = 0;
    }
    for (int i = 0; i < this->numPsyBins; i++) {
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
        int tableIdx = this->spreadTableIndex[i];
        std::pair<int, int> range = this->spreadTableRange[tableIdx];
        int baseIdx = tableIdx * this->numPsyBins;
        int startBin = std::max(0, i + range.first);
        int endBin = std::min(this->numPsyBins, i + range.second);
        for (int j = startBin; j < endBin; j++) {
            maskCurve[j] += this->spreadTable[baseIdx + this->numPsyBins / 2 + j - i] * magnitude;
        }
    }

    // for ultrasonic frequencies, skip the O(n^2) spread calculation and just copy the magnitude
    for (int i = this->numPsyBins; i < this->size / 2 + 1; i++) {
        float magnitude;
        if (i == this->size / 2) {
            magnitude = std::abs(spectrum[1]);
        } else {
            // although the negative frequencies are omitted because they are redundant,
            // the magnitude of the positive frequencies are not doubled.
            // Multiply the magnitude by 2 to simulate adding up the + and - frequencies.
            magnitude = abs(std::complex<float>(spectrum[2 * i], spectrum[2 * i + 1])) * 2;
        }
        maskCurve[i] = magnitude;
    }

    for (int i = 0; i < this->size / 2 + 1; i++) {
        maskCurve[i] = maskCurve[i] / this->marginCurve[i];
    }
}

void ShapingClipper::limitClipSpectrum(float* clipSpectrum, const float* maskCurve) {
    // bin 0
    float relativeDistortionLevel = std::abs(clipSpectrum[0]) / maskCurve[0];
    if (relativeDistortionLevel > 1.0) {
        clipSpectrum[0] /= relativeDistortionLevel;
    }
    // bin 1..N/2-1
    for (int i = 1; i < this->size / 2; i++) {
        float real = clipSpectrum[i * 2];
        float imag = clipSpectrum[i * 2 + 1];
        // although the negative frequencies are omitted because they are redundant,
        // the magnitude of the positive frequencies are not doubled.
        // Multiply the magnitude by 2 to simulate adding up the + and - frequencies.
        relativeDistortionLevel = abs(std::complex<float>(real, imag)) * 2 / maskCurve[i];
        if (relativeDistortionLevel > 1.0) {
            clipSpectrum[i * 2] /= relativeDistortionLevel;
            clipSpectrum[i * 2 + 1] /= relativeDistortionLevel;
        }
    }
    // bin N/2
    relativeDistortionLevel = std::abs(clipSpectrum[1]) / maskCurve[this->size / 2];
    if (relativeDistortionLevel > 1.0) {
        clipSpectrum[1] /= relativeDistortionLevel;
    }
}
