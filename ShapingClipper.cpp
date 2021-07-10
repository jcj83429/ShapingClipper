#include "ShapingClipper.h"
#include <algorithm>
#include <complex>

ShapingClipper::ShapingClipper(int sampleRate, int fftSize, float clipLevel){
  this->sampleFreq = sampleRate;
  this->size = fftSize;
  this->clipLevel = clipLevel;
  this->iterations = 6;
  this->adaptiveDistortionStrength = 1.0;
  this->overlap = fftSize/4;
  this->maskSpill = fftSize/64;
  this->pffft = pffft_new_setup(fftSize, PFFFT_REAL);

  this->window.resize(fftSize);
  this->invWindow.resize(fftSize);
  generateHannWindow();

  this->inFrame.resize(fftSize);
  this->outDistFrame.resize(fftSize);
  this->marginCurve.resize(fftSize/2 + 1);
  this->spreadTable.resize(((int64_t)fftSize / 2 + 1) * ((int64_t)fftSize / 2 + 1));

  generateMarginCurve();
  generateSpreadTable();
}

ShapingClipper::~ShapingClipper(){
  pffft_destroy_setup(this->pffft);
}

int ShapingClipper::getFeedSize(){
  return this->overlap;
}

int ShapingClipper::getDelay(){
  return this->size - this->overlap;
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

void ShapingClipper::feed(const float* inSamples, float* outSamples, bool diffOnly){
  //shift in/out buffers
  for(int i = 0; i < this->size - this->overlap; i++){
    this->inFrame[i] = this->inFrame[i + this->overlap];
    this->outDistFrame[i] = this->outDistFrame[i + this->overlap];
  }
  for(int i = 0; i < this->overlap; i++){
    this->inFrame[i + this->size - this->overlap] = inSamples[i];
    this->outDistFrame[i + this->size - this->overlap] = 0;
  }
  
  float peak;
  float *windowedFrame = (float*) alloca(sizeof(float) * this->size);
  float *clippingDelta = (float*) alloca(sizeof(float) * this->size);
  float *spectrumBuf = (float*) alloca(sizeof(float) * this->size);
  float *maskCurve = (float*) alloca(sizeof(float) * this->size/2 + 1);

  applyWindow(this->inFrame.data(), windowedFrame);
  pffft_transform_ordered(this->pffft, windowedFrame, spectrumBuf, NULL, PFFFT_FORWARD);
  calculateMaskCurve(spectrumBuf, maskCurve);

  //clear clippingDelta
  for(int i=0; i<this->size; i++)
    clippingDelta[i] = 0;


  //repeat clipping process a few times to get more clipping
  for(int i = 0; i < this->iterations; i++){
    // the last 1/3 rounds have boosted delta
    float deltaBoost = i >= this->iterations - this->iterations / 3 ? 2.0 : 1.0;
    clipToWindow(windowedFrame, clippingDelta, deltaBoost);
    pffft_transform_ordered(this->pffft, clippingDelta, spectrumBuf, NULL, PFFFT_FORWARD);

    limitClipSpectrum(spectrumBuf, maskCurve);
    
    pffft_transform_ordered(this->pffft, spectrumBuf, clippingDelta, NULL, PFFFT_BACKWARD);
    // see pffft.h
    for(int i = 0; i < this->size; i++) clippingDelta[i] /= this->size;

    peak = 0;
    for(int i=0; i<this->size; i++)
      peak = std::max<float>(peak, std::abs((windowedFrame[i] + clippingDelta[i]) * invWindow[i]));
    peak /= this->clipLevel;

    float maskCurveShift = std::max<float>(peak, 1.122); // 1.122 is 1dB
    maskCurveShift = 1.0 + (maskCurveShift - 1.0) * this->adaptiveDistortionStrength;

    //be less strict in the next iteration
    for(int i = 0; i < this->size / 2 + 1; i++)
      maskCurve[i] *= maskCurveShift;
  }

  applyWindow(clippingDelta, this->outDistFrame.data(), true); //overlap & add
  
  for(int i = 0; i < this->overlap; i++) {
    outSamples[i] = this->outDistFrame[i]/1.5;
    // 4 times overlap with squared hanning window results in 1.5 time increase in amplitude
    if(!diffOnly) {
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
    for (int i = 0; i < this->size / 2 + 1; i++) {
        float sum = 0;
        int baseIdx = i * (this->size / 2 + 1);
        // Calculate tent-shape function in log-log scale.
        // As an optimization, only consider bins close to i
        // This reduces the number of multiplications needed in calculateMaskCurve
        // The masking contribution at faraway bins is negligeable
        int startBin = i * 3 / 4;
        int endBin = std::min(this->size + 1, ((i + 1) * 4 + 2) / 3);
        for (int j = startBin; j < endBin; j++) {
            // add 0.5 so i=0 doesn't get log(0)
            float relIdxLog = std::abs(log((j + 0.5) / (i + 0.5)));
            float value;
            if (j >= i) {
                // mask up
                value = exp(-relIdxLog * 40);
            } else {
                // mask down
                value = exp(-relIdxLog * 80);
            }
            sum += value;
            this->spreadTable[baseIdx + j] = value;
        }
        // now normalize it
        for (int j = startBin; j < endBin; j++) {
            this->spreadTable[baseIdx + j] /= sum;
            //printf("%f ", spreadingFunction[baseIdx + j]);
        }
        //printf("\n\n");
    }
}

void ShapingClipper::generateMarginCurve(){
  // the normal curve trashes frequencies above 16khz (because I can't hear it...but some people might)
  int points[][2] = { {0,-10}, {80,-10}, {200,10}, {1000,20}, {4000,20}, {8000,15}, {16000,5}, {20000,-10} }; //normal

  int numPoints = 8;
  this->marginCurve[0] = points[0][1];
  
  int j = 0;
  for(int i = 0; i < numPoints-1; i++){
    while(j < this->size / 2 + 1 && j * this->sampleFreq / this->size < points[i+1][0]){
      marginCurve[j] = points[i][1] + (j - points[i][0] * this->size / this->sampleFreq) * (points[i+1][1] - points[i][1]) / ((points[i+1][0] - points[i][0]) * this->size / this->sampleFreq); // linear interpolation
      j++;
    }
  }
  while(j < this->size / 2 + 1){
    marginCurve[j] = points[numPoints-1][1];
    j++;
  }

  // convert margin curve to linear amplitude scale
  for(j = 0; j < this->size / 2 + 1; j++){
    marginCurve[j] = pow(10, marginCurve[j] / 20);
  }
}

void ShapingClipper::applyWindow(const float* inFrame, float* outFrame, const bool addToOutFrame){
  const float* window = this->window.data();
  for(int i = 0; i < this->size; i++){
    if(addToOutFrame)
      outFrame[i] += inFrame[i] * window[i];
    else
      outFrame[i] = inFrame[i] * window[i];
  }
}
  
void ShapingClipper::clipToWindow(const float* windowedFrame, float* clippingDelta, float deltaBoost){
  const float* window = this->window.data();
  for(int i = 0; i < this->size; i++){
    float limit = this->clipLevel * window[i];
    float effectiveValue = windowedFrame[i] + clippingDelta[i];
    if(effectiveValue > limit)
      clippingDelta[i] += (limit - effectiveValue)*deltaBoost;
    else if(effectiveValue < -limit)
      clippingDelta[i] += (-limit - effectiveValue)*deltaBoost;
  }
}

void ShapingClipper::calculateMaskCurve(const float* spectrum, float* maskCurve) {
    for (int i = 0; i < this->size / 2 + 1; i++) {
        maskCurve[i] = 0;
    }
    for (int i = 0; i < this->size / 2 + 1; i++) {
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
        int baseIdx = i * (this->size / 2 + 1);
        // As an optimization, only consider bins close to i
        // This reduces the number of multiplications needed in calculateMaskCurve
        // The masking contribution at faraway bins is negligeable
        int startBin = i * 3 / 4;
        int endBin = std::min(this->size + 1, ((i + 1) * 4 + 2) / 3);
        for (int j = startBin; j < endBin; j++) {
            maskCurve[j] += this->spreadTable[baseIdx + j] * magnitude;
        }
    }
    for (int i = 0; i < this->size / 2 + 1; i++) {
        maskCurve[i] = maskCurve[i] / this->marginCurve[i];
    }
}

void ShapingClipper::limitClipSpectrum(float *clipSpectrum, const float* maskCurve){
  // bin 0
  float relativeDistortionLevel = std::abs(clipSpectrum[0]) / maskCurve[0];
  if(relativeDistortionLevel > 1.0)
    clipSpectrum[0] /= relativeDistortionLevel;
  // bin 1..N/2-1
  for(int i = 1; i < this->size / 2; i++){
    float real = clipSpectrum[i*2];
    float imag = clipSpectrum[i*2 + 1];
    // although the negative frequencies are omitted because they are redundant,
    // the magnitude of the positive frequencies are not doubled.
    // Multiply the magnitude by 2 to simulate adding up the + and - frequencies.
    relativeDistortionLevel = abs(std::complex<float>(real, imag)) * 2 / maskCurve[i];
    if(relativeDistortionLevel > 1.0){
      clipSpectrum[i*2] /= relativeDistortionLevel;
      clipSpectrum[i*2 + 1] /= relativeDistortionLevel;
    }
  }
  // bin N/2
  relativeDistortionLevel = std::abs(clipSpectrum[1]) / maskCurve[this->size / 2];
  if(relativeDistortionLevel > 1.0)
    clipSpectrum[1] /= relativeDistortionLevel;
}
