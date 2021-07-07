#include "ShapingClipper.h"
#include <algorithm>
#include <complex>

ShapingClipper::ShapingClipper(int sampleRate, int fftSize, int clipLevel){
  this->sampleFreq = sampleRate;
  this->size = fftSize;
  this->clipLevel = clipLevel;
  this->overlap = fftSize/4;
  this->maskSpill = fftSize/64;
  this->pffft = pffft_new_setup(fftSize, PFFFT_REAL);

  this->window.resize(fftSize);
  this->invWindow.resize(fftSize);
  generateHannWindow();

  this->inFrame.resize(fftSize);
  this->outDistFrame.resize(fftSize);
  this->marginCurve.resize(fftSize/2 + 1);

  this->windowedFrame = new float[fftSize];
  this->clippingDelta = new float[fftSize];
  this->maskCurve = new float[fftSize/2 + 1];
  this->spectrumBuf = new float[fftSize];

  generateMarginCurve();
}

ShapingClipper::~ShapingClipper(){
  delete[] this->windowedFrame;
  delete[] this->clippingDelta;
  delete[] this->maskCurve;
  delete[] this->spectrumBuf;
  pffft_destroy_setup(this->pffft);
}

int ShapingClipper::getFeedSize(){
  return this->overlap;
}

int ShapingClipper::getDelay(){
  return this->size - this->overlap;
}

void ShapingClipper::feed(const double* inSamples, double* outSamples){
  //shift in/out buffers
  for(int i = 0; i < this->size - this->overlap; i++){
    this->inFrame[i] = this->inFrame[i + this->overlap];
    this->outDistFrame[i] = this->outDistFrame[i + this->overlap];
  }
  for(int i = 0; i < this->overlap; i++){
    this->inFrame[i + this->size - this->overlap] = inSamples[i];
    this->outDistFrame[i + this->size - this->overlap] = 0;
  }
  
  double peak;

  applyWindow(this->inFrame.data(), windowedFrame);
  pffft_transform_ordered(this->pffft, windowedFrame, spectrumBuf, NULL, PFFFT_FORWARD);
  calculateMaskCurve(spectrumBuf, maskCurve);

  //clear clippingDelta
  for(int i=0; i<this->size; i++)
    clippingDelta[i] = 0;


  //repeat clipping process a few times to get more clipping
  for(int i=0; i<6; i++){
    clipToWindow(windowedFrame, clippingDelta, (i >= 4 ? 2.0 : 1.0)); //last 2 rounds have boosted delta
    pffft_transform_ordered(this->pffft, clippingDelta, spectrumBuf, NULL, PFFFT_FORWARD);

    limitClipSpectrum(spectrumBuf, maskCurve);
    
    pffft_transform_ordered(this->pffft, spectrumBuf, clippingDelta, NULL, PFFFT_BACKWARD);
    // see pffft.h
    for(int i = 0; i < this->size; i++) clippingDelta[i] /= this->size;

    peak = 0;
    for(int i=0; i<this->size; i++)
      peak = std::max<float>(peak, std::abs((windowedFrame[i] + clippingDelta[i]) * invWindow[i]));

    float maskCurveShift = std::max<float>(peak, 1.122); // 1.122 is 1dB

    //be less strict in the next iteration
    for(int i = 0; i < this->size / 2 + 1; i++)
      maskCurve[i] *= maskCurveShift;
  }

  applyWindow(clippingDelta, this->outDistFrame.data(), true); //overlap & add
  
  for(int i = 0; i < this->overlap; i++)
    outSamples[i] = this->inFrame[i] + this->outDistFrame[i]/1.5;
    // 4 times overlap with squared hanning window results in 1.5 time increase in amplitude
}

void ShapingClipper::generateHannWindow() {
    double pi = acos(-1);
    for (int i = 0; i < this->size; i++) {
        float value = 0.5 * (1 - cos(2 * pi * i / this->size));
        this->window[i] = value;
        // 1/window to calculate unwindowed peak.
        this->invWindow[i] = value > 0.1 ? 1.0 / (value * clipLevel) : 0;
    }
}

void ShapingClipper::generateMarginCurve(){
  // the normal curve trashes frequencies above 16khz (because I can't hear it...but some people might)
  int points[][2] = {{0,-10}, {80,0}, {200,20}, {1000,20}, {6000,25}, {10000,25}, {16000,20}, {20000,0}}; //normal

  // the FM curve puts more distortion in the high frequencies to take advantage of pre/de-emphasis.
  // it also removes all distortion above 16khz as required by the FM stereo standard.
  // The FM curve may be outdated as it has not been tested with the current distortion shaping logic.
  //int points[][2] = {{0,-100}, {100,-20}, {200,0}, {1000,20}, {4000,20}, {10000,5}, {16000,-5}, {17000,1000}}; //FM

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

void ShapingClipper::calculateMaskCurve(const float *spectrum, float* maskCurve){
  const double maskSpillBaseVal = (this->size * 44100.0) / (256 * this->sampleFreq); //(size/256) * (48000/sampleFreq)
  //const int maskSpillBaseVal = 1 * (this->size / 256);
  const int maskSpillBaseFreq = 2000; //maskSpill = 1 at 2000 Hz
  float amp;
  int nextMaskSpillBand = maskSpillBaseFreq * this->size / this->sampleFreq;
  int maskSpill = maskSpillBaseVal;
  int maskSpillScale = 1;

  for(int j = 0; j < this->size / 2 + 1; j++)
    maskCurve[j] = 0;

  // handle bin 0 (DC)
  for(int j = 0; j < this->size / 64; j++)
    maskCurve[0+j] += std::abs(spectrum[0]) / (j*128/this->size + 1);

  for(int i = 1; i < this->size / 2; i++){
    // although the negative frequencies are omitted because they are redundant,
    // the magnitude of the positive frequencies are not doubled.
    // Multiply the magnitude by 2 to simulate adding up the + and - frequencies.
    amp = abs(std::complex<float>(spectrum[2*i], spectrum[2*i + 1])) * 2;
    maskCurve[i] += amp;

      // upward spill
      for(int j = 1; j < maskSpill; j++){
	int idx = i+j;
	if(idx <= this->size / 2)
	  maskCurve[idx] += amp / (4*j/maskSpill + 1);
      }
      // downward spill
      for(int j = 1; j < maskSpill / 2; j++){
	int idx = i-j;
	if(idx >= 0)
	  maskCurve[idx] += amp / (8*j/maskSpill + 1);
      }

      if(i >= nextMaskSpillBand){
	maskSpillScale++;
	maskSpill = maskSpillScale * maskSpillBaseVal;
	nextMaskSpillBand = maskSpillScale * maskSpillBaseFreq * this->size / this->sampleFreq;
      }
  }
  // top bin (nyquist frequency) is stored in array element 1
  maskCurve[this->size / 2] += std::abs(spectrum[1]);

  for(int i = 0; i < this->size / 2 + 1; i++)
    maskCurve[i] = maskCurve[i] / this->marginCurve[i];
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
