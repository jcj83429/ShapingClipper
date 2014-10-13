#include "ShapingClipper.h"

ShapingClipper::ShapingClipper(int sampleRate, int fftSize, int clipLevel){
  this->sampleFreq = sampleRate;
  this->size = fftSize;
  this->clipLevel = clipLevel;
  this->overlap = fftSize/4;
  this->maskSpill = fftSize/64;
  this->fft = Aquila::FftFactory::getFft(fftSize);

  this->window = new Aquila::HannWindow(fftSize);
  this->inFrame.resize(fftSize);
  this->outFrame.resize(fftSize);
  this->marginCurve.resize(fftSize/2 + 1);
  generateMarginCurve();
}

ShapingClipper::~ShapingClipper(){
  delete this->window;
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
    this->outFrame[i] = this->outFrame[i + this->overlap];
  }
  for(int i = 0; i < this->overlap; i++){
    this->inFrame[i + this->size - this->overlap] = inSamples[i];
    this->outFrame[i + this->size - this->overlap] = 0;
  }
  
  double windowedFrame[this->size], clippingDelta[this->size];

  applyWindow(this->inFrame.data(), windowedFrame);
  Aquila::SpectrumType origSpectrum = this->fft->fft(windowedFrame);
  double maskCurve[this->size/2 + 1];
  calculateMaskCurve(origSpectrum, maskCurve);

  for(int i=0; i<this->size; i++)
    clippingDelta[i] = 0;

  for(int i=0; i<4; i++){

    clipToWindow(windowedFrame, clippingDelta);  
    Aquila::SpectrumType clipSpectrum = this->fft->fft(clippingDelta);
   
    limitClipSpectrum(clipSpectrum, maskCurve);
    
    this->fft->ifft(clipSpectrum, clippingDelta);

  }

  for(int i=0; i<this->size; i++)
    windowedFrame[i] += clippingDelta[i];

  applyWindow(windowedFrame, this->outFrame.data(), true); //overlap & add
  
  for(int i = 0; i < this->overlap; i++)
    outSamples[i] = this->outFrame[i]/1.5;
    // 4 times overlap with hanning window results in 1.5 time increase in amplitude
}

void ShapingClipper::generateMarginCurve(){
  // the normal curve trashes frequencies above 16khz (because I can't hear it...but some people might)
  int points[][2] = {{0,0}, {500,40}, {2000,40}, {4000,35}, {6000,30}, {12000,30}, {16000,25}, {17000,-1000}}; //normal
  // the FM curve puts more distortion in the high frequencies to take advantage of pre/de-emphasis.
  // it also removes all distortion above 16khz as required by the FM stereo standard.
  //int points[][2] = {{0,-100}, {100,-50}, {200,40}, {3000,40}, {7000,20}, {12000,-30}, {16000,-100}, {17000,1000}}; //FM
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
}

void ShapingClipper::applyWindow(const double* inFrame, double* outFrame, const bool addToOutFrame){
  const double* window = this->window->toArray();
  for(int i = 0; i < this->size; i++){
    if(addToOutFrame)
      outFrame[i] += inFrame[i] * window[i];
    else
      outFrame[i] = inFrame[i] * window[i];
  }
}
  
void ShapingClipper::clipToWindow(const double* windowedFrame, double* clippingDelta){
  const double* window = this->window->toArray();
  for(int i = 0; i < this->size; i++){
    double limit = this->clipLevel * window[i];
    double effectiveValue = windowedFrame[i] + clippingDelta[i];
    if(effectiveValue > limit)
      clippingDelta[i] += limit - effectiveValue;
    else if(effectiveValue < -limit)
      clippingDelta[i] += -limit - effectiveValue;
  }
}

void ShapingClipper::calculateMaskCurve(const Aquila::SpectrumType &spectrum, double* maskCurve){
  for(int j = 0; j < this->size / 2 + 1; j++)
    maskCurve[j] = 0;

  for(int j = 0; j < this->maskSpill; j++)
    maskCurve[0+j] += abs(spectrum[0]) / (j*128/this->size + 1);
  for(int i = 1; i < this->size / 2; i++){
    // upward spill
    for(int j = 0; j < this->maskSpill; j++){
      int idx = i+j;
      idx = (idx > this->size / 2 ? this->size / 2 : idx);
      maskCurve[idx] += (abs(spectrum[i]) + abs(spectrum[this->size - i])) / (j*128/this->size + 1);
    }
    // downward spill
    for(int j = 1; j < this->maskSpill / 2; j++){
      int idx = i-j;
      idx = (idx < 0 ? 0 : idx);
      maskCurve[idx] += (abs(spectrum[i]) + abs(spectrum[this->size - i])) / (j*256/this->size + 1);
    }
  }
  maskCurve[this->size / 2] += abs(spectrum[this->size / 2]);
}

void ShapingClipper::limitClipSpectrum(Aquila::SpectrumType &clipSpectrum, const double* maskCurve){
  double* marginCurve = this->marginCurve.data(); // margin curve is already in dB
  double relativeDistortionLevel = Aquila::dB(abs(clipSpectrum[0])) - (Aquila::dB(maskCurve[0]) - marginCurve[0]);
  if(relativeDistortionLevel > 0)
    clipSpectrum[0] *= pow(10, -relativeDistortionLevel / 20);
  for(int i = 1; i < this->size / 2; i++){
    relativeDistortionLevel = Aquila::dB(abs(clipSpectrum[i]) + abs(clipSpectrum[this->size - i])) - (Aquila::dB(maskCurve[i]) - marginCurve[i]);
    if(relativeDistortionLevel > 0){
      clipSpectrum[i] *= pow(10, -relativeDistortionLevel / 20);
      clipSpectrum[this->size - i] *= pow(10, -relativeDistortionLevel / 20);
    }
  }
  relativeDistortionLevel = Aquila::dB(abs(clipSpectrum[this->size / 2])) - (Aquila::dB(maskCurve[this->size / 2]) - marginCurve[this->size / 2]);
  if(relativeDistortionLevel > 0)
    clipSpectrum[this->size / 2] *= pow(10, -relativeDistortionLevel / 20);
}
