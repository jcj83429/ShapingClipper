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

  //clear clippingDelta
  for(int i=0; i<this->size; i++)
    clippingDelta[i] = 0;

  //repeat clipping process a few times to get more clipping
  for(int i=0; i<4; i++){

    clipToWindow(windowedFrame, clippingDelta);
    Aquila::SpectrumType clipSpectrum = this->fft->fft(clippingDelta);
   
    limitClipSpectrum(clipSpectrum, maskCurve);
    
    this->fft->ifft(clipSpectrum, clippingDelta);

    //be less strict in the next iteration
    for(int i = 0; i < this->size / 2 + 1; i++)
      this->marginCurve[i] -= 2;
  }

  for(int i = 0; i < this->size / 2 + 1; i++)
    this->marginCurve[i] += 8;

  for(int i=0; i<this->size; i++)
    windowedFrame[i] += clippingDelta[i];

  //limitPeak(windowedFrame);

  applyWindow(windowedFrame, this->outFrame.data(), true); //overlap & add
  
  for(int i = 0; i < this->overlap; i++)
    outSamples[i] = this->outFrame[i]/1.5;
    // 4 times overlap with hanning window results in 1.5 time increase in amplitude
}

void ShapingClipper::generateMarginCurve(){
  // the normal curve trashes frequencies above 16khz (because I can't hear it...but some people might)
  int points[][2] = {{0,-10}, {80,0}, {200,20}, {1000,15}, {6000,20}, {10000,25}, {16000,20}, {17000,-1000}}; //normal

  // the FM curve puts more distortion in the high frequencies to take advantage of pre/de-emphasis.
  // it also removes all distortion above 16khz as required by the FM stereo standard.
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
  
void ShapingClipper::clipToWindow(const double* windowedFrame, double* clippingDelta, double deltaBoost){
  const double* window = this->window->toArray();
  for(int i = 0; i < this->size; i++){
    double limit = this->clipLevel * window[i];
    double effectiveValue = windowedFrame[i] + clippingDelta[i];
    if(effectiveValue > limit)
      clippingDelta[i] += (limit - effectiveValue)*deltaBoost;
    else if(effectiveValue < -limit)
      clippingDelta[i] += (-limit - effectiveValue)*deltaBoost;
  }
}

void ShapingClipper::calculateMaskCurve(const Aquila::SpectrumType &spectrum, double* maskCurve){
  const double maskSpillBaseVal = (this->size * 44100.0) / (256 * this->sampleFreq); //(size/256) * (48000/sampleFreq)
  //const int maskSpillBaseVal = 1 * (this->size / 256);
  const int maskSpillBaseFreq = 2000; //maskSpill = 1 at 2000 Hz
  double amp;
  int nextMaskSpillBand = maskSpillBaseFreq * this->size / this->sampleFreq;
  int maskSpill = maskSpillBaseVal;
  int maskSpillScale = 1;

  for(int j = 0; j < this->size / 2 + 1; j++)
    maskCurve[j] = 0;

  for(int j = 0; j < this->size / 64; j++)
    maskCurve[0+j] += abs(spectrum[0]) / (j*128/this->size + 1);

  for(int i = 1; i < this->size / 2; i++){
    amp = abs(spectrum[i]) + abs(spectrum[this->size - i]);
    //amp = amp / maskSpill;
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
  maskCurve[this->size / 2] += abs(spectrum[this->size / 2]);
}

void ShapingClipper::limitClipSpectrum(Aquila::SpectrumType &clipSpectrum, const double* maskCurve){
  double* marginCurve = this->marginCurve.data(); // margin curve is already in dB
  double relativeDistortionLevel = Aquila::dB(abs(clipSpectrum[0]) / maskCurve[0]) + marginCurve[0];
  if(relativeDistortionLevel > 0)
    clipSpectrum[0] *= pow(10, -relativeDistortionLevel / 20);
  for(int i = 1; i < this->size / 2; i++){
    relativeDistortionLevel = Aquila::dB((abs(clipSpectrum[i]) + abs(clipSpectrum[this->size - i])) / maskCurve[i]) + marginCurve[i];
    if(relativeDistortionLevel > 0){
      clipSpectrum[i] *= pow(10, -relativeDistortionLevel / 20);
      clipSpectrum[this->size - i] *= pow(10, -relativeDistortionLevel / 20);
    }
  }
  relativeDistortionLevel = Aquila::dB(abs(clipSpectrum[this->size / 2]) / maskCurve[this->size / 2]) + marginCurve[this->size / 2];
  if(relativeDistortionLevel > 0)
    clipSpectrum[this->size / 2] *= pow(10, -relativeDistortionLevel / 20);
}

void ShapingClipper::limitPeak(double* windowedFrame){
  const double* window = this->window->toArray();
  double multiplier = 1;
  for(int i = this->size / 4; i < this->size * 3 / 4; i++){
    double absVal = abs(windowedFrame[i]) / this->clipLevel;
    if(absVal > window[i]){
      double newMult = window[i]/absVal;
      if(newMult < multiplier){
	multiplier = newMult;
      }
    }
  }

  if(multiplier < 0.9999)
    for(int i = 0; i < this->size; i++)
      windowedFrame[i] *= multiplier;
}
