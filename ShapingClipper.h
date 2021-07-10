/**
 *  An audio clipper effect that limits the amount of distortion
 */

#include <vector>
#include "pffft.h"

class ShapingClipper
{
 public:
  
  /**
   *  sampleRate is used for psychoacoustic purposes only.
   *  fftSize should be a multiple of 4.
   *  clipLevel means symmetric clipping from -clipLevel to +clipLevel
   */
  ShapingClipper(int sampleRate, int fftSize, float clipLevel=16384);
  ~ShapingClipper();
  
  /**
   *  Put fftSize/4 samples in inSamples and get fftSize/4 samples out in
   *  outSamples.
   *  The output is delayed by fftSize*3/4 samples.
   */
  void feed(const float* inSamples, float* outSamples);

  /**
   *  Returns fftSize*3/4
   */
  int getDelay();

  /**
   *  Returns fftSize/4
   */
  int getFeedSize();
  float clipLevel;

 private:
  int size;
  int overlap;
  int maskSpill;
  PFFFT_Setup *pffft;
  float sampleFreq;

  /**
   *  inFrame: unmodified input audio
   *  outDistFrame: clipping distortion multiplied by 1.5. The 1.5 factor is due to overlap and add.
   *  marginCurve: see generateMarginCurve
   *  invWindow: inverse of the hanning window used by limitPeak
   */
  std::vector<float> inFrame, outDistFrame, marginCurve, window, invWindow;

  // these are buffers used by feed()
  float *windowedFrame, *clippingDelta, *maskCurve, *spectrumBuf;

  /**
   *  Generate the Hann window and inverse window.
   */
  void generateHannWindow();

  /**
   *  Generates a basic psychoacoustic curve.
   *  marginCurve represents the minimum ratio between
   *  the clean input and the clipping distortion at each frequency
   */
  void generateMarginCurve();

  /**
   *  Applies the window to the inFrame and store the result in the outFrame
   *  If addToOutFrame is true, the results is added to the outFrame instead
   */
  void applyWindow(const float* inFrame, float* outFrame, const bool addToOutFrame=false);

  /**
   *  Clips the windowedFrame to the window scaled by the clipLevel.
   *  The clipping distortion is multiplied by deltaBoost
   *  The existing values in clippingDelta is applied to the windowedFrame
   *  to get the effective sample values, taking previous clipping iterations
   *  into account.
   *  Should only be used with windowed input
   */
  void clipToWindow(const float* windowedFrame, float* clippingDelta, float deltaBoost=1.0);

  /**
   *  Calculates the original signal level considering psychoacoustic masking.
   *  Masking width scales with frequency.
   *  maskCurve is logarithmic.
   */
  void calculateMaskCurve(const float *spectrum, float* maskCurve);

  /**
   *  Based on the maskCurve and the marginCurve, limit the amount of
   *  distortion.
   */
  void limitClipSpectrum(float *clipSpectrum, const float* maskCurve);

};
