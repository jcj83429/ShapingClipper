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
  void feed(const float* inSamples, float* outSamples, bool diffOnly=false);

  /**
   *  Returns fftSize*3/4
   */
  int getDelay();

  /**
   *  Returns fftSize/4
   */
  int getFeedSize();

  /**
   *  Set the clipping level in sample value
   */
  void setClipLevel(float clipLevel);

  /**
   *  Set the number of clipping iterations.
   *  Setting iterations to 0 effectively works as a bypass.
   */
  void setIterations(int iterations);

  /**
   * Set the adaptive distortion strength.
   * The adaptive distortion strength affects how easily the clipper gives up
   * distortion control to reach the clip level target.
   */
  void setAdaptiveDistortionStrength(float strength);

 private:
  int size;
  int overlap;
  int numPsyBins;
  PFFFT_Setup *pffft;
  float sampleFreq;
  float clipLevel;
  float iterations;
  float adaptiveDistortionStrength;

  /**
   *  inFrame: unmodified input audio
   *  outDistFrame: clipping distortion multiplied by 1.5. The 1.5 factor is due to overlap and add.
   *  marginCurve: see generateMarginCurve
   *  window: the Hann window
   *  invWindow: inverse of the hanning window used by limitPeak
   *  spreadTable: a (size/2 + 1) by (size/2 + 1) matrix containing each bin's contribution to each bin's calculated masking level
   */
  std::vector<float> inFrame, outDistFrame, marginCurve, window, invWindow, spreadTable;

  /**
   *  Generate the Hann window and inverse window.
   */
  void generateHannWindow();

   /**
   *  Generate the spreadTable used by calculateMaskCurve
   */
  void generateSpreadTable();

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
   *  maskCurve is in linear scale.
   */
  void calculateMaskCurve(const float *spectrum, float* maskCurve);

  /**
   *  Based on the maskCurve and the marginCurve, limit the amount of
   *  distortion.
   */
  void limitClipSpectrum(float *clipSpectrum, const float* maskCurve);

};
