/**
 *  An audio clipper effect that limits the amount of distortion
 */

#include <vector>
#include "aquila/global.h"
#include "aquila/source/window/HannWindow.h"
#include "aquila/transform/FftFactory.h"
#include "aquila/functions.h"

class ShapingClipper
{
 public:
  
  /**
   *  sampleRate is used for psychoacoustic purposes only.
   *  fftSize should be a multiple of 4.
   *  clipLevel means symmetric clipping from -clipLevel to +clipLevel
   */
  ShapingClipper(int sampleRate, int fftSize, int clipLevel=16384);
  ~ShapingClipper();
  
  /**
   *  Put fftSize/4 samples in inSamples and get fftSize/4 samples out in
   *  outSamples.
   *  The output is delayed by fftSize*3/4 samples.
   */
  void feed(const double* inSamples, double* outSamples);

  /**
   *  Returns fftSize*3/4
   */
  int getDelay();

  /**
   *  Returns fftSize/4
   */
  int getFeedSize();
  int clipLevel;

 private:
  int size;
  int overlap;
  int maskSpill;
  std::shared_ptr<Aquila::Fft> fft;
  Aquila::FrequencyType sampleFreq;

  Aquila::HannWindow* window;
  std::vector<double> inFrame, outFrame, marginCurve, invWindow;

  /**
   *  Generates a basic psychoacoustic curve.
   *  marginCurve is in dB and represents the minimum level difference between
   *  the clean input and the clipping distortion at each frequency
   */
  void generateMarginCurve();

  /**
   *  Applies the window to the inFrame and store the result in the outFrame
   *  If addToOutFrame is true, the results is added to the outFrame instead
   */
  void applyWindow(const double* inFrame, double* outFrame, const bool addToOutFrame=false);

  /**
   *  Clips the windowedFrame to the window scaled by the clipLevel.
   *  The clipping distortion is multiplied by deltaBoost
   *  The existing values in clippingDelta is applied to the windowedFrame
   *  to get the effective sample values, taking previous clipping iterations
   *  into account.
   *  Should only be used with windowed input
   */
  void clipToWindow(const double* windowedFrame, double* clippingDelta, double deltaBoost=1.0);

  /**
   *  Calculates the original signal level considering psychoacoustic masking.
   *  Masking width scales with frequency.
   *  maskCurve is logarithmic.
   */
  void calculateMaskCurve(const Aquila::SpectrumType &spectrum, double* maskCurve);

  /**
   *  Based on the maskCurve and the marginCurve, limit the amount of
   *  distortion.
   */
  void limitClipSpectrum(Aquila::SpectrumType &clipSpectrum, const double* maskCurve);

  /**
   *  Scale the frame to the clipping threshold.
   *  Not sure if this will be useful.
   */
  void limitPeak(double* windowedFrame);
};
