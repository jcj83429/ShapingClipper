/**
 *  An audio clipper effect that limits the amount of distortion
 */

#include <vector>
#include "pffft.h"

class ShapingClipper
{
public:

    /**
     *  sampleRate is only6 used to generate the marginCurve.
     *  fftSize should be a multiple of 4.
     *  clipLevel means symmetric clipping from -clipLevel to +clipLevel
     */
    ShapingClipper(int sampleRate, int fftSize, float clipLevel = 16384);
    ~ShapingClipper();

    /**
     *  Put fftSize/4 samples in inSamples and get fftSize/4 samples out in outSamples.
     *  The output in outSamples corresponds to the input 3 calls ago.
     *  totalMarginShift is an optional parameter that can be used to get the total margin adjustment in the frame
     */
    void feed(const float* inSamples, float* outSamples, bool diffOnly = false, float* totalMarginShift = NULL);

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
     *  Set the adaptive distortion strength.
     *  The adaptive distortion strength affects how easily the clipper gives up
     *  distortion control to reach the clip level target.
     */
    void setAdaptiveDistortionStrength(float strength);

    /**
     *  Sets the marginCurve from a list of (Hz, dB) tuples.
     *  The curve is linearly interpolated between the points
     *  The first point must be at 0Hz.
     *  marginCurve represents the minimum ratio between
     *  the clean input and the clipping distortion at each frequency
     */
    void setMarginCurve(int points[][2], int numPoints);

private:
    int size;
    int overlap;
    int numPsyBins;
    PFFFT_Setup* pffft;
    float sampleFreq;
    float clipLevel;
    float iterations;
    float adaptiveDistortionStrength;

    /**
     *  inFrame: unmodified input audio
     *  outDistFrame: clipping distortion multiplied by 1.5. The 1.5 factor is due to overlap and add.
     *  marginCurve: see generateMarginCurve
     *  window: the Hann window
     *  invWindow: inverse of the Hann window used to calculate the unwindowed peak
     *  spreadTable: see generateSpreadTable
     */
    std::vector<float> inFrame, outDistFrame, marginCurve, window, invWindow, spreadTable;

    /**
     *  spreadtableIndex: for each bin, which entry in the spread table to use
     *                    (eacn entry is numPsyBins values)
     *  spreadTableRange: for each entry in the spread table, what is the range of bins (ie. -2 to +4) that are non-zero
     */
    std::vector<int> spreadTableIndex;
    std::vector<std::pair<int, int>> spreadTableRange;

    /**
     *  Generate the Hann window and inverse window.
     */
    void generateHannWindow();

    /**
    *  Generate the spreading functions used by calculateMaskCurve
    *  The spreadTable contains entries of size/2 values each.
    *  To save memory, only 2 entries are stored per octave, and each entry is shared by a range of bins.
    *  The spread scale proportionally with frequency.
    *  The spreadTable entries are normalized to add up to 1.
    #  Eacn entry is centred, meaning the numPsyBins/2'th value is the peak of the tent-shaped function.
    */
    void generateSpreadTable();

    /**
     *  Applies the window to the inFrame and store the result in the outFrame
     *  If addToOutFrame is true, the results is added to the outFrame instead
     */
    void applyWindow(const float* inFrame, float* outFrame, const bool addToOutFrame = false);

    /**
     *  Clips the windowedFrame to the window scaled by the clipLevel.
     *  The clipping distortion is multiplied by deltaBoost
     *  The existing values in clippingDelta is applied to the windowedFrame
     *  to get the effective sample values, taking previous clipping iterations
     *  into account.
     *  Should only be used with windowed input
     */
    void clipToWindow(const float* windowedFrame, float* clippingDelta, float deltaBoost = 1.0);

    /**
     *  Calculates the original signal level considering psychoacoustic masking.
     *  maskCurve is in linear scale.
     */
    void calculateMaskCurve(const float* spectrum, float* maskCurve);

    /**
     *  Limit the magnitude of each bin to the maskCurve
     */
    void limitClipSpectrum(float* clipSpectrum, const float* maskCurve);

};
