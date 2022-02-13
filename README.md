ShapingClipper
==============
A distortion-controlling audio clipper that uses FFT to shape the clipping distortion.

This effect is available as the Psychoacoustic Clipper in Calf Studio Gear or apsyclip in libavfilter (ffmpeg, mpv, etc)

## How it works
It divides the input into overlapping blocks. For each block, it does the following.
1. Calculate the input spectrum and the masking threshold ("maskCurve") based on the input spectrum.
2. Apply clipping and calculate the difference between the clipped signal and the original signal ("clippingDelta").
3. Convert the clippingDelta to the frequency domain (FFT), limit the magnitude of each frequency to the maskCurve, then convert it back to the time domain (IFFT).
4. Repeat 2 and 3 a few times to get the peaks down to the clipping level.
