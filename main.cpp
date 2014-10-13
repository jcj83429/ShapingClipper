#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "aquila/global.h"
#include "aquila/source/WaveFile.h"
#include "aquila/source/FramesCollection.h"
#include "aquila/source/window/HannWindow.h"
#include "aquila/transform/FftFactory.h"
#include "aquila/tools/TextPlot.h"
#include "aquila/functions.h"
#include "ShapingClipper.h"

int main(int argc, char* argv[])
{
    if(argc < 3){
      printf("Usage: %s infile.wav outfile.pcm", argv[0]);
      return 1;
    }
    
    Aquila::SignalSource wavL = Aquila::WaveFile(argv[1], Aquila::StereoChannel::LEFT);
    Aquila::SignalSource wavR = Aquila::WaveFile(argv[1], Aquila::StereoChannel::RIGHT);

    ShapingClipper clipperL(44100, 256, 32768*20/100);
    ShapingClipper clipperR(44100, 256, 32768*20/100);
    const int feedSize = clipperL.getFeedSize();

    double outBufL[feedSize], outBufR[feedSize];
    int16_t out16[feedSize*2];
    Aquila::FramesCollection frmsL(wavL, feedSize);
    Aquila::FramesCollection frmsR(wavR, feedSize);

    std::ofstream outfile;
    outfile.open(argv[2], std::ios::out | std::ios::app | std::ios::binary);

    int frameN = 0;
    while(frameN < frmsL.count()){
      auto frmArr = frmsL.frame(frameN).toArray();
      clipperL.feed(frmArr, outBufL);

      frmArr = frmsR.frame(frameN).toArray();
      clipperR.feed(frmArr, outBufR);

      if(frameN > 2){
	for(int i = 0; i < feedSize; i++){
	  out16[2*i] = (outBufL[i] > 32766 ? 32766 : (outBufL[i] < -32766 ? -32766 : outBufL[i]));
	  out16[2*i+1] = (outBufR[i] > 32766 ? 32766 : (outBufR[i] < -32766 ? -32766 : outBufR[i]));
	}
	
	outfile.write((const char*)out16, feedSize*2*2);
      }

      frameN++;
    }
    
    return 0;

    double dummy[feedSize];
    for(int i=0; i<feedSize; i++)
      dummy[i] = 0.0;

    for(int i=0; i<3; i++){
      clipperL.feed(dummy, outBufL);
      for(int j=0; j<3; j++)
	out16[i] = outBufL[i];
      outfile.write((const char*)out16, feedSize*2);
    }
    
    return 0;
}
