#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstdarg>
#include "aquila/global.h"
#include "aquila/source/WaveFile.h"
#include "aquila/source/FramesCollection.h"
#include "aquila/source/window/HannWindow.h"
#include "aquila/transform/FftFactory.h"
#include "aquila/tools/TextPlot.h"
#include "aquila/functions.h"
#include "ShapingClipper.h"

struct WaveFmt{
	uint8_t audio_format[2];
	uint8_t num_channels[2];
	uint8_t sample_rate[4];
	uint8_t byte_rate[4];
	uint8_t block_align[2];
	uint8_t bits_per_sample[2];
};

void die(char *format, ...){
	va_list args;
	va_start(format, args);
	vprintf(format, args);
	va_end(args);
	printf("\n");
	exit(1);
}

/*	seek wave to the data chunk and returns the data chunk's size */
uint32_t seekWaveToData(FILE* waveFile, struct WaveFmt *fmt){
	char waveHeader[13] = {0};
	uint32_t fileLength;
	long chunkStart;
	uint32_t chunkSize;

	fread(waveHeader, 12, 1, waveFile);
	fileLength = 8 + ((uint32_t)(uint8_t)waveHeader[4]) | 
		(((uint32_t)(uint8_t)waveHeader[5]) << 8) | 
		(((uint32_t)(uint8_t)waveHeader[6]) << 16) | 
		(((uint32_t)(uint8_t)waveHeader[7]) << 24);
	waveHeader[4] = 0;
	if(strcmp(&waveHeader[0], "RIFF") != 0 || strcmp(&waveHeader[8], "WAVE") != 0){
		die("not a valid wave file");
	}

	while(ftell(waveFile) < fileLength){
		fread(waveHeader, 8, 1, waveFile);
		chunkStart = ftell(waveFile);
		chunkSize = ((uint32_t)(uint8_t)waveHeader[4]) | 
			(((uint32_t)(uint8_t)waveHeader[5]) << 8) | 
			(((uint32_t)(uint8_t)waveHeader[6]) << 16) | 
			(((uint32_t)(uint8_t)waveHeader[7]) << 24);
		waveHeader[4] = 0;
		if(!strcmp(&waveHeader[0], "fmt ")){
			fread(fmt, 16, 1, waveFile);
			if(fmt->audio_format[0] != 1 || fmt->audio_format[1] != 0)
				die("only PCM audio is supported (audio format is 0x%2x%2x)", fmt->audio_format[1], fmt->audio_format[0]);
			if(fmt->bits_per_sample[0] % 8 != 0)
				die("only 8, 16, 24 bit audio supported (audio file is %d bits)", fmt->bits_per_sample[0]);
		}else if(!strcmp(&waveHeader[0], "data")){
			printf("data chunk found at %ld\n", chunkStart);
			return chunkSize;
		}
		fseek(waveFile, chunkStart + chunkSize, SEEK_SET);
	}

	die("data chunk not found!");
	return -1;
}

void writeLengthField(FILE* outFile, uint32_t dataLength){
	uint32_t sizeField;
	uint8_t sizeFieldBytes[4];
	sizeField = dataLength;
	for(int i = 0; i < 4; i++){
		sizeFieldBytes[i] = (uint8_t)(0xff & sizeField);
		sizeField >>= 8;
	}

	fwrite(sizeFieldBytes, 4, 1, outFile);
}

void writeWaveHeader(FILE* outFile, struct WaveFmt *fmt, uint32_t dataLength){
	fwrite("RIFF", 4, 1, outFile);
	writeLengthField(outFile, 4 + 24 /*fmt chunk*/ + 8 + dataLength);
	fwrite("WAVE", 4, 1, outFile);
	fwrite("fmt ", 4, 1, outFile);
	writeLengthField(outFile, 16);
	fwrite(fmt, 16, 1, outFile);
	fwrite("data", 4, 1, outFile);
	writeLengthField(outFile, dataLength);
}

void splitInterleavedSamples(uint8_t *inBufInterleaved, std::vector<double*> inChannelBufs, int channels, int samples, int bytesPerSample){
	uint32_t signBit = 0x80 << ((bytesPerSample-1) * 8);
	for(int i = 0; i < samples; i++){
		for(int c = 0; c < channels; c++){
			uint32_t uSample = 0;
			for(int b = 0; b < bytesPerSample; b++){
				uSample |= ((uint32_t)inBufInterleaved[i * channels * bytesPerSample + c * bytesPerSample + b]) << (8 * b);
			}
			if(uSample & signBit){ // sign extend
				for(int b = bytesPerSample; b < 4; b++){
					uSample |= 0xff << (b * 8);
				}
			}
			inChannelBufs[c][i] = (double)(int32_t)uSample;
		}
	}
}

void writeSamples(FILE* outFile, std::vector<double*> outChannelBufs, int channels, int samples, int bytesPerSample){
	int ioBlockSize = channels * samples * bytesPerSample;
	uint8_t *outBufInterleaved = new uint8_t[ioBlockSize];
	int32_t maxAbsSampleVal = (0x80 << ((bytesPerSample - 1) * 8)) - 2;
	for(int i = 0; i < samples; i++){
		for(int c = 0; c < channels; c++){
			int32_t sample = (int32_t) outChannelBufs[c][i];
			sample = std::min(std::max(sample, -maxAbsSampleVal), maxAbsSampleVal);
			for(int b = 0; b < bytesPerSample; b++){
				outBufInterleaved[i * channels * bytesPerSample + c * bytesPerSample + b] = (sample >> (8 * b)) & 0xff;
			}
		}
	}
	fwrite(outBufInterleaved, 1, ioBlockSize, outFile);
}

int main(int argc, char* argv[])
{
    if(argc < 3){
      printf("Usage: %s infile.wav outfile.wav", argv[0]);
      return 1;
    }

	FILE* inFile = fopen(argv[1], "rb");
	FILE* outFile = fopen(argv[2], "wb");
	int clipLevel = 50;

	if(inFile == NULL)
		die("cannot open input file");

	if(outFile == NULL)
		die("cannot open output file for writing");

	struct WaveFmt fmt;
	uint32_t dataLength = seekWaveToData(inFile, &fmt);
	int channels = fmt.num_channels[0];
	int sampleRate = ((uint32_t)fmt.sample_rate[0]) | 
		(((uint32_t)fmt.sample_rate[1]) << 8) | 
		(((uint32_t)fmt.sample_rate[2]) << 16) | 
		(((uint32_t)fmt.sample_rate[3]) << 24);
	int bytesPerSample = fmt.bits_per_sample[0] / 8;
	int32_t fullScale = 1 << (fmt.bits_per_sample[0] - 1);
	printf("%d ch, %d Hz, %d Bps\n", channels, sampleRate, bytesPerSample);

	std::vector<ShapingClipper*> clippers;
	for(int i = 0; i < channels; i++)
		clippers.push_back(new ShapingClipper(sampleRate, 256, 32768*50/100));
	const int feedSize = clippers[0]->getFeedSize();

	std::vector<double*> inChannelBufs;
	std::vector<double*> outChannelBufs;
	for(int i = 0; i < channels; i++){
		inChannelBufs.push_back(new double[feedSize]);
		outChannelBufs.push_back(new double[feedSize]);
	}

	int ioBlockSize = channels * feedSize * bytesPerSample;
	uint8_t *inBufInterleaved = new uint8_t[ioBlockSize];
	uint8_t *outBufInterleaved = new uint8_t[ioBlockSize];

	writeWaveHeader(outFile, &fmt, dataLength);

	int bytesRead;
	int count = 0;

	while((bytesRead = fread(inBufInterleaved, 1, ioBlockSize, inFile)) == ioBlockSize){
		splitInterleavedSamples(inBufInterleaved, inChannelBufs, channels, feedSize, bytesPerSample);

		for(int c = 0; c < channels; c++){
			clippers[c]->feed(inChannelBufs[c], outChannelBufs[c]);
		}

		if(count < 3){ // due to FFT delay, first 3 output blocks are empty
			count++;
			continue;
		}

		writeSamples(outFile, outChannelBufs, channels, feedSize, bytesPerSample);
	}

	// if last block is incomplete, it will be processed specially
	if(bytesRead > 0){
		for(int i = bytesRead; i < ioBlockSize; i++)
			inBufInterleaved[i] = 0;
		splitInterleavedSamples(inBufInterleaved, inChannelBufs, channels, feedSize, bytesPerSample);
		for(int c = 0; c < channels; c++){
			clippers[c]->feed(inChannelBufs[c], outChannelBufs[c]);
		}
		writeSamples(outFile, outChannelBufs, channels, feedSize, bytesPerSample);
	}

	int samplesRemaining = 2*feedSize + bytesRead/(channels * bytesPerSample);

	for(int c = 0; c < channels; c++){
		for(int i = 0; i < feedSize; i++){
			inChannelBufs[c][i] = 0;
		}
	}
	for(int i = 0; i < 3; i++){
		for(int c = 0; c < channels; c++){
			clippers[c]->feed(inChannelBufs[c], outChannelBufs[c]);
		}
		writeSamples(outFile, outChannelBufs, channels, (samplesRemaining > feedSize ? feedSize : samplesRemaining), bytesPerSample);
		samplesRemaining -= feedSize;
	}


	fclose(inFile);
	fclose(outFile);

	// don't bother freeing the buffers

	return 0;
}
