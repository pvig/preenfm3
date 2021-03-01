
#ifndef FX_BUS_
#define FX_BUS_

#include "SynthStateAware.h"

class FxBus : public SynthStateAware {
public:
	FxBus();
	virtual ~FxBus() {}
    void init(SynthState *synthState);

	void mixSumInit();
    void mixSum(float *inStereo, int timbreNum);
	void processBlock(int32_t *outBuff);
    float interpolation(int delayReadPosInt);
    float phaser();
    void vcf(int delayReadPosInt);


	float* getSampleBlock() {
	    return sampleBlock_;
	}

	const float* getSampleBlock() const {
	    return sampleBlock_;
	}

	//lfo
	float lfo1;
	float lfo1Inc = 0.0001f;
	float lfo2;
	float lfo2Inc = 0.011f;
	float lfo3;
	float lfo3Inc = 0.000371f;
	float lfo4;
	float lfo4Inc = 0.037f;

protected:

	float sampleBlock_[BLOCK_SIZE * 2];
	float *sample;

    int fxType = 0;
    float fxTime = 0.98;
    float fxFeedback = 0.5;
    float fxTone = 0.5;
    float fxDiffusion = 0;
    float fxWidth =  0;

    float earlySpeed = 2;

	static const int earlyEchoNum = 4;
	static const int earlyEchoBufferSize = 1024 * 2;

	static float earlyEchoBuffer[earlyEchoBufferSize];
    int delayWritePos = 0;
    float delayReadPos = 0;
    int delayReadPosInt = 0;

    int earlyEchoSampleCount = fxTime * earlyEchoBufferSize;

    float srMod = PREENFM_FREQUENCY / 1000;

    float earlyEchoTimeList[4] = {
		 2.7f * srMod,
		 4.23f * srMod,
		 6.87f * srMod,
		 10.7f * srMod
	};
    float earlyEchoGainList[4] = {
		-0.45f,
		-0.55f,
		-0.33f,
		-0.13f
	};
    float earlyEchoFeebackList[4] = {
		-0.025f,
		0.056f,
		-0.03f,
		-0.01f
	};

	static const int recirculatingEchoNum = 4;
	static const int recirculatingBufferSize = 4360 * 2;

	static float recirculatingBuffer[recirculatingBufferSize];
    int delayCircWritePos = 0;
    float delayCircReadPos = 0;
    int delayCircReadPosInt = 0;
    int recirculatingMaxSampleCount = fxTime * recirculatingBufferSize;

    float recirculatingTimeList[4] = {
		 60.0f * srMod,
		 71.9345f * srMod,
		 86.7545f * srMod,
		 95.945f * srMod
    };
    float recirculatingGainList[4] = {
		0.75f,
		0.95f,
		0.35f,
		0.15f
	};
    int recirculatingCrossList[4][3] = {
		{ 1, 1, 2},
		{-1, 0, 3},
		{ 1, 3, 0},
		{-1, 3, 1}
	};

    // Filter
    float v0L, v1L, v2L, v3L, v4L, v5L, v6L, v7L, v8L;
    float v0R, v1R, v2R, v3R, v4R, v5R, v6R, v7R, v8R;
	float coef1L;
	float coef2L;
	float coef3L;
	float coef4L;
	float coef1R;
	float coef2R;
	float coef3R;
	float coef4R;

	float vcfFreq;
	float vcfDiffusion;

	float _ly1L ;
	float _ly1R ;
	float _ly2L ;
	float _ly2R ;
	float _ly3L ;
	float _ly3R ;
	float _ly4L ;
	float _ly4R ;
	float _lx1L ;
	float _lx1R ;
	float _lx2L ;
	float _lx2R ;
	float _lx3L ;
	float _lx3R ;
	float _lx4L ;
	float _lx4R ;

};

#endif	// end FX_BUS_
