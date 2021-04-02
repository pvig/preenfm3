
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
    float forwardBufferInterpolation(float readPos, bool isRight);
    float forwardHermiteInterpolation(int readPos);
    float feedbackBufferInterpolation(float readPos, bool isRight);
    float feedbackHermiteInterpolation(int readPos);
    float feedMod();
    float fdbckMod();
    void vcf1(int readPos);
    void vcf2L(int readPos);
    void vcf2R(int readPos);

	float* getSampleBlock() {
	    return sampleBlock_;
	}

	const float* getSampleBlock() const {
	    return sampleBlock_;
	}

	//lfo
	float lfo1;
	float lfo1Inc = 0.00137521f;
	float lfo2;
	float lfo2Inc = 0.000113519845f;

protected:

	float sampleBlock_[BLOCK_SIZE * 2];
	float *sample;

    float fxTime = 0.98;
    float prevFxForward = 0;
    float prevFxFeedback = 0;
    float fxFeedforward = 0.5;
    float fxFeedback = 0.5;
    float fxTone = 0.25f;
    float fxDiffusion = 0.2f;
    float fxInputLevel = 0.5f;
    float fxMod =  0;
    float fxSpeed = 0;
	float fxLp = 0.8;
	float fxTimeShift = 0;

	static const int forwardSampleCount = 2048;
	static const int forwardBufferSize = forwardSampleCount * 2;
	static float forwardBuffer[forwardBufferSize];

    int forwardWritePos = 0;
    float forwardReadPos = 0;
    float forwardDelayLen = 0;
    float forwardFxTarget = 0;
    int forwardReadPosInt = 0;

	static const int feedbackSampleCount = 2048;
	static const int feedbackBufferSize = feedbackSampleCount * 2;

	static float feedbackBuffer[feedbackBufferSize];
    int feedbackWritePosL = 0;
    float feedbackReadPosL = 0;
    int feedbackReadPosIntL = 0;
    float feedbackDelayLenL = 0;
    float feedbackFxTargetL = 0;

    int feedbackWritePosR = 0;
    float feedbackReadPosR = 0;
    int feedbackReadPosIntR = 0;
    float feedbackDelayLenR = 0;
    float feedbackFxTargetR = 0;

    float inputGainCoef = 0;

    // Filter
    float v0L, v1L, v2L, v3L, v4L, v5L, v6L, v7L, v8L;
    float v0R, v1R, v2R, v3R, v4R, v5R, v6R, v7R, v8R;
    float f1L,f2L,f3L,f4L;
    float f1R,f2R,f3R,f4R;
	float coef1L;
	float coef2L;
	float coef3L;
	float coef4L;
	float coef1R;
	float coef2R;
	float coef3R;
	float coef4R;

	float inLpF, inHpF;

	const float invBLOCKSIZE= (1/BLOCK_SIZE);
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
