
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
    float forwardBufferInterpolation(int readPosInt);
    float feedbackBufferInterpolation(int readPosInt);
    float phaser();
    void vcf(int forwardReadPosInt);


	float* getSampleBlock() {
	    return sampleBlock_;
	}

	const float* getSampleBlock() const {
	    return sampleBlock_;
	}

	//lfo
	float lfo1;
	float lfo1Inc = 0.00015f;
	/*float lfo2;
	float lfo2Inc = 0.0011f;
	float lfo3;
	float lfo3Inc = 0.0000371f;
	float lfo4;
	float lfo4Inc = 0.0177f;*/

protected:

	float sampleBlock_[BLOCK_SIZE * 2];
	float *sample;

    float fxTime = 0.98;
    float prevFxTime;
    float fxFeedforward = 0.5;
    float fxFeedback = 0.5;
    float fxTone = 0.25f;
    float fxDiffusion = 0.2f;
    float fxInputLevel = 0.5f;
    float fxMod =  0;
    float fxSpeed = 0;

    const float bufferInc = 2;

	static const int forwardBufferSize = 4000 * 2;

	static float forwardBuffer[forwardBufferSize];
    int forwardWritePos = 0;
    float forwardReadPos = 0;
    float forwardDelayLen = 0;
    int forwardReadPosInt = 0;

    int forwardSampleCount = fxTime * forwardBufferSize;

    //float srMod = PREENFM_FREQUENCY / 5000;

	static const int feedbackBufferSize = 4000 * 2;

	static float feedbackBuffer[feedbackBufferSize];
    int feedbackWritePos = 0;
    float feedbackReadPos = 0;
    int feedbackReadPosInt = 0;
    float feedbackDelayLen = 0;
    float inputGainCoef = 0.07f;

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