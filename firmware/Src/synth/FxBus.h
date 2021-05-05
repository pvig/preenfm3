
#ifndef FX_BUS_
#define FX_BUS_

#include "SynthStateAware.h"

class FxBus : public SynthStateAware {
public:
	FxBus();
	virtual ~FxBus() {}
    void init(SynthState *synthState);

	void mixSumInit();
    void mixAdd(float *inStereo, int timbreNum);
	void processBlock(int32_t *outBuff);
    float feedbackHermiteInterpolation(int readPos);
    float fdbckMod();
    void vcf2L(int readPos);
    void vcf2R(int readPos);

	float* getSampleBlock() {
	    return sampleBlock_;
	}

	const float* getSampleBlock() const {
	    return sampleBlock_;
	}

protected:

	//lfo
	float lfo1, lfo1tri;
	float lfo1Inc = 0.00137521f;
	float lfo2;
	float lfo2Inc = 0.000113519845f;
	float lfoTremolo = 0, lfoTremoloSin = 0;
	float lfoTremoloInc = 0.00039845f;
	float tremoloEnvFollowAbs, tremoloEnvFollow = 0;

	float sampleBlock_[BLOCK_SIZE * 2];
	float *sample;

    float fxTime = 0.98, prevTime = -1;
    float prevFxFeedback = 0;
    float fxFeedback = 0.5;
    float fxTone = 0.25f;
    float fxDiffusion = 0.2f;
    float fxInputLevel = 0.5f, fxInputLevelAbs;
    float lfoDepth =  0;
    float fxSpeed = 0, prevSpeed = -1;
    float envMod, envModDepth, invtime = 1, invspeed = 1;
	float feedbackLp;
	float fxTremoloSpeed;
	float fxTremoloDepth;
	float fxCrossover;
	float envThreshold, envRelease, prevEnvThreshold = -1, prevEnvRelease = -1;
	float bounceLevel, prevBounce = -1, bouncingCv = 0;
	float timeCvControl = 0;
	float timeCv = 0, prevTimeCv = 0, timeCvSpeed = 0, prevtimeCvSpeed = 0, cvDelta;

	float fwL, fwR;
	float fbL, fbR;
	float combInR, combInL;
	float lpR, lpL;
	float lowcutR, lowcutL;
	float hpR, hpL;
	float inR, inL;
	float vcaR, vcaL;


	float envelope = 0;
	float blocksum = 0, envDest = 0, envM1, envM2;
	int envBlocknn = 0, envDetectSize = 32*32;

    float nodeL, nodeR, outL, outR;

	static const int feedbackSampleCount 	= 8192;
	static const int feedbackBufferSize 	= feedbackSampleCount * 2;
	static const int feedbackBufferSizeReadable = feedbackBufferSize - 6;

	static float feedbackBuffer[feedbackBufferSize];
    int feedbackWritePos 	= 0;
    float feedbackReadPosL 	= 0;
    float feedbackReadPosR 	= 0;
    int feedbackReadPosIntL = 0;
    int feedbackReadPosIntR = 0;
    int feedbackReadPosInt 	= 0;
    float feedbackDelayLen 	= 0;
    float feedbackFxTarget 	= 0;
	const float tremoloPanDepth 	= 0.07f;

    // Filter
    float v0L, v1L, v2L, v3L, v4L, v5L, v6L, v7L, v8L;
    float v0R, v1R, v2R, v3R, v4R, v5R, v6R, v7R, v8R;
    float f1L, f2L, f3L, f4L;
    float f1R, f2R, f3R, f4R;
	float coef1L;
	float coef2L;
	float coef3L;
	float coef4L;
	float coef1R;
	float coef2R;
	float coef3R;
	float coef4R;

	float inLpF, harmTremoloCutF;

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
