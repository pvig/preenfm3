
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
    float delay1HermiteInterpolation(int readPos);
    float delay1Interpolation(float readPos);
    float delay2Interpolation(float readPos);

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
	float lfo2tri, lfo2btri;
	float lfo2, lfo2b;
	float lfo2Inc = 0.000113519845f;
	float lfoTremolo = 0, lfoTremoloSin = 0;
	float lfoTremoloInc = 0.00039845f;
	float tremoloEnvFollowAbs, tremoloEnvFollow = 0;

	float sampleBlock_[BLOCK_SIZE * 2];
	float *sample;

    float fxTime = 0.98, prevTime = -1, fxTimeLinear;
    float prevfeedbackGain = 0;
    float feedbackGain = 0.5;
    float fxTone = 0.25f;
    float fxDiffusion = 0.2f;
    float fxInputLevel = 0.5f, fxInputLevelAbs;
    float lfoDepth =  0;
    float fxSpeed = 0, prevSpeed = -1;
    float envMod, envModDepth, invtime = 1, invspeed = 1;
	float delay1Lp;
	float fxTremoloSpeed;
	float fxTremoloDepth;
	float fxCrossover;
	float pingpongFactor;
	float envThreshold, envRelease, prevEnvThreshold = -1, prevEnvRelease = -1;
	float bounceLevel, prevBounce = -1, bouncingCv = 0;
	float timeCvControl = 0;
	float timeCv = 0, prevTimeCv = 0, timeCvSpeed = 0, prevtimeCvSpeed = 0, cvDelta;
	float spread, ratio;

	float fbPoint, fbPoint2;
	float combInR, combInL;
	float lpR, lpL;
	float lowcutR, lowcutL;
	float hpR, hpL;
	float inR, inL;
	float loopHp, hpVal;
	float vcaR, vcaL;


	float envelope = 0;
	float blocksum = 0, envDest = 0, envM1, envM2;
	int envBlocknn = 0, envDetectSize = 32*32;

    float nodeL, nodeR, outL, outR;

	static const int delay1BufferSize 	= 4000;
	static const int delay1BufferSizeM1	= delay1BufferSize - 1;
	static float delay1Buffer[delay1BufferSize];
    int delay1WritePos 	= 0;
    float delay1ReadPos 	= 0;
    int delay1ReadPosInt = 0;
    float delay1DelayLen 	= 0;
    float delay1FxTarget 	= 0;

	static const int delay2BufferSize 	= 3500;
	static const int delay2BufferSizeM1	= delay2BufferSize - 1;
	static float delay2Buffer[delay2BufferSize];
    int delay2WritePos 	= 0;
    int delay2ReadPosInt;
    float delay2ReadPos 	= 0;
    float delay2DelayLen 	= 0;
    float delay2FxTarget 	= 0;

	const float tremoloPanDepth 	= 0.07f;

	// tap delay input

	static const int tapDelayBufferSize 	= 4096;
	static float tapDelayBuffer[tapDelayBufferSize];
    int tapDelayWritePos 		= 0;

	static const int tapCount 	= 6;
	static int tapDelayPos[tapCount];
	static float tapDelayAmp[tapCount];

	// diffuser

	static const int diffuserBufferLen1 = 229;
	static const int diffuserBufferLen2 = 173;
	static const int diffuserBufferLen3 = 611;
	static const int diffuserBufferLen4 = 447;

	static float diffuserBuffer1[diffuserBufferLen1];
	static float diffuserBuffer2[diffuserBufferLen2];
	static float diffuserBuffer3[diffuserBufferLen3];
	static float diffuserBuffer4[diffuserBufferLen4];
    int diffuserWritePos1 	= 0;
    int diffuserWritePos2 	= 0;
    int diffuserWritePos3 	= 0;
    int diffuserWritePos4 	= 0;

    int diffuserReadPos1;
    int diffuserReadPos2;
    int diffuserReadPos3;
    int diffuserReadPos4;

    float diffuserCoef1 = 0.5f, diffuserCoef1b;
    float diffuserCoef2 = 0.55f, diffuserCoef2b;

    int rand1 = 0,rand2 = 0,rand3 = 0;

    // Filter
    float v0L, v1L, v2L, v3L, v4L, v5L, v6L, v7L, v8L;
    float v0R, v1R, v2R, v3R, v4R, v5R, v6R, v7R, v8R;

    float lowL = 0, highL = 0, bandL = 0;
    float lowR = 0, highR = 0, bandR = 0;

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
