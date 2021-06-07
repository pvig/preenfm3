
#ifndef FX_BUS_
#define FX_BUS_

#include "SynthStateAware.h"

class FxBus : public SynthStateAware {
public:
	FxBus();
	virtual ~FxBus() {}
    void init(SynthState *synthState);

	void mixSumInit();
    void mixAdd(float *inStereo, int timbreNum, int32_t *outBuffer);
	void processBlock(int32_t *outBuff);
    float delay1HermiteInterpolation(int readPos);
    float delay1Interpolation(float readPos);
    float delay2Interpolation(float readPos);
    float delay3Interpolation(float readPos);
    float delay4Interpolation(float readPos);
    float predelayInterpolation(float readPos);
    float diffuser1Interpolation(float readPos);
    float diffuser1CubicInterpolation(float readPos);
    float diffuser2Interpolation(float readPos);
    float diffuser3Interpolation(float readPos);
    float diffuser3CubicInterpolation(float readPos);
    float diffuser4Interpolation(float readPos);
    float * crossFade(float t);

	float* getSampleBlock() {
	    return sampleBlock_;
	}

	const float* getSampleBlock() const {
	    return sampleBlock_;
	}

protected:
	#define _dattorroSampleRateMod 1.6100265448f //PREENFM_FREQUENCY / 29761.0f

	//lfo
	float lfo1, lfo1tri;
	float lfo1btri, lfo1b;
	float lfo1Inc = 0.000137521f;
	float lfo2tri, lfo2btri;
	float lfo2, lfo2b;
	const float lfo2Inc = 0.00001041666667f;
	float lfo2IncModSampleInc = 0;
	float lfo2IncMod;
	float lfo2ModVal;
	int   lfo2ChangeCounter = 0;
	const int   lfo2ChangePeriod = 17733;
	const int   lfo2ChangePeriodInv = 1 / lfo2ChangePeriod;
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
    float inputDiffusion, prevInputDiffusion;
    float decayDiffusion, prevDecayDiffusion;
    float damping;
    float predelayMixLevel = 0.5f;
    float predelayMixAttn = predelayMixLevel * 0.75;
    float lfoDepth =  0;
    float fxSpeed = 0, prevSpeed = -1;
    float envMod, envModDepth, invtime = 1, invspeed = 1;
	float loopLpf, loopLpf2, loopHpf,  monoInHpf;
	float fxTremoloSpeed;
	float fxTremoloDepth;
	float fxCrossover;
	float pingpongFactor;
	float envThreshold, envRelease, prevEnvThreshold = -1, prevEnvRelease = -1;
	float bounceLevel, prevBounce = -1, bouncingCv = 0;
	float timeCvControl = 0, timeCvControlb = 0;
	float timeCv = 0, prevTimeCv = 0, timeCvSpeed = 0, prevtimeCvSpeed = 0, cvDelta;
	float spread, ratio;
	const float decoupler1 = 0.0047f;
	const float decoupler2 = 0.0033f;
	float loopDecoupler = 0;
	float loopDecoupler2 = 0;
	float loopDecouplerModVal;
	int   loopDecouplerChangeCounter = 0;
	const int   loopDecouplerChangePeriod = 15733;

	float fbPoint, fbPoint2;
	float combInR, combInL;
	float tremoloOut;
	float lpR, lpL;
	float lowcutR, lowcutL;
	float hpR, hpL;
	float inR, inL;
	float loopHp, hpVal;
	float vca, vcaR, vcaL;


	float envelope = 0;
	float blocksum = 0, envDest = 0, envM1, envM2;
	int envBlocknn = 0, envDetectSize = 32*32;

    float nodeL, nodeR, outL, outR;

	static const int delay1BufferSize 	= 4453 * _dattorroSampleRateMod;
	static const int delay1BufferSizeM1	= delay1BufferSize - 1;
	static float delay1Buffer[delay1BufferSize];
    int delay1WritePos 	= 0;
    float delay1ReadPos 	= 0;
    int delay1ReadPosInt = 0;
    float delay1DelayLen 	= 0;
    float delay1FxTarget 	= 0;

	static const int delay2BufferSize 	= 3720 * _dattorroSampleRateMod;
	static const int delay2BufferSizeM1	= delay2BufferSize - 1;
	static float delay2Buffer[delay2BufferSize];
    int delay2WritePos 	= 0;
    int delay2ReadPos;
    float delay2DelayLen 	= 0;
    float delay2FxTarget 	= 0;

	static const int delay3BufferSize 	= 4217 * _dattorroSampleRateMod;
	static const int delay3BufferSizeM1	= delay3BufferSize - 1;
	static float delay3Buffer[delay3BufferSize];
    int delay3WritePos 	= 0;
    int delay3ReadPosInt;
    float delay3ReadPos 	= 0;
    float delay3DelayLen 	= 0;
    float delay3FxTarget 	= 0;

	static const int delay4BufferSize 	= 3163 * _dattorroSampleRateMod;
	static const int delay4BufferSizeM1	= delay4BufferSize - 1;
	static float delay4Buffer[delay4BufferSize];
    int delay4WritePos 	= 0;
    int delay4ReadPos;
    float delay4DelayLen 	= 0;
    float delay4FxTarget 	= 0;


	//pre delay

	static const int predelayBufferSize 	= 16000;
	static const int predelayBufferSizeM1 	= predelayBufferSize - 1;
	static float predelayBuffer[predelayBufferSize];
    int predelaySize 			= predelayBufferSize;
    int predelaySizeM1 			= predelayBufferSize - 1;
    int predelayWritePos 		= 0;
    float predelayReadPos 		= 0;

	// input diffuser

	static const int inputBufferLen1 = 141 * _dattorroSampleRateMod;
	static const int inputBufferLen2 = 107 * _dattorroSampleRateMod;
	static const int inputBufferLen3 = 379 * _dattorroSampleRateMod;
	static const int inputBufferLen4 = 227 * _dattorroSampleRateMod;
	static const int inputBufferLen1M1 = inputBufferLen1 - 1;
	static const int inputBufferLen2M1 = inputBufferLen2 - 1;
	static const int inputBufferLen3M1 = inputBufferLen3 - 1;
	static const int inputBufferLen4M1 = inputBufferLen4 - 1;

	static float inputBuffer1[inputBufferLen1];
	static float inputBuffer2[inputBufferLen2];
	static float inputBuffer3[inputBufferLen3];
	static float inputBuffer4[inputBufferLen4];
    int inputWritePos1 	= 0;
    int inputWritePos2 	= 0;
    int inputWritePos3 	= 0;
    int inputWritePos4 	= 0;

    int inputReadPos1;
    int inputReadPos2;
    int inputReadPos3;
    int inputReadPos4;

    float inputCoef1 = 0.7f, inputCoef1b;
    float inputCoef2 = 0.625f, inputCoef2b;

	// diffuser decay

	static const int diffuserBufferLen1 = 672 	* _dattorroSampleRateMod;
	static const int diffuserBufferLen2 = 1800 	* _dattorroSampleRateMod;
	static const int diffuserBufferLen3 = 908 	* _dattorroSampleRateMod;
	static const int diffuserBufferLen4 = 2656 	* _dattorroSampleRateMod;
	static const int diffuserBufferLen1M1 = diffuserBufferLen1 - 1;
	static const int diffuserBufferLen1M4 = diffuserBufferLen1 - 4;
	static const int diffuserBufferLen2M1 = diffuserBufferLen2 - 1;
	static const int diffuserBufferLen3M1 = diffuserBufferLen3 - 1;
	static const int diffuserBufferLen3M4 = diffuserBufferLen3 - 4;
	static const int diffuserBufferLen4M1 = diffuserBufferLen4 - 1;

	static float diffuserBuffer1[diffuserBufferLen1];
	static float diffuserBuffer2[diffuserBufferLen2];
	static float diffuserBuffer3[diffuserBufferLen3];
	static float diffuserBuffer4[diffuserBufferLen4];
    int diffuserWritePos1 	= 0;
    int diffuserWritePos2 	= 0;
    int diffuserWritePos3 	= 0;
    int diffuserWritePos4 	= 0;

    float diffuserReadPos1;
    int diffuserReadPos2;
    float diffuserReadPos3;
    int diffuserReadPos4;

    float diffuserCoef1 = 0.75f, diffuserCoef1b;
    float diffuserCoef2 = 0.65f, diffuserCoef2b;

	float monoIn, diff1Out, diff2Out, diff3Out, diff4Out;
    float ap1In, ap2In, ap3In, ap4In;
	float ap1Out, ap2Out, ap3Out, ap4Out;
	float feedbackInL = 0, feedbackInR = 0;

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

	const int kl1 = 266;// * _dattorroSampleRateMod;
	const int kl2 = 2974;// * _dattorroSampleRateMod;
	const int kl3 = 1913;// * _dattorroSampleRateMod;
	const int kl4 = 1996;// * _dattorroSampleRateMod;
	const int kl5 = 1990;// * _dattorroSampleRateMod;
	const int kl6 = 187;// * _dattorroSampleRateMod;
	const int kl7 = 1066;// * _dattorroSampleRateMod;
    const long _kLeftTaps[7] = {kl1 , kl2, kl3, kl4, kl5, kl6, kl7};


	const int kr1 = 353;// * _dattorroSampleRateMod;
	const int kr2 = 3627;// * _dattorroSampleRateMod;
	const int kr3 = 1228;// * _dattorroSampleRateMod;
	const int kr4 = 2673;// * _dattorroSampleRateMod;
	const int kr5 = 2111;// * _dattorroSampleRateMod;
	const int kr6 = 335;// * _dattorroSampleRateMod;
	const int kr7 = 121;// * _dattorroSampleRateMod;

    //const long _kRightTaps[7] = {266, 2974, 1913, 1996, 1990, 187, 1066};
    const long _kRightTaps[7] = {kr1, kr2, kr3, kr4, kr5, kr6, kr7};
};

#endif	// end FX_BUS_
