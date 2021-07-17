
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
    float delayInterpolation(float readPos, float buffer[], int bufferLenM1);
    void lfoProcess(float *lfo, float *lfotri, float *lfoInc);
    float delayAllpassInterpolation(float readPos, float buffer[], int bufferLenM1, float prevVal);
	float* getSampleBlock() {
	    return sampleBlock_;
	}
	const float* getSampleBlock() const {
	    return sampleBlock_;
	}

protected:
	#define _dattorroSampleRateMod PREENFM_FREQUENCY / 29761.0f

	float headRoomMultiplier = 1;
	float headRoomDivider = 1;
	float sampleMultipler;

	float totalSent = 0;
	bool isActive = true;
	int silentBlockCount = 0;

	//lfo
	float lfo1 = 0;
	float lfo1tri;
	float lfo1Inc = 0.00037f;

	float lfo2 = 0;
	float lfo2tri;
	float lfo2Inc = lfo1Inc * 1.07f;

	float lfo3 = 1;
	float lfo3tri;
	float lfo3Inc = lfo1Inc;

	float lfo4 = 1;
	float lfo4tri;
	float lfo4Inc = lfo2Inc;


	float sampleBlock_[BLOCK_SIZE * 2];
	float *sample;

	int presetNum = -1, prevPresetNum = -1;
    float fxTime = 0.98, prevTime = -1, fxTimeLinear, prevFxTimeLinear = -1;
    float prevdecayVal = 0;
    float decayVal = 0.5, prevDecayVal = -1, decayFdbck = 0;
    const float decayMaxVal = 0.96f;
    float sizeParam, prevSizeParam, nextSizeParam;
    float inputDiffusion, prevInputDiffusion;
    float diffusion, prevDiffusion;
    float damping;
    float predelayMixLevel = 0.5f;
    float predelayMixAttn = predelayMixLevel * 0.75;
    float lfoDepth =  0;
    float lfoSpeed = 0, lfoSpeedLinear = 0, prevLfoSpeedLinear = -1;
	float loopLpf, loopLpf2, tiltInput, prevTilt = -1;
	float lfoLpCoef1 = 4999, lfoLpCoef2 = 0.0002f;
	float envThreshold, envRelease, prevEnvThreshold = -1, prevEnvRelease = -1;
	float envDecayMod = 0;
	float timeCvControl1 = 0, timeCvControl2 = 0, timeCvControl3 = 0, timeCvControl4 = 0;
	float timeCv = 0, prevTimeCv = 0, timeCvSpeed = 0, prevtimeCvSpeed = 0, cvDelta;

	float loopHpf = 0.04f;
	float inHpf = 0.1f;

	float inR, inL;

    float outL, outR;

    float int1 = 0, int2 = 0, int3 = 0, int4 = 0;

	static const int delay1BufferSize 	= 4453 * _dattorroSampleRateMod;
	static const int delay1BufferSizeM1	= delay1BufferSize - 1;
	static float delay1Buffer[delay1BufferSize];
    int delay1WritePos 	= 0;
    int delay1ReadPos 	= 0;
    float delay1DelayLen 	= 0, delay1ReadLen;
    float delay1FxTarget 	= 0;

	static const int delay2BufferSize 	= 3720 * _dattorroSampleRateMod;
	static const int delay2BufferSizeM1	= delay2BufferSize - 1;
	static float delay2Buffer[delay2BufferSize];
    int delay2WritePos 	= 0;
    int delay2ReadPos;
    float delay2DelayLen 	= 0, delay2ReadLen;
    float delay2FxTarget 	= 0;

	static const int delay3BufferSize 	= 4217 * _dattorroSampleRateMod;
	static const int delay3BufferSizeM1	= delay3BufferSize - 1;
	static float delay3Buffer[delay3BufferSize];
    int delay3WritePos 	= 0;
    int delay3ReadPos 	= 0;
    float delay3DelayLen 	= 0, delay3ReadLen;
    float delay3FxTarget 	= 0;

	static const int delay4BufferSize 	= 3163 * _dattorroSampleRateMod;
	static const int delay4BufferSizeM1	= delay4BufferSize - 1;
	static float delay4Buffer[delay4BufferSize];
    int delay4WritePos 	= 0;
    int delay4ReadPos;
    float delay4DelayLen 	= 0, delay4ReadLen;
    float delay4FxTarget 	= 0;

    const float dcBlockerCoef1 = 0.999f;
    const float dcBlockerCoef2 = 0.995f;
    const float dcBlockerCoef3 = 0.992f;

	//pre delay

	static const int predelayBufferSize 	= 16000;
	static const int predelayBufferSizeM1 	= predelayBufferSize - 1;
	static float predelayBuffer[predelayBufferSize];
    int predelaySize 			= predelayBufferSize;
    int predelayWritePos 		= 0;
    float predelayReadPos 		= 0;

	// input diffuser

	static const int inputBufferLen1 = 141 * _dattorroSampleRateMod;
	static const int inputBufferLen2 = 107 * _dattorroSampleRateMod;
	static const int inputBufferLen3 = 379 * _dattorroSampleRateMod;
	static const int inputBufferLen4 = 227 * _dattorroSampleRateMod;
	static float inputBuffer1[inputBufferLen1];
	static float inputBuffer2[inputBufferLen2];
	static float inputBuffer3[inputBufferLen3];
	static float inputBuffer4[inputBufferLen4];
	float inputBuffer1ReadLen = inputBufferLen1;
	float inputBuffer2ReadLen = inputBufferLen2;
	float inputBuffer3ReadLen = inputBufferLen3;
	float inputBuffer4ReadLen = inputBufferLen4;
    int inputWritePos1 	= 0;
    int inputWritePos2 	= 0;
    int inputWritePos3 	= 0;
    int inputWritePos4 	= 0;

    int inputReadPos1;
    int inputReadPos2;
    int inputReadPos3;
    int inputReadPos4;

    float inputCoef1 = 0.7f;
    float inputCoef2 = 0.625f;

	// diffuser decay

	static const int diffuserBufferLen1 = 672 	* _dattorroSampleRateMod;
	static const int diffuserBufferLen2 = 1800 	* _dattorroSampleRateMod;
	static const int diffuserBufferLen3 = 908 	* _dattorroSampleRateMod;
	static const int diffuserBufferLen4 = 2656 	* _dattorroSampleRateMod;
	static const int diffuserBufferLen1M1 = diffuserBufferLen1 - 1;
	static const int diffuserBufferLen2M1 = diffuserBufferLen2 - 1;
	static const int diffuserBufferLen3M1 = diffuserBufferLen3 - 1;
	static const int diffuserBufferLen4M1 = diffuserBufferLen4 - 1;
	float diffuserBuffer1ReadLen = diffuserBufferLen1M1;
	float diffuserBuffer2ReadLen = diffuserBufferLen2M1;
	float diffuserBuffer3ReadLen = diffuserBufferLen3M1;
	float diffuserBuffer4ReadLen = diffuserBufferLen4M1;

	float diffuserBuffer1ReadLen_b = diffuserBuffer1ReadLen * 0.1333f;//0.5f;
	float diffuserBuffer2ReadLen_b = diffuserBuffer2ReadLen * 0.3237f;//0.6180339887f;
	float diffuserBuffer3ReadLen_b = diffuserBuffer3ReadLen * 0.17f;//0.5f;
	float diffuserBuffer4ReadLen_b = diffuserBuffer4ReadLen * 0.5237f;//0.6180339887f;

	static float diffuserBuffer1[diffuserBufferLen1];
	static float diffuserBuffer2[diffuserBufferLen2];
	static float diffuserBuffer3[diffuserBufferLen3];
	static float diffuserBuffer4[diffuserBufferLen4];
    int diffuserWritePos1 	= 0;
    int diffuserWritePos2 	= 0;
    int diffuserWritePos3 	= 0;
    int diffuserWritePos4 	= 0;

    float diffuserReadPos1;
    float diffuserReadPos2;
    float diffuserReadPos3;
    float diffuserReadPos4;

    float diffuserCoef1 = 0.75f;
    float diffuserCoef2 = 0.65f;

	float monoIn, preDelayOut = 0, diff1Out, diff2Out, diff3Out, diff4Out;
    float ap1In, ap2In, ap3In, ap4In;
	float ap1Out, ap2Out, ap3Out, ap4Out;
	float feedbackInL = 0, feedbackInR = 0;

    // Filter
    float v0L, v1L, v2L, v3L, v4L, v5L, v6L, v7L, v8L;
    float v0R, v1R, v2R, v3R, v4R, v5R, v6R, v7R, v8R;
    float dcBlock1a, dcBlock1b, dcBlock2a, dcBlock2b, dcBlock3a, dcBlock3b, dcBlock4a, dcBlock4b, dcBlock5a, dcBlock5b;
    float hp_y0, hp_y1, hp_x1;
    float hp_a0, hp_a1, hp_b1;
	float inLpF;

	const int kl1 = 266 	* _dattorroSampleRateMod;
	const int kl2 = 2974 	* _dattorroSampleRateMod;
	const int kl3 = 1713 	* _dattorroSampleRateMod;
	const int kl4 = 1996 	* _dattorroSampleRateMod;
	const int kl5 = 1990 	* _dattorroSampleRateMod;
	const int kl6 = 187 	* _dattorroSampleRateMod;
	const int kl7 = 1066 	* _dattorroSampleRateMod;
    const long _kLeftTaps[7] = {kl1 , kl2, kl3, kl4, kl5, kl6, kl7};


	const int kr1 = 353 	* _dattorroSampleRateMod;
	const int kr2 = 3627 	* _dattorroSampleRateMod;
	const int kr3 = 1228 	* _dattorroSampleRateMod;
	const int kr4 = 2673 	* _dattorroSampleRateMod;
	const int kr5 = 2111 	* _dattorroSampleRateMod;
	const int kr6 = 335 	* _dattorroSampleRateMod;
	const int kr7 = 121 	* _dattorroSampleRateMod;
    const long _kRightTaps[7] = {kr1, kr2, kr3, kr4, kr5, kr6, kr7};

};

#endif	// end FX_BUS_
