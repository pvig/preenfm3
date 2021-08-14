/*
 * FxBus.cpp
 *
 *  Created on: Feb 14, 2021
 *      Author: patvig
 */

#include "FxBus.h"

extern float diatonicScaleFrequency[];

inline float fold(float x4) {
    // https://www.desmos.com/calculator/ge2wvg2wgj
    // x4 = x / 4
    return (fabsf(x4 + 0.25f - roundf(x4 + 0.25f)) - 0.25f);
}
inline
int modulo(int d, int max) {
  return unlikely(d >= max) ? d - max : d;
}
inline
float modulo2(float readPos, int bufferLen) {
    if( unlikely(readPos < 0) )
    	readPos += bufferLen;
    return readPos;
}
inline
float clamp(float d, float min, float max) {
  const float t = unlikely(d < min) ? min : d;
  return unlikely(t > max) ? max : t;
}
inline
float sqrt3(const float x)
{
  union
  {
    int i;
    float x;
  } u;

  u.x = x;
  u.i = (1 << 29) + (u.i >> 1) - (1 << 22);
  return u.x;
}
inline
float fastroot(float f,int n)
{
 long *lp,l;
 lp=(long*)(&f);
 l=*lp;l-=0x3F800000l;l>>=(n-1);l+=0x3F800000l;
 *lp=l;
 return f;
}
inline
float sigmoid(float x)
{
    return x * (1.5f - 0.5f * x * x);
}
//***------------***------------***------------***------------***----- FxBus -------***------------***------------***------------

float FxBus::delay1Buffer[delay1BufferSize] __attribute__((section(".ram_d1")));
float FxBus::delay2Buffer[delay2BufferSize] __attribute__((section(".ram_d1")));
float FxBus::delay3Buffer[delay3BufferSize] __attribute__((section(".ram_d2")));
float FxBus::delay4Buffer[delay4BufferSize] __attribute__((section(".ram_d2")));

float FxBus::predelayBuffer[predelayBufferSize] __attribute__((section(".ram_d1")));

float FxBus::inputBuffer1[inputBufferLen1] __attribute__((section(".ram_d2")));
float FxBus::inputBuffer2[inputBufferLen2] __attribute__((section(".ram_d2")));
float FxBus::inputBuffer3[inputBufferLen3] __attribute__((section(".ram_d2")));
float FxBus::inputBuffer4[inputBufferLen4] __attribute__((section(".ram_d2")));

float FxBus::diffuserBuffer1[diffuserBufferLen1] __attribute__((section(".ram_d2")));
float FxBus::diffuserBuffer2[diffuserBufferLen2] __attribute__((section(".ram_d2")));
float FxBus::diffuserBuffer3[diffuserBufferLen3] __attribute__((section(".ram_d2")));
float FxBus::diffuserBuffer4[diffuserBufferLen4] __attribute__((section(".ram_d2")));

FxBus::FxBus() {}

void FxBus::init(SynthState *synthState) {
    this->synthState_ = synthState;

	for (int s = 0; s < predelayBufferSize; s++) {
    	predelayBuffer[s] = 0;
    }
    for (int s = 0; s < delay1BufferSize; s++) {
    	delay1Buffer[ s ] = 0;
    }
    for (int s = 0; s < delay2BufferSize; s++) {
    	delay2Buffer[ s ] = 0;
    }
    for (int s = 0; s < delay3BufferSize; s++) {
    	delay3Buffer[ s ] = 0;
    }
    for (int s = 0; s < delay4BufferSize; s++) {
    	delay4Buffer[ s ] = 0;
    }

    for (int s = 0; s < inputBufferLen1; s++) {
    	inputBuffer1[ s ] = 0;
    }
    for (int s = 0; s < inputBufferLen2; s++) {
    	inputBuffer2[ s ] = 0;
    }
    for (int s = 0; s < inputBufferLen3; s++) {
    	inputBuffer3[ s ] = 0;
    }
    for (int s = 0; s < inputBufferLen4; s++) {
    	inputBuffer4[ s ] = 0;
    }

    for (int s = 0; s < diffuserBufferLen1; s++) {
    	diffuserBuffer1[ s ] = 0;
    }
    for (int s = 0; s < diffuserBufferLen2; s++) {
    	diffuserBuffer2[ s ] = 0;
    }
    for (int s = 0; s < diffuserBufferLen3; s++) {
    	diffuserBuffer3[ s ] = 0;
    }
    for (int s = 0; s < diffuserBufferLen4; s++) {
    	diffuserBuffer4[ s ] = 0;
    }

    // Init FX variables

    v0L = 0;
    v1L = 0;
    v0R = 0;
    v1R = 0;
    v2L = 0;
    v3L = 0;
    v2R = 0;
    v3R = 0;

    v4L = 0;
    v5L = 0;
    v4R = 0;
	v5R = 0;
    v6L = 0;
    v6R = 0;
    v7L = 0;
	v7R = 0;

	dcBlock1a = 0;
	dcBlock1b = 0;
	dcBlock2a = 0;
	dcBlock2b = 0;
	dcBlock3a = 0;
	dcBlock3b = 0;
	dcBlock4a = 0;
	dcBlock4b = 0;
	dcBlock5a = 0;
	dcBlock5b = 0;

	hp_y0 = 0;
	hp_y1 = 0;
	hp_x1 = 0;

}
/**
 * init before timbres summing
 */
void FxBus::mixSumInit() {

	if(!isActive) {
		return;
	}

    float temp, sizeParamInpt, sizeSqrt;
	sample = getSampleBlock();

    for (int s = 0; s < 8; s++) {
    	*(sample++) = 0;
    	*(sample++) = 0;
    	*(sample++) = 0;
    	*(sample++) = 0;
    	*(sample++) = 0;
    	*(sample++) = 0;
    	*(sample++) = 0;
    	*(sample++) = 0;
    }

    // assign user param to variables :

    // ------ page 1

	presetNum = synthState_->fullState.masterfxConfig[GLOBALFX_PRESETNUM];

	if (prevPresetNum != presetNum)
	{
		if (presetNum < 15)
		{
			int size = presetNum * 0.333333333f;
			int brightness = presetNum % 3;

			switch (size)
			{
			case 0:
				//synthState_->fullState.masterfxConfig[GLOBALFX_INPUTWIDTH] = 0.4f;

				synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHBASE] = 0.57f;
				//synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHSPREAD] = 0.23f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFOSPEED] = 0.1f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFODEPTH] = 0.9;

				switch (presetNum)
				{
				case 0:
					synthState_->fullState.masterfxConfig[GLOBALFX_SIZE] = 0.1f;
					synthState_->fullState.masterfxConfig[GLOBALFX_DECAY] = 0.04f;
					synthState_->fullState.masterfxConfig[GLOBALFX_DIFFUSION] = 0.6f;
					break;
				case 1:
					synthState_->fullState.masterfxConfig[GLOBALFX_SIZE] = 0.13f;
					synthState_->fullState.masterfxConfig[GLOBALFX_DECAY] = 0.02f;
					synthState_->fullState.masterfxConfig[GLOBALFX_DIFFUSION] = 0.9f;
					break;
				case 2:
					synthState_->fullState.masterfxConfig[GLOBALFX_SIZE] = 0.23f;
					synthState_->fullState.masterfxConfig[GLOBALFX_DECAY] = 0.08f;
					synthState_->fullState.masterfxConfig[GLOBALFX_DIFFUSION] = 0.6f;
					break;
				default:
					break;
				}
				break;
			case 1:
				//synthState_->fullState.masterfxConfig[GLOBALFX_INPUTWIDTH] = 0.4f;
				synthState_->fullState.masterfxConfig[GLOBALFX_SIZE] = 0.30f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DECAY] = 0.2f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DIFFUSION] = 0.66f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFOSPEED] = 0.55f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFODEPTH] = 0.12f;
				synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHBASE] = 0.52f;
				//synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHSPREAD] = 0.23f;
				break;
			case 2:
				//synthState_->fullState.masterfxConfig[GLOBALFX_INPUTWIDTH] = 0.4f;
				synthState_->fullState.masterfxConfig[GLOBALFX_SIZE] = 0.465f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DECAY] = 0.3f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DIFFUSION] = 0.7f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFOSPEED] = 0.57f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFODEPTH] = 0.21;
				synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHBASE] = 0.4f;
				//synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHSPREAD] = 0.23f;
				break;
			case 3:
				//synthState_->fullState.masterfxConfig[GLOBALFX_INPUTWIDTH] = 0.42f;
				synthState_->fullState.masterfxConfig[GLOBALFX_SIZE] = 0.775f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DECAY] = 0.46f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DIFFUSION] = 0.76f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFOSPEED] = 0.57f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFODEPTH] = 0.21;
				synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHBASE] = 0.48f;
				//synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHSPREAD] = 0.23f;
				break;
			case 4:
				//synthState_->fullState.masterfxConfig[GLOBALFX_INPUTWIDTH] = 0.43f;
				synthState_->fullState.masterfxConfig[GLOBALFX_SIZE] = 0.97f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DECAY] = 0.72f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DIFFUSION] = 0.93f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFOSPEED] = 0.57f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFODEPTH] = 0.21;
				synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHBASE] = 0.196f;
				//synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHSPREAD] = 0.23f;
				break;
			default:
				break;
			}

			switch (brightness)
			{
			case 0:
				synthState_->fullState.masterfxConfig[GLOBALFX_DAMPING] = 0.4f;
				//synthState_->fullState.masterfxConfig[GLOBALFX_INPUTWIDTH] *= 0.9f;
				break;
			case 1:
				synthState_->fullState.masterfxConfig[GLOBALFX_DAMPING] = 0.62f * (1 - synthState_->fullState.masterfxConfig[GLOBALFX_SIZE] * 0.06f);
				//synthState_->fullState.masterfxConfig[GLOBALFX_INPUTWIDTH] *= 0.95f;
				break;
			case 2:
				synthState_->fullState.masterfxConfig[GLOBALFX_DAMPING] = 0.9f * (1 - synthState_->fullState.masterfxConfig[GLOBALFX_SIZE] * 0.1f);
				break;
			default:
				break;
			}

		}
		else
		{
			switch (presetNum)
			{
			case 15:
				//freeze
				synthState_->fullState.masterfxConfig[GLOBALFX_DAMPING] = 0.96f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFOSPEED] = 0.87f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFODEPTH] = 0.05f;
				synthState_->fullState.masterfxConfig[GLOBALFX_SIZE] = 1;
				synthState_->fullState.masterfxConfig[GLOBALFX_DECAY] = 0.943f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DIFFUSION] = 0.05f;
				synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHBASE] = 0.38f;
				//synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHSPREAD] = 0.68f;
				break;
			case 16:
				//hall
				synthState_->fullState.masterfxConfig[GLOBALFX_DAMPING] = 0.57f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFOSPEED] = 0.28f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFODEPTH] = 0.27f;
				synthState_->fullState.masterfxConfig[GLOBALFX_SIZE] = 0.33f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DECAY] = 0.49f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DIFFUSION] = 0.66f;
				synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHBASE] = 0.1055f;
				//synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHSPREAD] = 0.22f;
				break;
			case 17:
				//cave
				synthState_->fullState.masterfxConfig[GLOBALFX_DAMPING] = 0.38f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFOSPEED] = 0.37f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFODEPTH] = 0.33f;
				synthState_->fullState.masterfxConfig[GLOBALFX_SIZE] = 0.94f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DECAY] = 0.43f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DIFFUSION] = 0.8f;
				synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHBASE] = 0.22f;
				//synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHSPREAD] = 0.46f;
				break;
			case 18:
				//apartment
				synthState_->fullState.masterfxConfig[GLOBALFX_DAMPING] = 0.44f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFOSPEED] = 0.33f;
				synthState_->fullState.masterfxConfig[GLOBALFX_LFODEPTH] = 0.33f;
				synthState_->fullState.masterfxConfig[GLOBALFX_SIZE] = 0.6f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DECAY] = 0.22f;
				synthState_->fullState.masterfxConfig[GLOBALFX_DIFFUSION] = 0.72f;
				synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHBASE] = 0.31f;
				//synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHSPREAD] = 0.48f;
				break;
			default:
				break;
			}
		}
	}
	prevPresetNum = presetNum;

    inputWidth 		= 	synthState_->fullState.masterfxConfig[GLOBALFX_INPUTWIDTH] * 0.6f;

	float tilt = synthState_->fullState.masterfxConfig[GLOBALFX_INPUTBASE];
	if (prevTilt != tilt || prevInputWidth != inputWidth)
	{
		tiltInput = (tilt * tilt * 0.6f);
		inHpf = clamp(tiltInput, 0, 1);
		inLpF = clamp(tiltInput + (inputWidth * inputWidth) , 0, 1);
	}
	prevTilt = tilt;
	prevInputWidth = inputWidth;

	fxTimeLinear = synthState_->fullState.masterfxConfig[GLOBALFX_PREDELAYTIME];
	if (prevFxTimeLinear == fxTimeLinear)
	{
		prevTime = clamp(fxTimeLinear, 0.0003f, 0.9996f);
		prevTime *= prevTime * prevTime;
		fxTime = fxTime * 0.9f + prevTime * 0.1f;
		predelaySize = fxTimeLinear * predelayBufferSizeM1;
	}
	prevFxTimeLinear = fxTimeLinear;

	predelayMixLevel = synthState_->fullState.masterfxConfig[GLOBALFX_PREDELAYMIX];
	predelayMixAttn = predelayMixLevel * (1 - (predelayMixLevel * predelayMixLevel * 0.1f));

	lfoSpeedLinear = synthState_->fullState.masterfxConfig[GLOBALFX_LFOSPEED];
	if (prevLfoSpeedLinear == lfoSpeedLinear)
	{
		temp = lfoSpeedLinear;
		temp *= temp * temp;
		lfoSpeed = lfoSpeed * 0.9f + temp * 0.1f;
	}
	prevLfoSpeedLinear = lfoSpeedLinear;

	temp = synthState_->fullState.masterfxConfig[GLOBALFX_LFODEPTH];
	temp = temp * (1 - lfoSpeedLinear * 0.5f) * (1 - sizeParam * 0.5f);
	lfoDepth = lfoDepth * 0.9f + temp * 0.1f;

	// ------ page 2

    nextSizeParam 	= clamp(synthState_->fullState.masterfxConfig[GLOBALFX_SIZE], 0.03f, 1);
	sizeParam 		= sizeParam * 0.99f + nextSizeParam * 0.01f;

	diffusion 		= 	synthState_->fullState.masterfxConfig[GLOBALFX_DIFFUSION];
	damping 		= 	synthState_->fullState.masterfxConfig[GLOBALFX_DAMPING] * 0.76f;
	damping 		*= 	damping;
	loopLpf 		= 	0.05f + damping * 0.95f;
    decayVal 		= 	synthState_->fullState.masterfxConfig[ GLOBALFX_DECAY ];
    notchBase 		= 	synthState_->fullState.masterfxConfig[ GLOBALFX_NOTCHBASE ];

    //------- some process

    if(diffusion != prevDiffusion) {
		diffuserCoef1 	= 	-(0.01f + diffusion * 0.67f);
		diffuserCoef2 	= 	-(0.01f + diffusion * 0.5f);
		prevDiffusion = diffusion;
    }

    if(prevDecayVal != decayVal) {
    	decayFdbck =  fastroot(decayVal, 3) * decayMaxVal;

    	float decayValSquare = decayVal * decayVal;
        envRelease 			= 	0.005f + decayValSquare * 0.6f;

    	headRoomMultiplier = (1 + (1 - decayValSquare) * 0.75f) * 4 * 0.65f;
    	headRoomDivider = 0.125f;
        sampleMultipler = headRoomMultiplier * (float) 0x7fffff;
    }
    prevDecayVal = decayVal;

    notchSpread = synthState_->fullState.masterfxConfig[GLOBALFX_NOTCHSPREAD];
    if(notchBase != prevNotchBase || notchSpread != prevNotchSpread) {
        float offset = notchSpread * notchSpread * 0.17f;
        float range = 2;

        const float windowMin = 0.005f, windowMax = 0.99f;

        f1L = clamp(fold((notchBase - offset) * range) * 0.5f, windowMin, windowMax);
        f2L = clamp(fold((notchBase + offset) * range) * 0.5f, windowMin, windowMax);
        f3L = clamp(fold((notchBase - (offset * 2)) * range) * 0.5f, windowMin, 1);
        f4L = clamp(fold((notchBase + (offset * 2)) * range) * 0.5f, windowMin, 1);

        coef1L = (1.0f - f1L) / (1.0f + f1L);
        coef2L = (1.0f - f2L) / (1.0f + f2L);
        coef3L = (1.0f - f3L) / (1.0f + f3L);
        coef4L = (1.0f - f4L) / (1.0f + f4L);
    }
    prevNotchBase = notchBase;
    prevNotchSpread = notchSpread;

	if(sizeParam != prevSizeParam) {
		sizeSqrt = sqrt3(sizeParam);

		sizeParamInpt = 0.1f + sizeSqrt * 0.9f;
		inputBuffer1ReadLen = inputBufferLen1 * sizeParamInpt;
		inputBuffer2ReadLen = inputBufferLen2 * sizeParamInpt;
		inputBuffer3ReadLen = inputBufferLen3 * sizeParamInpt;
		inputBuffer4ReadLen = inputBufferLen4 * sizeParamInpt;

		delay1ReadLen = 1 + delay1BufferSizeM1 * (1 - sizeParam);
		delay2ReadLen = 1 + delay2BufferSizeM1 * (1 - sizeParam);
		delay3ReadLen = 1 + delay3BufferSizeM1 * (1 - sizeParam);
		delay4ReadLen = 1 + delay4BufferSizeM1 * (1 - sizeParam);

		diffuserBuffer1ReadLen = -1 + diffuserBufferLen1M1 * sizeParam;
		diffuserBuffer2ReadLen = -1 + diffuserBufferLen2M1 * sizeParam;
		diffuserBuffer3ReadLen = -1 + diffuserBufferLen3M1 * sizeParam;
		diffuserBuffer4ReadLen = -1 + diffuserBufferLen4M1 * sizeParam;

		diffuserBuffer1ReadLen_b = diffuserBuffer1ReadLen * (0.1333f 	* 1.5f);
		diffuserBuffer2ReadLen_b = diffuserBuffer2ReadLen * (0.2237f 	* 1.5f);
		diffuserBuffer3ReadLen_b = diffuserBuffer3ReadLen * (0.18f 		* 1.5f);
		diffuserBuffer4ReadLen_b = diffuserBuffer4ReadLen * (0.4237f 	* 1.5f);
	}
	prevSizeParam = sizeParam;

	//

	totalSent = 0;
}

/**
 * add timbre block to bus mix
 */
void FxBus::mixAdd(float *inStereo, int timbreNum) {

	if(synthState_->mixerState.instrumentState_[timbreNum].send > 0) {

		const float level = fastroot(synthState_->mixerState.instrumentState_[timbreNum].send, 3) * headRoomDivider;
		totalSent += synthState_->mixerState.instrumentState_[timbreNum].send;

		sample = getSampleBlock();
		for (int s = 0; s < 8; s++) {
			*(sample++) += *inStereo++ * level;
			*(sample++) += *inStereo++ * level;
			*(sample++) += *inStereo++ * level;
			*(sample++) += *inStereo++ * level;
			*(sample++) += *inStereo++ * level;
			*(sample++) += *inStereo++ * level;
			*(sample++) += *inStereo++ * level;
			*(sample++) += *inStereo++ * level;
		}

	}
}

/**
 * process fx on bus mix
 */
void FxBus::processBlock(int32_t *outBuff) {
	sample = getSampleBlock();

	isActive = totalSent > 0 || silentBlockCount < 64;

	if(!isActive) {
		return;
	}

	float outSum = 0;
	for (int s = 0; s < BLOCK_SIZE; s++) {

        // --- audio in

    	inR = *(sample);
    	inL = *(sample + 1);

    	monoIn = (inR + inL);

		// allpass / notch

    	lowL = coef1L * (lowL + monoIn) - bandL;
    	bandL = monoIn;
    	lowL2 = coef2L * (lowL2 + lowL) - bandL2;
    	bandL2 = lowL;
    	lowL3 = coef3L * (lowL3 + lowL2) - bandL3;
    	bandL3 = lowL2;
        lowL4 = coef4L * (lowL4 + lowL3) - bandL4;
        bandL4 = lowL3;

        monoIn += lowL4;

		dcBlock4a = monoIn - dcBlock4b + dcBlockerCoef1 * dcBlock4a;			// dc blocker
		dcBlock4b = monoIn;

        monoIn = dcBlock4a;

        // --- cut high

        v6R += inLpF * v7R;						// lowpass
        v7R += inLpF * (monoIn - v6R - v7R);
        v6L += inLpF * v7L;						// lowpass
        v7L += inLpF * (v6R - v6L - v7L);

        monoIn = v6L;

        // --- hi pass

        v0R += inHpf * v1R;						// hipass
        v1R += inHpf * ( monoIn - v0R - v1R);
        v0L += inHpf * v1L;						// hipass
        v1L += inHpf * ( v0R - v0L - v1L);

        monoIn -= v0L;


    	//--- pre delay

    	predelayBuffer[predelayWritePos] = monoIn;

    	predelayReadPos = modulo2(predelayWritePos - predelaySize, predelayBufferSize);
    	preDelayOut = delayInterpolation(predelayReadPos, predelayBuffer, predelayBufferSizeM1);
    	monoIn = predelayMixAttn * preDelayOut + (1 - predelayMixAttn) * monoIn;

    	// --- input diffuser

        // ---- diffuser 1

        inputReadPos1 	= modulo2(inputWritePos1 - inputBuffer1ReadLen, inputBufferLen1);
        float in_apSum1 = monoIn + inputBuffer1[inputReadPos1] * inputCoef1;
        diff1Out 		= monoIn - in_apSum1 * inputCoef1;
        inputBuffer1[inputWritePos1] 		= in_apSum1;

        // ---- diffuser 2

    	inputReadPos2 	= modulo2(inputWritePos2 - inputBuffer2ReadLen, inputBufferLen2);
        float in_apSum2 = diff1Out + inputBuffer2[inputReadPos2] * inputCoef1;
        diff2Out 		= diff1Out - in_apSum2 * inputCoef1;
        inputBuffer2[inputWritePos2] 		= in_apSum2;

        // ---- diffuser 3

    	inputReadPos3 	= modulo2(inputWritePos3 - inputBuffer3ReadLen, inputBufferLen3);
        float in_apSum3 = diff2Out + inputBuffer3[inputReadPos3] * inputCoef2;
        diff3Out 		= diff2Out - in_apSum3 * inputCoef2;
        inputBuffer3[inputWritePos3] 		= in_apSum3;

        // ---- diffuser 4

    	inputReadPos4 	= modulo2(inputWritePos4 - inputBuffer4ReadLen, inputBufferLen4);
        float in_apSum4 = diff3Out + inputBuffer4[inputReadPos4] * inputCoef2;
        diff4Out 		= diff3Out - in_apSum4 * inputCoef2;
        inputBuffer4[inputWritePos4] 		= in_apSum4;

		// ---- dc blocker

		dcBlock1a = diff4Out - dcBlock1b + dcBlockerCoef2 * dcBlock1a;
		dcBlock1b = diff4Out;

		monoIn = dcBlock1a;

        // ---- ap 1

        ap1In = monoIn + feedbackInL * decayFdbck;

		diffuserReadPos1 = modulo2(diffuserWritePos1 - timeCvControl1, diffuserBufferLen1);

		float inSum1 	= 	ap1In + int1 * diffuserCoef1;
    	ap1Out 			= 	int1 - inSum1 * diffuserCoef1;
        diffuserBuffer1[diffuserWritePos1] 		= inSum1;
    	int1 = delayAllpassInterpolation(diffuserReadPos1, diffuserBuffer1, diffuserBufferLen1M1, int1);

    	// ----------------------------------------------< inject in delay1
    	delay1Buffer[ delay1WritePos ] 		= ap1Out;
    	// ----------------------------------------------> read delay1
		ap2In = delay1Buffer[ delay1ReadPos ];

        // ---------------------------------------------------- filter

        v4R += loopLpf * v5R;						// lowpass
        v5R += loopLpf * ( ap2In - v4R - v5R);

		dcBlock3a = v4R - dcBlock3b + dcBlockerCoef1 * dcBlock3a;			// dc blocker
		dcBlock3b = v4R;

    	//dcBlock5a = dcBlock3a - dcBlock5b + dcBlockerCoef3 * dcBlock5a;			// dc blocker
    	//dcBlock5b = dcBlock3a;

		ap2In = dcBlock3a * decayFdbck;

    	// ---- ap 2

		diffuserReadPos2 = modulo2(diffuserWritePos2 - timeCvControl2, diffuserBufferLen2);

		float inSum2 	= 	ap2In + int2 * diffuserCoef2;
    	ap2Out 			= 	int2 - inSum2 * diffuserCoef2;
        diffuserBuffer2[diffuserWritePos2] 		= inSum2;
    	int2 = delayAllpassInterpolation(diffuserReadPos2, diffuserBuffer2, diffuserBufferLen2M1, int2);

    	// ----------------------------------------------< inject in delay2
    	delay2Buffer[ delay2WritePos ] 		= ap2Out;
    	// ----------------------------------------------> read delay2

    	ap2Out = delay2Buffer[ delay2ReadPos ];

        // ---- ap 3

		ap3In = monoIn + feedbackInR * decayFdbck;

		diffuserReadPos3 = modulo2(diffuserWritePos3 - timeCvControl3, diffuserBufferLen3);

		float inSum3 	= 	ap3In + int3 * diffuserCoef1;
    	ap3Out 			= 	int3 - inSum3 * diffuserCoef1;
        diffuserBuffer3[diffuserWritePos3] 		= inSum3;
    	int3 = delayAllpassInterpolation(diffuserReadPos3, diffuserBuffer3, diffuserBufferLen3M1, int3);

    	// ----------------------------------------------< inject in delay3
    	delay3Buffer[ delay3WritePos ] 		= ap3Out;
    	// ----------------------------------------------> read delay3
		ap4In = delay3Buffer[ delay3ReadPos ];

        // ---------------------------------------------------- filter

        v4L += loopLpf * v5L;						// lowpass
        v5L += loopLpf * ( ap4In - v4L - v5L);

        //v4L += loopLpf * v5L;						// lowpass
        //v5L += loopLpf * ( ap4In - v4L - v5L);

		dcBlock2a = v4L - dcBlock2b + dcBlockerCoef2 * dcBlock2a;			// dc blocker
		dcBlock2b = v4L;

        ap4In = dcBlock2a * decayFdbck;				// decay

        // ---- ap 4

        diffuserReadPos4 = modulo2(diffuserWritePos4 - timeCvControl4, diffuserBufferLen4);

    	float inSum4 	= 	ap4In + int4 * diffuserCoef2;
    	ap4Out 			= 	int4 - inSum4 * diffuserCoef2;
        diffuserBuffer4[diffuserWritePos4] 		= inSum4;
       	int4 = delayAllpassInterpolation(diffuserReadPos4, diffuserBuffer4, diffuserBufferLen4M1, int4);

    	// ----------------------------------------------< inject in delay4
    	delay4Buffer[ delay4WritePos ] 		= ap4Out;
    	// ----------------------------------------------> read delay4
    	ap4Out = delay4Buffer[ delay4ReadPos ];

    	// ================================================  mix out

    	feedbackInL = ap4Out;
    	feedbackInR = ap2Out;

        outL = ap1Out;
        outL += delay1Buffer[ 		modulo(_kLeftTaps[0] + delay1WritePos, 		delay1BufferSize)		];
        outL += delay1Buffer[ 		modulo(_kLeftTaps[1] + delay1WritePos, 		delay1BufferSize)		];
        outL -= diffuserBuffer2[ 	modulo(_kLeftTaps[2] + diffuserWritePos2, 	diffuserBufferLen2)		];
        outL += delay2Buffer[ 		modulo(_kLeftTaps[3] + delay2WritePos, 		delay2BufferSize)		];
        outL -= delay3Buffer[ 		modulo(_kLeftTaps[4] + delay3WritePos, 		delay3BufferSize)		];
        outL -= diffuserBuffer4[ 	modulo(_kLeftTaps[5] + diffuserWritePos4, 	diffuserBufferLen4)		];
        outL -= delay4Buffer[ 		modulo(_kLeftTaps[6] + delay4WritePos, 		delay4BufferSize)		];

        outR = ap3Out;
        outR += delay3Buffer[ 		modulo(_kRightTaps[0] + delay3WritePos, 	delay3BufferSize) 		];
        outR += delay3Buffer[ 		modulo(_kRightTaps[1] + delay3WritePos, 	delay3BufferSize) 		];
        outR -= diffuserBuffer4[ 	modulo(_kRightTaps[2] + diffuserWritePos4, 	diffuserBufferLen4) 	];
        outR += delay4Buffer[ 		modulo(_kRightTaps[3] + delay4WritePos, 	delay4BufferSize)	 	];
        outR -= delay1Buffer[ 		modulo(_kRightTaps[4] + delay1WritePos, 	delay1BufferSize)	 	];
        outR -= diffuserBuffer2[ 	modulo(_kRightTaps[5] + diffuserWritePos2, 	diffuserBufferLen2) 	];
        outR -= delay2Buffer[ 		modulo(_kRightTaps[6] + delay2WritePos, 	delay2BufferSize)	 	];

    	*(outBuff++) += (int32_t) ( outL * sampleMultipler);
    	*(outBuff++) += (int32_t) ( outR * sampleMultipler);

    	outSum += fabsf(outL + outR);
    	// ================================================ index increment

    	sample += 2;

    	inputWritePos1		= modulo(inputWritePos1 + 1 , inputBufferLen1);
    	inputWritePos2		= modulo(inputWritePos2 + 1 , inputBufferLen2);
    	inputWritePos3		= modulo(inputWritePos3 + 1 , inputBufferLen3);
    	inputWritePos4		= modulo(inputWritePos4 + 1 , inputBufferLen4);

    	predelayWritePos	= modulo(predelayWritePos + 1 , predelayBufferSize);

        diffuserWritePos1	= modulo(diffuserWritePos1 + 1 , diffuserBufferLen1);
        diffuserWritePos2	= modulo(diffuserWritePos2 + 1 , diffuserBufferLen2);
        diffuserWritePos3	= modulo(diffuserWritePos3 + 1 , diffuserBufferLen3);
        diffuserWritePos4	= modulo(diffuserWritePos4 + 1 , diffuserBufferLen4);

        delay1WritePos		= modulo(delay1WritePos + 1 , delay1BufferSize);
        delay2WritePos		= modulo(delay2WritePos + 1 , delay2BufferSize);
        delay3WritePos		= modulo(delay3WritePos + 1 , delay3BufferSize);
        delay4WritePos		= modulo(delay4WritePos + 1 , delay4BufferSize);

    	delay1ReadPos 		= modulo(delay1WritePos + delay1ReadLen, delay1BufferSize);
    	delay2ReadPos 		= modulo(delay2WritePos + delay2ReadLen, delay2BufferSize);
    	delay3ReadPos 		= modulo(delay3WritePos + delay3ReadLen, delay3BufferSize);
    	delay4ReadPos 		= modulo(delay4WritePos + delay4ReadLen, delay4BufferSize);

    	// --- lfo increment

    	lfoProcess(&lfo1, &lfo1tri, &lfo1Inc);
    	lfoProcess(&lfo2, &lfo2tri, &lfo2Inc);
    	lfoProcess(&lfo3, &lfo3tri, &lfo3Inc);
    	lfoProcess(&lfo4, &lfo4tri, &lfo4Inc);

     	// -------- mods :

    	timeCvControl1 = clamp(diffuserBufferLen1 - diffuserBuffer1ReadLen  	+  	lfo1  	* lfoDepth * diffuserBuffer1ReadLen_b, -diffuserBufferLen1M1, diffuserBufferLen1);
    	timeCvControl2 = clamp(diffuserBufferLen2 - diffuserBuffer2ReadLen  	+ 	lfo2 	* lfoDepth * diffuserBuffer2ReadLen_b, -diffuserBufferLen2M1, diffuserBufferLen2);
    	timeCvControl3 = clamp(diffuserBufferLen3 - diffuserBuffer3ReadLen  	+  	lfo3	* lfoDepth * diffuserBuffer3ReadLen_b, -diffuserBufferLen3M1, diffuserBufferLen3);
    	timeCvControl4 = clamp(diffuserBufferLen4 - diffuserBuffer4ReadLen  	+ 	lfo4	* lfoDepth * diffuserBuffer4ReadLen_b, -diffuserBufferLen4M1, diffuserBufferLen4);
    }

	silentBlockCount = (outSum < 0.00001f) ? silentBlockCount + 1 : 0;
}
float FxBus::delayAllpassInterpolation(float readPos, float buffer[], int bufferLenM1, float prevVal) {
	//v[n] = VoiceL[i + 1] + (1 - frac)  * VoiceL[i] - (1 - frac)  * v[n - 1]
	int readPosInt = readPos;
	float y0 = buffer[readPosInt];
	float y1 = buffer[(unlikely(readPosInt >= bufferLenM1) ? readPosInt - bufferLenM1 + 1 : readPosInt + 1)];
	//float y1 = buffer[((readPosInt == 0 ) ? bufferLenM1: readPosInt - 1)];
	float x = readPos - floorf(readPos);
    return y1 + (1 - x) * (y0 - prevVal);
}
float FxBus::delayInterpolation(float readPos, float buffer[], int bufferLenM1) {
	int readPosInt = readPos;
	float y0 = buffer[readPosInt];
	float y1 = buffer[(unlikely(readPosInt == 0) ? bufferLenM1 : readPosInt - 1)];
	float x = readPos - floorf(readPos);
    return y0 + x * (y1 - y0);
}

void FxBus::lfoProcess(float *lfo, float *lfotri, float *lfoInc) {
	*lfotri += *lfoInc * lfoSpeed;
	if(unlikely(*lfotri >= 1)) {
		*lfotri = 1;
		*lfoInc = -*lfoInc;
	}
	if(unlikely(*lfotri <= 0)) {
		*lfotri = 0;
		*lfoInc = -*lfoInc;
	}
	*lfo = (*lfo * lfoLpCoef1 + *lfotri ) * lfoLpCoef2;
}
