/*
 * FxBus.cpp
 *
 *  Created on: Feb 14, 2021
 *      Author: patvig
 */

#include "FxBus.h"

#define filterWindowMin 0.01f
#define filterWindowMax 0.99f

extern float diatonicScaleFrequency[];
extern float sinTable[];

inline
float clamp(float d, float min, float max) {
  const float t = unlikely(d < min) ? min : d;
  return unlikely(t > max) ? max : t;
}
inline float getQuantizedTime(float t, float maxSize)
{
	return t * maxSize;

	/*int q = t * 127;
	return clamp((PREENFM_FREQUENCY / diatonicScaleFrequency[127 - q]), 0, maxSize);*/
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
// laurent de soras
inline float hermite4(float frac_pos, float xm1, float x0, float x1, float x2)
{
   const float    c     = (x1 - xm1) * 0.5f;
   const float    v     = x0 - x1;
   const float    w     = c + v;
   const float    a     = w + v + (x2 - x0) * 0.5f;
   const float    b_neg = w + a;

   return ((((a * frac_pos) - b_neg) * frac_pos + c) * frac_pos + x0);
}
inline
float sigmoid(float x)
{
    return x * (1.5f - 0.5f * x * x);
}
inline
float tri2sin(float x) {
	return x * x * (3 - x - x);
}
//***------------***------------***------------***------------***----- FxBus -------***------------***------------***------------

float FxBus::delay1Buffer[delay1BufferSize] __attribute__((section(".ram_d1")));
float FxBus::delay2Buffer[delay2BufferSize] __attribute__((section(".ram_d1")));
float FxBus::delay3Buffer[delay3BufferSize] __attribute__((section(".ram_d1")));
float FxBus::delay4Buffer[delay4BufferSize] __attribute__((section(".ram_d1")));
float FxBus::tapDelayBuffer[tapDelayBufferSize] __attribute__((section(".ram_d1")));
int FxBus::tapDelayPos[tapCount] __attribute__((section(".ram_d1")));
float FxBus::tapDelayAmp[tapCount] __attribute__((section(".ram_d1")));

float FxBus::predelayBuffer[predelayBufferSize] __attribute__((section(".ram_d2")));

float FxBus::inputBuffer1[inputBufferLen1] __attribute__((section(".ram_d2")));
float FxBus::inputBuffer2[inputBufferLen2] __attribute__((section(".ram_d2")));
float FxBus::inputBuffer3[inputBufferLen3] __attribute__((section(".ram_d2")));
float FxBus::inputBuffer4[inputBufferLen4] __attribute__((section(".ram_d2")));

float FxBus::diffuserBuffer1[diffuserBufferLen1] __attribute__((section(".ram_d2")));
float FxBus::diffuserBuffer2[diffuserBufferLen2] __attribute__((section(".ram_d2")));
float FxBus::diffuserBuffer3[diffuserBufferLen3] __attribute__((section(".ram_d2")));
float FxBus::diffuserBuffer4[diffuserBufferLen4] __attribute__((section(".ram_d2")));


extern float noise[32];

FxBus::FxBus() {
	lfo1 = 0;
}

void FxBus::init(SynthState *synthState) {
    this->synthState_ = synthState;

	for (int s = 0; s < predelayBufferSize; s++) {
    	predelayBuffer[s] = 0;
    }
	for (int s = 0; s < tapDelayBufferSize; s++) {
		tapDelayBuffer[ s ] = 0;
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

	_ly1L = 0;
	_ly1R = 0;
	_ly2L = 0;
	_ly2R = 0;
	_ly3L = 0;
	_ly3R = 0;
	_ly4L = 0;
	_ly4R = 0;
	_lx1L = 0;
	_lx1R = 0;
	_lx2L = 0;
	_lx2R = 0;
	_lx3L = 0;
	_lx3R = 0;
	_lx4L = 0;
	_lx4R = 0;

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

	monoInHpf = 0.05f;
	loopLpf = 0.43f;
	loopLpf2 = loopLpf - 0.02f;
	loopHpf = 0.1f;
	harmTremoloCutF = 0.58f;

	inLpF = 0.55f;

	diffuserCoef1b =  diffuserCoef1 * (1  - (diffuserCoef1 * diffuserCoef1));
	diffuserCoef2b =  diffuserCoef2 * (1  - (diffuserCoef2 * diffuserCoef2));

	inputCoef1b =  (1  - (inputCoef1 * inputCoef1));
	inputCoef2b =  (1  - (inputCoef2 * inputCoef2));

	tapDelayPos[0] 	= 187;
	tapDelayPos[1] 	= 266;
	tapDelayPos[2] 	= 1066;
	tapDelayPos[3] 	= 1913;
	tapDelayPos[4] 	= 1996;
	tapDelayPos[5] 	= 1990;
	tapDelayPos[6] 	= 2974;

	tapDelayAmp[0] 	= 0.6f;
	tapDelayAmp[1] 	= 0.45f;
	tapDelayAmp[2] 	= -0.4f;
	tapDelayAmp[3] 	= 0.3f;
	tapDelayAmp[4] 	= -0.2f;
	tapDelayAmp[5] 	= -0.1f;
	tapDelayAmp[6] 	= -0.05f;

}

/**
 * init before timbres summing
 */
void FxBus::mixSumInit() {
    float temp;
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

    // ------ page 1

    fxTimeLinear	= 	synthState_->fullState.masterfxConfig[ GLOBALFX_TIME ];
	prevTime 		=  	clamp( fxTimeLinear, 0.0003f, 0.9996f);
	prevTime		*= 	prevTime * prevTime;
	fxTime 			= 	fxTime * 0.9f + prevTime * 0.1f;

	delay1FxTarget = 	getQuantizedTime(fxTime, delay1BufferSize);
	delay1DelayLen = 	(delay1DelayLen 	+ (delay1FxTarget - delay1DelayLen)	* 0.01f);

	delay3FxTarget = 	delay3BufferSizeM1;
	delay3DelayLen = 	(delay3DelayLen 	+ (delay3FxTarget - delay3DelayLen)	* 0.01f);

	predelaySize 	= 	fxTimeLinear * predelayBufferSizeM1;

	prevSpeed 		= 	synthState_->fullState.masterfxConfig[ GLOBALFX_LFOSPEED ];

	temp 			= 	prevSpeed;
	temp 			*=	temp * temp;
	fxSpeed 		= 	fxSpeed * 0.9f + temp * 0.1f;

    feedbackGain 		= 	synthState_->fullState.masterfxConfig[ GLOBALFX_FBACK ];

    predelayMixLevel 	= 	synthState_->fullState.masterfxConfig[ GLOBALFX_PREDELAYMIX ];
    predelayMixAttn 	= 	predelayMixLevel * (1 - (predelayMixLevel * predelayMixLevel * 0.1f));

	prevEnvThreshold = 	synthState_->fullState.masterfxConfig[ GLOBALFX_ENVTHRESHOLD ];
	envThreshold	= 	envThreshold * 0.9f + prevEnvThreshold * 0.1f;

	prevBounce 		= 	synthState_->fullState.masterfxConfig[ GLOBALFX_BOUNCE ] + 0.005f ;
	bounceLevel	 	= 	bounceLevel * 0.9f + prevBounce * 0.1f;

	prevEnvRelease 	= 	synthState_->fullState.masterfxConfig[ GLOBALFX_ENVRELEASE ] + 0.005f ;
	prevEnvRelease	*= 	prevEnvRelease * prevEnvRelease;
	envRelease	 	= 	envRelease * 0.9f + prevEnvRelease * 0.1f;

    // ------ page 2

    temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_LFODEPTH ];
	invspeed 		= 	1 - sqrt3(fxTimeLinear);
    temp 			= 	temp * temp * (0.01f + invspeed * 0.9f);
    lfoDepth 		= 	lfoDepth * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_ENVMOD] * 0.5f;
    temp 			*= 	temp * (0.4f + invspeed * 0.6f);
    envModDepth 	= 	envModDepth * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_TREMOLOSPEED];
	fxTremoloSpeed	= 	fxTremoloSpeed * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_TREMOLODEPTH ];
	invspeed 		= 	1 - sqrt3(fxTremoloSpeed);
    temp 			= 	temp * (0.8f + invspeed * 0.2f) ;
	fxTremoloDepth	= 	fxTremoloDepth * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_TREMOLOENVFOLLOW] * 0.99f;
    temp 			= 	temp * (0.8f + invspeed * 0.2f);
    tremoloEnvFollow= 	tremoloEnvFollow * 0.9f + temp * 0.1f;
    tremoloEnvFollowAbs = fabsf(tremoloEnvFollow);

    // ------ env follow

	envBlocknn += BLOCK_SIZE;

	if(envBlocknn > envDetectSize) {
		if( blocksum > envThreshold * 0.5f) {
			// attack
			envDest = 1;
			envM1 = 249;
			envM2 = 0.004f;
			//envelope = 0;// restart env
		} else if(envDest == 1) {
			// release
			envDest = 0;
			envM1 = envRelease * 800000;
			envM2 = (1 / (envM1 + 1));
		}
	    blocksum = 0;
	    envBlocknn = 0;
	}

	if (lfo2ChangeCounter++ > lfo2ChangePeriod) {
		float nextLfo2ModVal = noise[0] * noise[0] * 0.125f;
		float delta = nextLfo2ModVal - lfo2ModVal;
		lfo2ChangeCounter = 0;
		lfo2ModVal = nextLfo2ModVal;
		lfo2IncModSampleInc = delta * lfo2ChangePeriodInv;
	}

	if (loopDecouplerChangeCounter++ > loopDecouplerChangePeriod) {
		loopDecouplerChangeCounter = 0;
		loopDecouplerModVal = 0.8f + (noise[0] * noise[0]) * 0.2f;
		loopDecoupler = decoupler1 * loopDecouplerModVal;
		loopDecoupler2 = decoupler2 * loopDecouplerModVal;
	}

	inLpF = 0.5f - (fxTimeLinear + fxTimeLinear * predelayMixLevel) * 0.15f;

}

/**
 * add timbre block to bus mix
 */
void FxBus::mixAdd(float *inStereo, int timbreNum, int32_t *outBuffer) {
	const float level = synthState_->mixerState.instrumentState_[timbreNum].send * 0.01f; // divide for more headroom

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

/**
 * process fx on bus mix
 */
void FxBus::processBlock(int32_t *outBuff) {
	sample = getSampleBlock();
    const float sampleMultipler = 80 * (float) 0x7fffff; // fx level

	float tremoloEnvFollowMod, tremoloEnvFollowModAttn;
	float tremoloMod = 0, inLpMod;
	float tapAccumulator;
	float decouple1, decouple2, decouple3;
	float loopOutL, loopOutR;
	int tapPos;

	for (int s = 0; s < BLOCK_SIZE; s++) {

		decouple1 = lfo2 * loopDecoupler;
		decouple2 = lfo2b * loopDecoupler2;
		//decouple3 = lfo2b * loopDecoupler2;

    	delay1ReadPos = delay1WritePos + 1;
    	while( unlikely(delay1ReadPos < 0) )
    		delay1ReadPos += delay1BufferSize;
    	while( unlikely(delay1ReadPos >= delay1BufferSizeM1) )
    		delay1ReadPos -= delay1BufferSizeM1;


    	delay2ReadPos = delay2WritePos + 1;
    	if( unlikely(delay2ReadPos < 0) )
    		delay2ReadPos += delay2BufferSize;
    	if( unlikely(delay2ReadPos >= delay2BufferSizeM1) )
    		delay2ReadPos -= delay2BufferSizeM1;

    	delay3ReadPos = delay3WritePos + 1;
    	while( unlikely(delay3ReadPos < 0) )
    		delay3ReadPos += delay3BufferSize;
    	while( unlikely(delay3ReadPos >= delay3BufferSizeM1) )
    		delay3ReadPos -= delay3BufferSizeM1;

    	delay4ReadPos = delay4WritePos + 1;
    	if( unlikely(delay4ReadPos < 0) )
    		delay4ReadPos += delay4BufferSize;
    	if( unlikely(delay4ReadPos >= delay4BufferSizeM1) )
    		delay4ReadPos -= delay4BufferSizeM1;

        // --- audio in

    	inR = *(sample);
    	inL = *(sample + 1);

    	monoIn = (inR + inL) * 0.5f;

    	//--- pre delay

    	predelayReadPos = predelayWritePos - predelaySize;
    	if( unlikely(predelayReadPos < 0) )
    		predelayReadPos += predelayBufferSize;
    	if( unlikely(predelayReadPos >= predelayBufferSizeM1) )
    		predelayReadPos -= predelayBufferSizeM1;

    	predelayBuffer[predelayWritePos] = monoIn;

    	//crossFade()
    	monoIn = predelayMixAttn * predelayInterpolation(predelayReadPos) + (1-predelayMixAttn) * monoIn;


        // --- cut high

    	inLpMod = inLpF + envelope * 0.05f;

        v6L += inLpMod * v7L;
        v7L += inLpMod * (monoIn - v6L - v7L);
        v6L += inLpMod * v7L;
        v7L += inLpMod * (monoIn - v6L - v7L);

        monoIn = v6L;
    	// --- audio in > harmonic tremolo

        tremoloEnvFollowMod 	= 1 + envelope * tremoloEnvFollow;
        tremoloEnvFollowModAttn = (tremoloEnvFollowMod + (1 - tremoloEnvFollowAbs));

        _lx3L += harmTremoloCutF * _ly3L; // low pass filter
        _ly3L += harmTremoloCutF * ((monoIn) - _lx3L - _ly3L);

        lpL = _lx3L;
        hpL = inL - lpL;

    	tremoloMod = ((lfoTremoloSin 		* fxTremoloDepth + 1 - fxTremoloDepth) * 2 - 1) * tremoloEnvFollowModAttn;
    	monoIn = lpL * tremoloMod + hpL * (1 - tremoloMod);

    	// --- enveloppe calculation

        blocksum 	+= fabsf(monoIn);
        envelope 	= (envelope * envM1 + envDest) * envM2;
        envMod 		= envModDepth * envelope;

    	// initial tap

        tapDelayBuffer[ tapDelayWritePos++ ]  = monoIn;

        if(tapDelayWritePos >= tapDelayBufferSize) {
        	tapDelayWritePos = 0;
        }

        tapAccumulator = monoIn;

        for(int aa=0; aa<tapCount; aa++) {
            tapPos = (tapDelayWritePos + tapDelayPos[aa]);
            if(tapPos >= tapDelayBufferSize)
            	tapPos -= tapDelayBufferSize;
            tapAccumulator += tapDelayBuffer[ tapPos ] * tapDelayAmp[aa];
        }

        monoIn = (tapAccumulator);

        // --- hi pass

        v0R += monoInHpf * v1R;						// hipass
        v1R += monoInHpf * ( monoIn - v0R - v1R);
        monoIn -= v0R;

    	// --- input diffuser

        // ---- diffuser 1

        inputReadPos1 = inputWritePos1 + 1;
    	if( unlikely(inputReadPos1 >= inputBufferLen1M1) )
    		inputReadPos1 -= inputBufferLen1M1;
    	inputBuffer1[inputWritePos1] = monoIn - inputCoef1 * inputBuffer1[inputReadPos1];
    	diff1Out = monoIn * inputCoef1 + 	inputBuffer1[inputReadPos1] * inputCoef1b;

        // ---- diffuser 2

    	inputReadPos2 = inputWritePos2 + 1;
    	if( unlikely(inputReadPos2 >= inputBufferLen2M1) )
    		inputReadPos2 -= inputBufferLen2M1;
    	inputBuffer2[inputWritePos2] = diff1Out - inputCoef1 * inputBuffer2[inputReadPos2];
    	diff2Out = diff1Out * inputCoef1 + 	inputBuffer2[inputReadPos2] * inputCoef1b;

        // ---- diffuser 3

    	inputReadPos3 = inputWritePos3 + 1;
    	if( unlikely(inputReadPos3 >= inputBufferLen3M1) )
    		inputReadPos3 -= inputBufferLen3M1;
    	inputBuffer3[inputWritePos3] = diff2Out - inputCoef2 * inputBuffer3[inputReadPos3];
    	diff3Out = diff2Out * inputCoef2 + 	inputBuffer3[inputReadPos3] * inputCoef2b;

        // ---- diffuser 4

    	inputReadPos4 = inputWritePos4 + 1;
    	if( unlikely(inputReadPos4 >= inputBufferLen4M1) )
    		inputReadPos4 -= inputBufferLen4M1;
    	inputBuffer4[inputWritePos4] = diff3Out - inputCoef2 * inputBuffer4[inputReadPos4];
    	diff4Out = diff3Out * inputCoef2 + 	inputBuffer4[inputReadPos4] * inputCoef2b;

    	monoIn = diff4Out;

    	float envFeedbackMod = feedbackGain * clamp((1 + envMod), -1, 1);
        // ---- ap 1

        ap1In = monoIn + feedbackInL * envFeedbackMod;

        diffuserReadPos1 = diffuserWritePos1 - (diffuserBufferLen1 + timeCvControl + lfo2 * 24);
    	while( (diffuserReadPos1 < 0) )
    		diffuserReadPos1 += diffuserBufferLen1;
    	while( (diffuserReadPos1 >= diffuserBufferLen1M1) )
    		diffuserReadPos1 -= diffuserBufferLen1M1;

    	float int1 = diffuser1Interpolation(diffuserReadPos1);
    	//float int1 = diffuser1HermiteInterpolation((int) diffuserReadPos1);

        diffuserBuffer1[diffuserWritePos1] 		= ap1In + diffuserCoef1 * int1;
        ap1Out = -ap1In * diffuserCoef1 + 	int1 * diffuserCoef1b;

    	// ----------------------------------------------< inject in delay1
    	delay1Buffer[ delay1WritePos ] 		= ap1Out;


    	// ----------------------------------------------> read delay1
		ap2In = delay1Interpolation(delay1ReadPos);

        // --- filter

        v2L += loopHpf * v3L;						// hipass
        v3L += loopHpf * ( ap2In - v2L - v3L);
        v4L += loopLpf * v5L;						// lowpass
        v5L += loopLpf * ( (ap2In - v3L) - v4L - v5L);

        ap2In = v4L * envFeedbackMod;				// decay

    	// ---- ap 2

        diffuserReadPos2 = diffuserWritePos2 + 1;
    	if( diffuserReadPos2 >= diffuserBufferLen2M1 )
    		diffuserReadPos2 -= diffuserBufferLen2M1;

        diffuserBuffer2[diffuserWritePos2] 		= ap2In - diffuserCoef2 * diffuserBuffer2[diffuserReadPos2];
        ap2Out = ap2In * diffuserCoef2 + 	diffuser2Interpolation(diffuserReadPos2) * diffuserCoef2b;

    	// ----------------------------------------------< inject in delay2
    	delay2Buffer[ delay2WritePos ] 		= ap2Out;

    	// ----------------------------------------------> read delay2
    	ap2Out = delay2Interpolation(delay2ReadPos);

        // ---- ap 3

		ap3In = monoIn + feedbackInR * envFeedbackMod;

        diffuserReadPos3 = diffuserWritePos3 - (diffuserBufferLen3 + timeCvControl + lfo2b * 48);
        //diffuserReadPos3 = diffuserWritePos3 + 1 + lfo2b;

    	while( (diffuserReadPos3 < 0) )
    		diffuserReadPos3 += diffuserBufferLen3;
    	while( (diffuserReadPos3 >= diffuserBufferLen3M1) )
    		diffuserReadPos3 -= diffuserBufferLen3M1;

    	float int3 = diffuser3Interpolation(diffuserReadPos3);
    	//float int3 = diffuser3HermiteInterpolation((int) diffuserReadPos3);

        diffuserBuffer3[diffuserWritePos3] 		=  ap3In + diffuserCoef1 * int3;
        ap3Out = -ap3In * diffuserCoef1 + 	int3 * diffuserCoef1b;

    	// ----------------------------------------------< inject in delay3
    	delay3Buffer[ delay3WritePos ] 		= ap3Out;

    	// ----------------------------------------------> read delay3
		ap4In = delay3Interpolation(delay3ReadPos);

        // --- filter
        v2R += loopHpf * v3R;						// hipass
        v3R += loopHpf * (ap4In - v2R - v3R);
        //v4R += loopLpf2 * v5R;						// lowpass
        //v5R += loopLpf2 * ((ap4In-v2R) - v4R - v5R);

        ap4In = (ap4In-v2R) * envFeedbackMod;				// decay

        // ---- ap 4

        diffuserReadPos4 = diffuserWritePos4 + 1;
    	if( diffuserReadPos4 >= diffuserBufferLen4M1 )
    		diffuserReadPos4 -= diffuserBufferLen4M1;

        diffuserBuffer4[diffuserWritePos4] 		= ap4In - diffuserCoef2 * diffuserBuffer4[diffuserReadPos4];
        ap4Out = ap4In * diffuserCoef2 + 	diffuser4Interpolation(diffuserReadPos4) * diffuserCoef2b;

    	// ----------------------------------------------< inject in delay4
    	delay4Buffer[ delay4WritePos ] 		= ap4Out;

    	// ----------------------------------------------> read delay4
    	ap4Out = delay4Interpolation(delay4ReadPos);

    	// ================================================  mix out

    	//float attnA = 0.66f + lfo2 * 0.33f;
    	//float attnB = lfo2b * 0.33f;
    	feedbackInL = ap4Out;// * attnA + ap2In * attnB;
    	feedbackInR = ap2Out;// * attnA + ap4In * attnB;

    	outL = ap1In;// + ap2Out * 0.25f;
    	outR = ap3In;// + ap4Out * 0.25f;

    	*(outBuff++) += (int32_t) ( outL * sampleMultipler);
    	*(outBuff++) += (int32_t) ( outR * sampleMultipler);

    	// ================================================ index increment

    	sample += 2;

    	inputWritePos1		++;
    	inputWritePos2		++;
    	inputWritePos3		++;
    	inputWritePos4		++;
    	if( inputWritePos1 >= inputBufferLen1 )
    		inputWritePos1 = 0;
    	if( inputWritePos2 >= inputBufferLen2 )
    		inputWritePos2 = 0;
    	if( inputWritePos3 >= inputBufferLen3 )
    		inputWritePos3 = 0;
    	if( inputWritePos4 >= inputBufferLen4 )
    		inputWritePos4 = 0;

    	predelayWritePos ++;
    	if( predelayWritePos >= predelayBufferSize )
    		predelayWritePos = 0;

        diffuserWritePos1 	++;
        diffuserWritePos2 	++;
        diffuserWritePos3 	++;
        diffuserWritePos4 	++;
    	if( diffuserWritePos1 >= diffuserBufferLen1 )
    		diffuserWritePos1 = 0;
    	if( diffuserWritePos2 >= diffuserBufferLen2 )
    		diffuserWritePos2 = 0;
    	if( diffuserWritePos3 >= diffuserBufferLen3 )
    		diffuserWritePos3 = 0;
    	if( diffuserWritePos4 >= diffuserBufferLen4 )
    		diffuserWritePos4 = 0;

    	delay1WritePos 	++;
    	if( delay1WritePos >= delay1BufferSize )
    		delay1WritePos = 0;

    	delay2WritePos 	++;
    	if( delay2WritePos >= delay2BufferSize )
    		delay2WritePos = 0;

    	delay3WritePos 	++;
    	if( delay3WritePos >= delay3BufferSize )
    		delay3WritePos = 0;

    	delay4WritePos 	++;
    	if( delay4WritePos >= delay4BufferSize )
    		delay4WritePos = 0;

    	// --- lfo increment

    	lfo1tri -= lfo1Inc * fxSpeed;
    	if(lfo1tri >= 1) {
    		lfo1tri = 1;
    		lfo1Inc = -lfo1Inc;
    	}
    	if(lfo1tri <= -1) {
    		lfo1tri = -1;
    		lfo1Inc = -lfo1Inc;
    	}
    	lfo1 = sigmoid(lfo1tri);

    	lfo2IncMod = lfo2IncMod + lfo2IncModSampleInc;

    	lfo2tri += lfo2IncMod;
    	if(lfo2tri >= 1) {
    		lfo2tri = 1;
    		lfo2IncMod = -lfo2IncMod;
    	}
    	if(lfo2tri <= 0) {
    		lfo2tri = 0;
    		lfo2IncMod = -lfo2IncMod;
    	}
    	lfo2btri = lfo2tri + 0.5f;
    	if(lfo2btri >= 1) {
    		lfo2btri = -1;
    	}

    	lfo2 = tri2sin(lfo2tri);
    	lfo2b = tri2sin(lfo2btri);

    	lfoTremolo		+= lfoTremoloInc * fxTremoloSpeed * tremoloEnvFollowMod;
    	lfoTremoloSin	= tri2sin(lfoTremolo);

    	if(lfoTremolo >= 1 ) {
    		lfoTremolo = 1;
    		lfoTremoloInc = -lfoTremoloInc;
    	}
    	if(lfoTremolo <= 0) {
    		lfoTremolo = 0;
    		lfoTremoloInc = -lfoTremoloInc;
    	}

    	// bounce

    	prevTimeCv 	= timeCv;
    	timeCv 		= (lfo1 * lfoDepth + envMod * 0.1f) * delay3BufferSize;
    	cvDelta = -(prevTimeCv - timeCv) * 0.07f;
    	timeCvSpeed += cvDelta;
    	timeCvSpeed *= 0.9f;
    	bouncingCv 	+= timeCvSpeed;

    	timeCvControl = bounceLevel * bouncingCv + (1 - bounceLevel) * timeCv;
    }
}
float FxBus::predelayInterpolation(float readPos) {
	int readPosInt = (int) readPos;
	float y0 = predelayBuffer[readPosInt];
	float y1 = predelayBuffer[readPosInt + 1];

	float x = readPos - floorf(readPos);

	return y0 * (1 - x) + y1 * x;
}
float FxBus::diffuser1Interpolation(float readPos) {
	int readPosInt = (int) readPos;
	float y0 = diffuserBuffer1[readPosInt];
	float y1 = diffuserBuffer1[readPosInt + 1];

	float x = readPos - floorf(readPos);

	return y0 * (1 - x) + y1 * x;
}
float FxBus::diffuser1HermiteInterpolation(int readPos) {

	float y0 = diffuserBuffer1[(readPos - 1 + diffuserBufferLen1) % (diffuserBufferLen1) ];
	float y1 = diffuserBuffer1[readPos];
	float y2 = diffuserBuffer1[(readPos + 1)% (diffuserBufferLen1)];
	float y3 = diffuserBuffer1[(readPos + 2)% (diffuserBufferLen1)];

	float x = y2 - y1;
	x -= floorf(x);

    return hermite4(x, y0, y1, y2, y3);
}
float FxBus::diffuser2Interpolation(float readPos) {
	int readPosInt = (int) readPos;
	float y0 = diffuserBuffer2[readPosInt];
	float y1 = diffuserBuffer2[readPosInt + 1];

	float x = readPos - floorf(readPos);

	return y0 * (1 - x) + y1 * x;
}
float FxBus::diffuser3Interpolation(float readPos) {
	int readPosInt = (int) readPos;
	float y0 = diffuserBuffer3[readPosInt];
	float y1 = diffuserBuffer3[readPosInt + 1];

	float x = readPos - floorf(readPos);

	return y0 * (1 - x) + y1 * x;
}
float FxBus::diffuser3HermiteInterpolation(int readPos) {

	float y0 = diffuserBuffer3[(readPos - 1 + diffuserBufferLen3) % (diffuserBufferLen3) ];
	float y1 = diffuserBuffer3[readPos];
	float y2 = diffuserBuffer3[(readPos + 1)% (diffuserBufferLen3)];
	float y3 = diffuserBuffer3[(readPos + 2)% (diffuserBufferLen3)];

	float x = y2 - y1;
	x -= floorf(x);

    return hermite4(x, y0, y1, y2, y3);
}
float FxBus::diffuser4Interpolation(float readPos) {
	int readPosInt = (int) readPos;
	float y0 = diffuserBuffer4[readPosInt];
	float y1 = diffuserBuffer4[readPosInt + 1];

	float x = readPos - floorf(readPos);

	return y0 * (1 - x) + y1 * x;
}

float FxBus::delay1Interpolation(float readPos) {
	int readPosInt = (int) readPos;
	float y0 = delay1Buffer[readPosInt];
	float y1 = delay1Buffer[readPosInt + 1];

	float x = readPos - floorf(readPos);

    return y0 * (1 - x) + y1 * x;
}
float FxBus::delay2Interpolation(float readPos) {
	int readPosInt = (int) readPos;
	float y0 = delay2Buffer[readPosInt];
	float y1 = delay2Buffer[readPosInt + 1];

	float x = readPos - floorf(readPos);

    return y0 * (1 - x) + y1 * x;
}
float FxBus::delay3Interpolation(float readPos) {
	int readPosInt = (int) readPos;
	float y0 = delay3Buffer[readPosInt];
	float y1 = delay3Buffer[readPosInt + 1];

	float x = readPos - floorf(readPos);

    return y0 * (1 - x) + y1 * x;
}
float FxBus::delay4Interpolation(float readPos) {
	int readPosInt = (int) readPos;
	float y0 = delay4Buffer[readPosInt];
	float y1 = delay4Buffer[readPosInt + 1];

	float x = readPos - floorf(readPos);

    return y0 * (1 - x) + y1 * x;
}

float * crossFade(float t) {
	static float volumes[2];
	volumes[0] = sqrt3(0.5f * (1 + t));
	volumes[1] = sqrt3(0.5f * (1 - t));
	return volumes;
}


