/*
 * FxBus.cpp
 *
 *  Created on: Feb 14, 2021
 *      Author: patvig
 */

#include "FxBus.h"

#define excursion 18

extern float diatonicScaleFrequency[];
extern float sinTable[];

inline float fold(float x4) {
    // https://www.desmos.com/calculator/ge2wvg2wgj
    // x4 = x / 4
    return (fabsf(x4 + 0.25f - roundf(x4 + 0.25f)) - 0.25f);
}
inline
float expf_fast(float a) {
  //https://github.com/ekmett/approximate/blob/master/cbits/fast.c
  union { float f; int x; } u;
  u.x = (int) (12102203 * a + 1064866805);
  return u.f;
}
inline
int modulo(int d, int max) {
  return unlikely(d >= max) ? d - max : d;
}
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


extern float noise[32];

FxBus::FxBus() {
	lfo1 = 0;
}

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

	inHpf = 0.1f;
	loopHpf = 0.1f;
	loopLpf2 = 0.95f;

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

    feedbackGain 		= 	synthState_->fullState.masterfxConfig[ GLOBALFX_FBACK ] * 0.91f;
    envRelease 			= 	0.005f + feedbackGain * feedbackGain * 0.7f;

    predelayMixLevel 	= 	synthState_->fullState.masterfxConfig[ GLOBALFX_PREDELAYMIX ];
    predelayMixAttn 	= 	predelayMixLevel * (1 - (predelayMixLevel * predelayMixLevel * 0.1f));

	prevBounce 		= 	synthState_->fullState.masterfxConfig[ GLOBALFX_BOUNCE ] + 0.005f ;
	bounceLevel	 	= 	bounceLevel * 0.9f + prevBounce * 0.1f;

	/*prevEnvRelease 	= 	synthState_->fullState.masterfxConfig[ GLOBALFX_ENVRELEASE ] + 0.005f ;
	prevEnvRelease	*= 	prevEnvRelease * prevEnvRelease;
	envRelease	 	= 	envRelease * 0.9f + prevEnvRelease * 0.1f;*/

	inputDiffusion 	= synthState_->fullState.masterfxConfig[GLOBALFX_INPUTDIFFUSION];
	if(inputDiffusion != prevInputDiffusion) {
		inputCoef1 = 0.1f + inputDiffusion * 0.8f;
		inputCoef2 = 0.1f + inputDiffusion * 0.712f;
		inputCoef1b =  (1  - (inputCoef1 * inputCoef1));
		inputCoef2b =  (1  - (inputCoef2 * inputCoef2));
	}
	prevInputDiffusion = inputDiffusion;

	decayDiffusion 	= synthState_->fullState.masterfxConfig[GLOBALFX_DECAYDIFFUSION];
	if(decayDiffusion != prevDecayDiffusion) {
		diffuserCoef1 = 0.1f + decayDiffusion * 0.8f;
		diffuserCoef2 = 0.1f + decayDiffusion * 0.6f;
		diffuserCoef1b =  (1  - (diffuserCoef1 * diffuserCoef1));
		diffuserCoef2b =  (1  - (diffuserCoef2 * diffuserCoef2));
	}
	prevDecayDiffusion = decayDiffusion;

	damping 		= synthState_->fullState.masterfxConfig[GLOBALFX_INPUTDAMPING] * 0.6f;
	damping 		*= damping;
	loopLpf 		= 0.05f + damping * 0.95f;

    // ------ page 2

	speedLinear 	= 	synthState_->fullState.masterfxConfig[ GLOBALFX_LFOSPEED ];
	temp 			= 	speedLinear;
	temp 			*=	temp * temp;
	fxSpeed 		= 	fxSpeed * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_LFODEPTH ] * 0.7f;
    temp 			*= 	1 - speedLinear * 0.9f;
    lfoDepth 		= 	lfoDepth * 0.9f + temp * 0.1f;

    temp = 	synthState_->fullState.masterfxConfig[ GLOBALFX_ENVTHRESHOLD];
	prevEnvThreshold = temp * temp * headRoomDivider * 640;
	envThreshold	= 	envThreshold * 0.9f + prevEnvThreshold * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_ENVMOD] * 0.85f;
    envModDepth 	= 	envModDepth * 0.9f + temp * 0.1f;
	envModDepthNeg	=	fabsf(envModDepth);

    temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_INPUTTILT ] * 0.75f;
    temp			*= 	temp;
    tiltInput		= 	tiltInput * 0.9f + temp * 0.1f;

    inHpf			= 	clamp(tiltInput * fold(tiltInput * 20) , 0, 1);
    inLpF			=	clamp(tiltInput + 0.02f, 0, 1);

    temp = 	synthState_->fullState.masterfxConfig[ GLOBALFX_ENVFEEDBACK] * 0.25f;
	envFeedback	= 	envFeedback * 0.9f + temp * 0.1f;


    // ------ env follow

	envBlocknn += BLOCK_SIZE;

	if(envBlocknn > envDetectSize) {
		if( blocksum > envThreshold) {
			// attack
			envDest = 1;
			envM1 = 999;
			envM2 = 0.001f;
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

	/*if (lfo2ChangeCounter++ > lfo2ChangePeriod) {
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
	}*/

}

/**
 * add timbre block to bus mix
 */
void FxBus::mixAdd(float *inStereo, int timbreNum) {
	const float level = sqrt3(synthState_->mixerState.instrumentState_[timbreNum].send) * headRoomDivider; // divide for more headroom

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
    const float sampleMultipler = headRoomMultiplier * (float) 0x7fffff;

	float inLpMod;

	for (int s = 0; s < BLOCK_SIZE; s++) {

        // --- audio in

    	inR = *(sample);
    	inL = *(sample + 1);

    	monoIn = (inR + inL);

        // --- cut high

    	inLpMod = inLpF + envelope * 0.05f;

        v6R += inLpMod * v7R;
        v7R += inLpMod * (monoIn - v6R - v7R);
        v6L += inLpMod * v7L;
        v7L += inLpMod * (v6R - v6L - v7L);

        monoIn = v6L;

        // --- hi pass

        v0R += inHpf * v1R;						// hipass
        v1R += inHpf * ( monoIn - v0R - v1R);
        v0R += inHpf * v1R;						// hipass
        v1R += inHpf * ( monoIn - v0R - v1R);

        monoIn -= v0R;

    	// --- enveloppe calculation

        blocksum 	+= 	fabsf(monoIn);
        envelope 	= 	(envelope * envM1 + envDest) * envM2;
      	float envprep = (envDest == 1)? envelope*envelope : envelope;
        envMod 		= 	(envMod * 99 + (envModDepth < 0 ? envModDepthNeg * (1 - envprep) : envModDepth * envprep)) * 0.01f;

    	//--- pre delay

    	predelayReadPos = predelayWritePos - predelaySize;
    	if( unlikely(predelayReadPos < 0) )
    		predelayReadPos += predelayBufferSize;
    	if( unlikely(predelayReadPos >= predelayBufferSizeM1) )
    		predelayReadPos -= predelayBufferSizeM1;

    	predelayBuffer[predelayWritePos] = monoIn;

    	monoIn = predelayMixAttn * delayInterpolation(predelayReadPos, predelayBuffer, predelayBufferSizeM1) + (1 - predelayMixAttn) * monoIn;

    	// --- input diffuser

        // ---- diffuser 1

        inputReadPos1 = modulo(inputWritePos1 + 1, inputBufferLen1);
    	inputBuffer1[inputWritePos1] = monoIn - inputCoef1 * inputBuffer1[inputReadPos1];
    	diff1Out = monoIn * inputCoef1 + 	inputBuffer1[inputReadPos1] * inputCoef1b;

        // ---- diffuser 2

    	inputReadPos2 = modulo(inputWritePos2 + 1, inputBufferLen2);
    	inputBuffer2[inputWritePos2] = diff1Out - inputCoef1 * inputBuffer2[inputReadPos2];
    	diff2Out = diff1Out * inputCoef1 + 	inputBuffer2[inputReadPos2] * inputCoef1b;

        // ---- diffuser 3

    	inputReadPos3 = modulo(inputWritePos3 + 1, inputBufferLen3);
    	inputBuffer3[inputWritePos3] = diff2Out - inputCoef2 * inputBuffer3[inputReadPos3];
    	diff3Out = diff2Out * inputCoef2 + 	inputBuffer3[inputReadPos3] * inputCoef2b;

        // ---- diffuser 4

    	inputReadPos4 = modulo(inputWritePos4 + 1, inputBufferLen4);
    	inputBuffer4[inputWritePos4] = diff3Out - inputCoef2 * inputBuffer4[inputReadPos4];
    	diff4Out = diff3Out * inputCoef2 + 	inputBuffer4[inputReadPos4] * inputCoef2b;

    	monoIn = diff4Out;

        v4R = monoIn - v5R + 0.999f * v4R;			// dc blocker
        v5R = monoIn;
        monoIn = v4R;

    	float decayVal = clamp((feedbackGain + envFeedback * envelope), 0, 1); // ------ decay

        // ---- ap 1

        ap1In = monoIn + feedbackInL * decayVal * (0.95f + lfo2 * 0.05f);

        diffuserReadPos1 = diffuserWritePos1 + 1 - (timeCvControl1);
        while( (diffuserReadPos1 >= diffuserBufferLen1) )
    		diffuserReadPos1 -= diffuserBufferLen1;
        while( (diffuserReadPos1 < 0) )
    		diffuserReadPos1 += diffuserBufferLen1;

    	float int1 = delayInterpolation(diffuserReadPos1, diffuserBuffer1, diffuserBufferLen1M1);
    	//float int1 = diffuserBuffer1[(int) diffuserReadPos1];

        diffuserBuffer1[diffuserWritePos1] 		= ap1In + diffuserCoef1 * int1;
        ap1Out 	= 	-ap1In * diffuserCoef1 	+ 	int1 * diffuserCoef1b;

    	// ----------------------------------------------< inject in delay1
    	delay1Buffer[ delay1WritePos ] 		= ap1Out;
    	// ----------------------------------------------> read delay1
		ap2In = delay1Buffer[ (int) delay1ReadPos ];

        // ---------------------------------------------------- filter
        v4L += loopLpf * v5L;						// lowpass
        v5L += loopLpf * ( ap2In - v4L - v5L);
        ap2In = v4L;

        ap2In *= decayVal;				// decay

    	// ---- ap 2

        diffuserReadPos2 = diffuserWritePos2 + 1 - (timeCvControl2);
    	while( (diffuserReadPos2 >= diffuserBufferLen2) )
    		diffuserReadPos2 -= diffuserBufferLen2;
        while( (diffuserReadPos2 < 0) )
    		diffuserReadPos2 += diffuserBufferLen2;

    	//float int2 = diffuser2Interpolation(diffuserReadPos2);
    	float int2 = delayInterpolation(diffuserReadPos2, diffuserBuffer2, diffuserBufferLen2M1);

        diffuserBuffer2[diffuserWritePos2] 		= ap2In - diffuserCoef2 * int2;
        ap2Out = ap2In * diffuserCoef2 + 	int2 * diffuserCoef2b;

    	// ----------------------------------------------< inject in delay2
    	delay2Buffer[ delay2WritePos ] 		= ap2Out;
    	// ----------------------------------------------> read delay2

    	ap2Out = delay2Buffer[ (int) delay2ReadPos ];

        // ---- ap 3

		ap3In = monoIn + feedbackInR * decayVal * (0.95f + lfo2b * 0.05f);

        diffuserReadPos3 = diffuserWritePos3 + 1 - (timeCvControl3);

        while( (diffuserReadPos3 >= diffuserBufferLen3 ) )
    		diffuserReadPos3 -= diffuserBufferLen3;
        while( (diffuserReadPos3 < 0) )
    		diffuserReadPos3 += diffuserBufferLen3;

    	//float int3 = diffuser3Interpolation(diffuserReadPos3);
    	float int3 = delayInterpolation(diffuserReadPos3, diffuserBuffer3, diffuserBufferLen3M1);
    	//float int3 = diffuserBuffer3[(int) diffuserReadPos3];

        diffuserBuffer3[diffuserWritePos3] 		=  ap3In + diffuserCoef1 * int3;
        ap3Out = -ap3In * diffuserCoef1 + 	int3 * diffuserCoef1b;

    	// ----------------------------------------------< inject in delay3
    	delay3Buffer[ delay3WritePos ] 		= ap3Out;
    	// ----------------------------------------------> read delay3
		ap4In = delay3Buffer[ (int) delay3ReadPos ];

        // ---------------------------------------------------- filter

        v2R = ap4In - v3R + 0.999f * v2R;			// dc blocker
        v3R = ap4In;
        ap4In = v2R;

		ap4In *= decayVal;				// decay

        // ---- ap 4

        diffuserReadPos4 = diffuserWritePos4 + 1 - (timeCvControl4);

        while( (diffuserReadPos4 >= diffuserBufferLen4 ) )
        	diffuserReadPos4 -= diffuserBufferLen4;
        while( (diffuserReadPos4 < 0) )
        	diffuserReadPos4 += diffuserBufferLen4;

    	//float int4 = diffuser4Interpolation(diffuserReadPos4);
    	float int4 = delayInterpolation(diffuserReadPos4, diffuserBuffer4, diffuserBufferLen4M1);

        diffuserBuffer4[diffuserWritePos4] 		= ap4In - diffuserCoef2 * int4;
        ap4Out = ap4In * diffuserCoef2 + 	int4 * diffuserCoef2b;

    	// ----------------------------------------------< inject in delay4
    	delay4Buffer[ delay4WritePos ] 		= ap4Out;
    	// ----------------------------------------------> read delay4
    	ap4Out = delay4Buffer[ (int) delay4ReadPos ];

    	// ================================================  mix out

    	feedbackInL = ap4Out;
    	feedbackInR = ap2Out;

    	//outL = ap1Out;// + ap1Out - ap2In - ap2Out ;
        outL = ap1Out;
        outL += delay1Buffer[ 		modulo(_kLeftTaps[0] + delay1WritePos, delay1BufferSize	)		];
        outL += delay1Buffer[ 		modulo(_kLeftTaps[1] + delay1WritePos, delay1BufferSize	)		];
        outL -= diffuserBuffer2[ 	modulo(_kLeftTaps[2] + diffuserWritePos2, diffuserBufferLen2)	];
        outL += delay2Buffer[ 		modulo(_kLeftTaps[3] + delay2WritePos, delay2BufferSize)		];
        outL -= delay3Buffer[ 		modulo(_kLeftTaps[4] + delay3WritePos, delay3BufferSize)		];
        outL -= diffuserBuffer4[ 	modulo(_kLeftTaps[5] + diffuserWritePos4, diffuserBufferLen4)	];
        outL -= delay4Buffer[ 		modulo(_kLeftTaps[6] + delay4WritePos, delay4BufferSize)		];

    	//outR = ap3Out;// + ap3Out - ap4In - ap4Out ;
        outR = ap3Out;
        outR += delay3Buffer[ 		modulo(_kRightTaps[0] + delay3WritePos, delay3BufferSize) 		];
        outR += delay3Buffer[ 		modulo(_kRightTaps[1] + delay3WritePos, delay3BufferSize) 		];
        outR -= diffuserBuffer4[ 	modulo(_kRightTaps[2] + diffuserWritePos4, diffuserBufferLen4) 	];
        outR += delay4Buffer[ 		modulo(_kRightTaps[3] + delay4WritePos, delay4BufferSize)	 	];
        outR -= delay1Buffer[ 		modulo(_kRightTaps[4] + delay1WritePos, delay1BufferSize)	 	];
        outR -= diffuserBuffer2[ 	modulo(_kRightTaps[5] + diffuserWritePos2, diffuserBufferLen2) 	];
        outR -= delay2Buffer[ 		modulo(_kRightTaps[6] + delay2WritePos, delay2BufferSize)	 	];

    	*(outBuff++) += (int32_t) ( outL * sampleMultipler);
    	*(outBuff++) += (int32_t) ( outR * sampleMultipler);

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

    	delay1ReadPos = modulo(delay1WritePos + 1, delay1BufferSize);
    	delay2ReadPos = modulo(delay2WritePos + 1, delay2BufferSize);
    	delay3ReadPos = modulo(delay3WritePos + 1, delay3BufferSize);
    	delay4ReadPos = modulo(delay4WritePos + 1, delay4BufferSize);

    	// --- lfo increment

    	// ------- lfo 1

    	lfo1tri -= lfo1Inc * fxSpeed;
    	if(lfo1tri >= 1) {
    		lfo1tri = 1;
    		lfo1Inc = -lfo1Inc;
    	}
    	if(lfo1tri <= 0) {
    		lfo1tri = 0;
    		lfo1Inc = -lfo1Inc;
    	}
    	lfo1 = tri2sin(lfo1tri);

    	lfo1btri = lfo1tri + 0.5f;
    	if(lfo1btri >= 1) {
    		lfo1btri -= 1;
    	}

    	lfo1b = tri2sin(lfo1btri);

    	// ------- lfo 2

    	lfo2tri += lfo2Inc * fxSpeed;
    	if(lfo2tri >= 1) {
    		lfo2tri = 1;
    		lfo2Inc = -lfo2Inc;
    	}
    	if(lfo2tri <= 0) {
    		lfo2tri = 0;
    		lfo2Inc = -lfo2Inc;
    	}
    	lfo2btri = lfo2tri + 0.5f;
    	if(lfo2btri >= 1) {
    		lfo2btri -= 1;
    	}

    	lfo2 = tri2sin(lfo2tri);
    	lfo2b = tri2sin(lfo2btri);

    	// -------- mods :

    	float envModdiv = envMod;// * envMod;//* 0.5f;
    	timeCvControl1 = (lfo1  * lfoDepth) * diffuserBufferLen1;
    	timeCvControl2 = (lfo2b  * lfoDepth + envModdiv) * diffuserBufferLen2;
    	timeCvControl3 = (lfo2 * lfoDepth) * diffuserBufferLen3;
    	timeCvControl4 = (lfo1b * lfoDepth + envModdiv) * diffuserBufferLen4;
    }
}
float FxBus::predelayInterpolation(float readPos) {
	int readPosInt = (int) readPos;
	float y0 = predelayBuffer[readPosInt];
	float y1 = predelayBuffer[readPosInt + 1];
	float x = readPos - floorf(readPos);
    return y0 + x * (y1 - y0);
}
float FxBus::delayInterpolation(float readPos, float buffer[], int bufferLenM1) {
	int readPosInt = (int) readPos;
	float y0 = buffer[readPosInt];
	int next = readPosInt - 1;
	if(next<0)
		next = bufferLenM1;
	float y1 = buffer[next];
	float x = readPos - floorf(readPos);
    return y0 + x * (y1 - y0);
}
