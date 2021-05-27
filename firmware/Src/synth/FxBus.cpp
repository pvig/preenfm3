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
float FxBus::predelayBuffer[predelayBufferSize] __attribute__((section(".ram_d1")));
float FxBus::tapDelayBuffer[tapDelayBufferSize] __attribute__((section(".ram_d1")));
int FxBus::tapDelayPos[tapCount] __attribute__((section(".ram_d1")));
float FxBus::tapDelayAmp[tapCount] __attribute__((section(".ram_d1")));
float FxBus::diffuserBuffer1[diffuserBufferLen1] __attribute__((section(".ram_d1")));
float FxBus::diffuserBuffer2[diffuserBufferLen2] __attribute__((section(".ram_d1")));
float FxBus::diffuserBuffer3[diffuserBufferLen3] __attribute__((section(".ram_d1")));
float FxBus::diffuserBuffer4[diffuserBufferLen4] __attribute__((section(".ram_d1")));


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

	monoInHpf = 0.1443f;
	loopLpf = 0.348f;
	harmTremoloCutF = 0.38f;

	inLpF = 0.23f;
	loopHp = 0.1f;

	diffuserCoef1b =  (1  - (diffuserCoef1 * diffuserCoef1));
	diffuserCoef2b =  (1  - (diffuserCoef2 * diffuserCoef2));

	/*tapDelayPos[0] 	= 228;
	tapDelayPos[1] 	= 2046;
	tapDelayPos[2] 	= 1900;
	tapDelayPos[3] 	= 1500;
	tapDelayPos[4] 	= 722;
	tapDelayPos[5] 	= 322;*/

	/*tapDelayPos[0] 	= 228;
	tapDelayPos[1] 	= 2046;
	tapDelayPos[2] 	= 1500;
	tapDelayPos[3] 	= 1300;
	tapDelayPos[4] 	= 722;
	tapDelayPos[5] 	= 322;*/

	tapDelayPos[0] 	= 142;
	tapDelayPos[1] 	= 107;
	tapDelayPos[2] 	= 379;
	tapDelayPos[3] 	= 277;
	tapDelayPos[4] 	= 228;
	tapDelayPos[5] 	= 322;

	tapDelayAmp[0] 	= 0.5f;
	tapDelayAmp[1] 	= 0.5f;
	tapDelayAmp[2] 	= -0.5f;
	tapDelayAmp[3] 	= 0.5f;
	tapDelayAmp[4] 	= -0.5f;
	tapDelayAmp[5] 	= -0.5f;

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
	prevTime 		=  	clamp( fxTimeLinear, 0.0003f, 0.9997f);
	prevTime		*= 	prevTime * prevTime;
	fxTime 			= 	fxTime * 0.9f + prevTime * 0.1f;

	delay1FxTarget = 	getQuantizedTime(fxTime, delay1BufferSize);
	delay1DelayLen = 	(delay1DelayLen 	+ (delay1FxTarget - delay1DelayLen)	* 0.01f);// * (1 + (noise[0] - 0.5f) * 0.000035f);

	delay3FxTarget = 	getQuantizedTime(fxTime, delay3BufferSize);
	delay3DelayLen = 	(delay3DelayLen 	+ (delay3FxTarget - delay3DelayLen)	* 0.01f);// * (1 + (noise[0] - 0.5f) * 0.000035f);

	predelaySize 	= fxTime * predelayBufferSize;

	prevSpeed 		= 	synthState_->fullState.masterfxConfig[ GLOBALFX_LFOSPEED ];

	temp 			= 	prevSpeed;
	temp 			*=	temp * temp;
	fxSpeed 		= 	fxSpeed * 0.9f + temp * 0.1f;

    feedbackGain 		= 	synthState_->fullState.masterfxConfig[ GLOBALFX_FBACK ];
    fxInputLevel 	= 	synthState_->fullState.masterfxConfig[ GLOBALFX_INPUTLEVEL ];
    fxInputLevelAbs	=	1 - fabsf(fxInputLevel);

    feedback1 = feedbackGain;
    feedback2 = feedbackGain * 0.99f;


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

	//rand1 = lfo2  ;
	//rand2 = lfo2b ;
	//rand3 = 2 * (int)(fabsf(noise[2]) * 2);

	if (lfo2ChangeCounter++ > lfo2ChangePeriod) {
		lfo2ChangeCounter = 0;
		lfo2ModVal = noise[0] * noise[0] * 0.5f;
		lfo2IncMod = lfo2Inc * (1 + lfo2ModVal);
	}

	if (loopDecouplerChangeCounter++ > loopDecouplerChangePeriod) {
		loopDecouplerChangeCounter = 0;
		loopDecouplerModVal = noise[0] * noise[0];
		loopDecoupler = decoupler1 * loopDecouplerModVal;
		loopDecoupler2 = decoupler2 * loopDecouplerModVal;
	}

	float tapDelayMod = 0.6f + fxTime * 0.4f;
	tapDelayPos[0] 	= (int)(228 * tapDelayMod)&0xfffffffe;
	tapDelayPos[1] 	= (int)(2046 * tapDelayMod)&0xfffffffe;
	tapDelayPos[2] 	= (int)(1900 * tapDelayMod)&0xfffffffe;
	tapDelayPos[3] 	= (int)(1500 * tapDelayMod)&0xfffffffe;
	tapDelayPos[4] 	= (int)(722 * tapDelayMod)&0xfffffffe;
	tapDelayPos[5] 	= (int)(322 * tapDelayMod)&0xfffffffe;
}

/**
 * add timbre block to bus mix
 */
void FxBus::mixAdd(float *inStereo, int timbreNum) {
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
    const float sampleMultipler = 40.5f * (float) 0x7fffff; // fx level

	float tremoloEnvFollowMod, tremoloEnvFollowModAttn;
	float tremoloMod = 0, inLpMod;
	float tapAccumulator;
	float declouple1, declouple2, declouple3;
	float loopOutL, loopOutR;
	int tapPos;

	for (int s = 0; s < BLOCK_SIZE; s++) {

		declouple1 = lfo2 * loopDecoupler;
		declouple2 = lfo2 * loopDecoupler2;
		declouple3 = lfo2b * loopDecoupler2;

    	delay1ReadPos = delay1WritePos - (delay1DelayLen + timeCvControl);
    	while( delay1ReadPos < 0 )
    		delay1ReadPos += delay1BufferSize;
    	while( delay1ReadPos >= delay1BufferSizeM1 )
    		delay1ReadPos -= delay1BufferSizeM1;


    	delay2ReadPos = delay2WritePos - (delay2DelayLen * (1 + declouple1) );;
    	while( delay2ReadPos < 0 )
    		delay2ReadPos += delay2BufferSize;
    	while( delay2ReadPos >= delay2BufferSizeM1 )
    		delay2ReadPos -= delay2BufferSizeM1;

    	delay3ReadPos = delay3WritePos - (delay3DelayLen + lfo2);
    	while( delay3ReadPos < 0 )
    		delay3ReadPos += delay3BufferSize;
    	while( delay3ReadPos >= delay3BufferSizeM1 )
    		delay3ReadPos -= delay3BufferSizeM1;

    	delay4ReadPos = delay4WritePos - (delay4DelayLen + lfo2b );
    	while( delay4ReadPos < 0 )
    		delay4ReadPos += delay4BufferSize;
    	while( delay4ReadPos >= delay4BufferSizeM1 )
    		delay4ReadPos -= delay4BufferSizeM1;

        // --- audio in

    	inR = *(sample);
    	inL = *(sample + 1);

        // --- cut high

    	inLpMod = inLpF + envelope * 0.05f;

        v6L += inLpMod * v7L;
        v7L += inLpMod * (inL - v6L - v7L);
        v6L += inLpMod * v7L;
        v7L += inLpMod * (inL - v6L - v7L);

        v6R += inLpMod * v7R;
        v7R += inLpMod * (inR - v6R - v7R);
        v6R += inLpMod * v7R;
        v7R += inLpMod * (inR - v6R - v7R);

        nodeL = v6L;
        nodeR = v6R;

    	// --- audio in > harmonic tremolo

        tremoloEnvFollowMod 	= 1 + envelope * tremoloEnvFollow;
        tremoloEnvFollowModAttn = (tremoloEnvFollowMod + (1 - tremoloEnvFollowAbs));

        _lx3L += harmTremoloCutF * _ly3L; // low pass filter
        _ly3L += harmTremoloCutF * ((v6L) - _lx3L - _ly3L);

        lpL = _lx3L;
        hpL = inL - lpL;

        /*_lx3R += harmTremoloCutF * _ly3R;
        _ly3R += harmTremoloCutF * (v6R - _lx3R - _ly3R);

        lpR = _lx3R;
        hpR = inR - lpR;*/

    	//tremoloModL = ((lfoTremoloSin 		* fxTremoloDepth + 1 - fxTremoloDepth) * 2 - 1) * tremoloEnvFollowModAttn;
    	//tremoloModR = (((1 - lfoTremoloSin) * fxTremoloDepth + 1 - fxTremoloDepth) * 2 - 1) * tremoloEnvFollowModAttn;

        //combInL = lpL * tremoloModL + hpL * (1 - tremoloModL);
        //combInR = lpR * tremoloModR + hpR * (1 - tremoloModR);

    	tremoloMod = ((lfoTremoloSin 		* fxTremoloDepth + 1 - fxTremoloDepth) * 2 - 1) * tremoloEnvFollowModAttn;
        tremoloOut = lpL * tremoloMod + hpL * (1 - tremoloMod);

    	// --- enveloppe calculation

        blocksum 	+= fabsf(tremoloOut);
        envelope 	= (envelope * envM1 + envDest) * envM2;
        envMod 		= envModDepth * envelope;

    	// initial tap

        tapDelayBuffer[ tapDelayWritePos++ ]  = tremoloOut;

        if(tapDelayWritePos >= tapDelayBufferSize) {
        	tapDelayWritePos = 0;
        }

        tapAccumulator = 0;

        for(int aa=0; aa<tapCount; aa++) {
            tapPos = (tapDelayWritePos + tapDelayPos[aa])&(tapDelayBufferSizeM1);
            tapAccumulator += tapDelayBuffer[ tapPos ] * tapDelayAmp[aa];
        }

        monoIn = (tapAccumulator) * (tremoloMod * 2 - 1) * 0.5f;

        // --- hi pass

        v0R += monoInHpf * v1R;						// hipass
        v1R += monoInHpf * ( monoIn - v0R - v1R);
        monoIn -= v0L;

    	//--- pre delay

    	/*predelayBuffer[predelayWritePos] = monoIn;
    	predelayReadPos = predelayWritePos - predelaySize;
    	monoIn = predelayBuffer[predelayReadPos];

    	while( predelayReadPos < 0 )
    		predelayReadPos += predelayBufferSize;

    	predelayWritePos += 2;
    	if( predelayWritePos >= predelayBufferSize )
    		predelayWritePos = 0;*/

    	// --- input tremolo

        // diffuser

        // ---- ap 1

        ap1In = monoIn + clamp( ap4Out * feedback1 * (1 + envMod + declouple1), -1, 1);

        diffuserReadPos1 = diffuserWritePos1 + 1 + declouple3;
    	if( diffuserReadPos1 >= diffuserBufferLen1M1 )
    		diffuserReadPos1 -= diffuserBufferLen1M1;

    	float int1 = diffuser1Interpolation(diffuserReadPos1);
        diffuserBuffer1[diffuserWritePos1] 		= ap1In - diffuserCoef1 * int1;
        ap1Out = ap1In * diffuserCoef1 + 	int1 * diffuserCoef1b;

    	// ----------------------------------------------< inject in delay1
    	delay1Buffer[ delay1WritePos ] 		= ap1Out;


    	// ----------------------------------------------> read delay1
		ap2In = delay1Interpolation(delay1ReadPos);

        // --- low pass

        v2L += loopLpf * v3L;						// lowpass
        v3L += loopLpf * ( ap2In - v2L - v3L);
        v2L += loopLpf * v3L;						// lowpass
        v3L += loopLpf * ( ap2In - v2L - v3L);
        v2L += loopLpf * v3L;						// lowpass
        v3L += loopLpf * ( ap2In - v2L - v3L);

        ap2In = v2L;

    	// ---- ap 2

        diffuserReadPos2 = diffuserWritePos2 + 1;
    	if( diffuserReadPos2 >= diffuserBufferLen2M1 )
    		diffuserReadPos2 -= diffuserBufferLen2M1;
        diffuserBuffer2[diffuserWritePos2] 		= ap2In - diffuserCoef1 * diffuserBuffer2[diffuserReadPos2];

        ap2Out = ap2In * diffuserCoef1 + 	diffuserBuffer2[diffuserReadPos2] * diffuserCoef1b;

    	// ----------------------------------------------< inject in delay2
    	delay2Buffer[ delay2WritePos ] 		= ap2Out;

    	// ----------------------------------------------> read delay2
    	ap2Out = delay2Interpolation(delay2ReadPos);

        // ---- ap 3

		ap3In = monoIn + clamp( ap2Out * feedback2 * (1 + envMod + lfo2b * loopDecoupler), -1, 1);

        diffuserReadPos3 = diffuserWritePos3 + 1 + declouple2;
    	if( diffuserReadPos3 >= diffuserBufferLen3M1 )
    		diffuserReadPos3 -= diffuserBufferLen3M1;

    	float int3 = diffuser3Interpolation(diffuserReadPos3);
        diffuserBuffer3[diffuserWritePos3] 		=  ap3In - diffuserCoef2 * int3;
        ap3Out = ap3In * diffuserCoef2 + 	int3 * diffuserCoef2b;

    	// ----------------------------------------------< inject in delay3
    	delay3Buffer[ delay3WritePos ] 		= ap3Out;

    	// ----------------------------------------------> read delay3
		ap4In = delay3Interpolation(delay3ReadPos);

        // --- low pass
        v2R += loopLpf * v3R;						// lowpass
        v3R += loopLpf * (ap4In - v2R - v3R);
        v2R += loopLpf * v3R;						// lowpass
        v3R += loopLpf * (ap4In - v2R - v3R);
        v2R += loopLpf * v3R;						// lowpass
        v3R += loopLpf * (ap4In - v2R - v3R);

        ap4In = v2R;

        // ---- ap 4

        diffuserReadPos4 = diffuserWritePos4 + 1;
    	if( diffuserReadPos4 >= diffuserBufferLen4M1 )
    		diffuserReadPos4 -= diffuserBufferLen4M1;
        diffuserBuffer4[diffuserWritePos4] 		= ap4In - diffuserCoef2 * diffuserBuffer4[diffuserReadPos4];

        ap4Out = ap4In * diffuserCoef2 + 	diffuserBuffer4[diffuserReadPos4] * diffuserCoef2b;

    	// ----------------------------------------------< inject in delay4
    	delay4Buffer[ delay4WritePos ] 		= ap4Out;

    	// ----------------------------------------------> read delay4
    	ap4Out = delay4Interpolation(delay4ReadPos);

        // increment

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


    	// ================================================  mix out

    	loopOutL = ap1Out - ap2In + ap2Out;
    	loopOutR = ap3Out - ap4In + ap4Out;

    	outL = fxInputLevel * nodeL + loopOutL;
        outR = fxInputLevel * nodeR + loopOutR;

    	*(outBuff++) += (int32_t) ( outL * sampleMultipler);
    	*(outBuff++) += (int32_t) ( outR * sampleMultipler);

    	sample += 2;

    	// --- write index increment

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
    	timeCv 		= (lfo1 * lfoDepth + envMod * 0.1f) * delay1FxTarget;
    	cvDelta = -(prevTimeCv - timeCv) * 0.07f;
    	timeCvSpeed += cvDelta;
    	timeCvSpeed *= 0.9f;
    	bouncingCv 	+= timeCvSpeed;

    	timeCvControl = bounceLevel * bouncingCv + (1 - bounceLevel) * timeCv;
    }
}

float FxBus::delay1HermiteInterpolation(int readPos) {

	float y0 = delay1Buffer[(readPos + delay1BufferSize - 1) & (delay1BufferSize - 1) ];
	float y1 = delay1Buffer[readPos];
	float y2 = delay1Buffer[(readPos + 1)& (delay1BufferSize - 1)];
	float y3 = delay1Buffer[(readPos + 2)& (delay1BufferSize - 1)];

	float x = y2 - y1;
	x -= floorf(x);

    return hermite4(x, y0, y1, y2, y3);
}
float FxBus::diffuser1Interpolation(float readPos) {
	int readPosInt = (int) readPos;
	float y0 = diffuserBuffer1[readPosInt];
	float y1 = diffuserBuffer1[readPosInt + 1];

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



