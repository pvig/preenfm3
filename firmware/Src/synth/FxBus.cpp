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
float fold(float x4) {
    // https://www.desmos.com/calculator/ge2wvg2wgj
    // x4 = x / 4
    return (fabsf(x4 + 0.25f - roundf(x4 + 0.25f)) - 0.25f);
}
inline
float tri2sin(float x) {
	return x * x * (3 - x - x);
}
//***------------***------------***------------***------------***----- FxBus -------***------------***------------***------------

float FxBus::feedbackBuffer[FxBus::feedbackBufferSize] __attribute__((section(".ram_d1")));

FxBus::FxBus() {
	lfo1 = 0;
}

void FxBus::init(SynthState *synthState) {
    this->synthState_ = synthState;

    for (int s = 0; s < FxBus::feedbackBufferSize; s++) {
    	feedbackBuffer[ s ] = 0;
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

	vcfFreq = 0.88f;
	vcfDiffusion = 0.175f;

	feedbackLp = 0.095f;
	fxCrossover = 0.05f;
	harmTremoloCutF = 0.424f;

	inLpF = 0.56f;
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

	temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_TIME ];
	prevTime 		=  	clamp( temp, 0.0003f, 0.9997f);
	prevTime		*= 	prevTime * prevTime;
	fxTime 			= 	fxTime * 0.9f + prevTime * 0.1f;

	feedbackFxTarget = 	getQuantizedTime(fxTime, feedbackBufferSize);
	feedbackDelayLen = 	feedbackDelayLen 	+ (feedbackFxTarget - feedbackDelayLen)	* 0.01f;

	prevSpeed 		= 	synthState_->fullState.masterfxConfig[ GLOBALFX_LFOSPEED ];

	temp 			= 	prevSpeed;
	temp 			*=	temp * temp;
	fxSpeed 		= 	fxSpeed * 0.9f + temp * 0.1f;

    fxFeedback 		= 	synthState_->fullState.masterfxConfig[ GLOBALFX_FBACK ];
    fxInputLevel 	= 	synthState_->fullState.masterfxConfig[ GLOBALFX_INPUTLEVEL ];
    fxInputLevelAbs	=	1 - fabsf(fxInputLevel);

	prevEnvThreshold = 	synthState_->fullState.masterfxConfig[ GLOBALFX_ENVTHRESHOLD ];
	envThreshold	= 	envThreshold * 0.9f + prevEnvThreshold * 0.1f;

	prevBounce 		= 	synthState_->fullState.masterfxConfig[ GLOBALFX_BOUNCE ] + 0.005f ;
	bounceLevel	 	= 	bounceLevel * 0.9f + prevBounce * 0.1f;

	prevEnvRelease 	= 	synthState_->fullState.masterfxConfig[ GLOBALFX_ENVRELEASE ] + 0.005f ;
	prevEnvRelease	*= 	prevEnvRelease * prevEnvRelease;
	envRelease	 	= 	envRelease * 0.9f + prevEnvRelease * 0.1f;

    // ------ page 2

    temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_LFODEPTH ];
	invspeed 		= 	1 - sqrt3(fxTime * temp);
    temp 			*= 	temp * temp * (0.05f + invspeed * 0.95f);
    lfoDepth 		= 	lfoDepth * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_ENVMOD] * 0.99f;
    temp 			*= 	temp * (0.4f + invspeed * 0.6f);
    envModDepth 	= 	envModDepth * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_TREMOLOSPEED];
	fxTremoloSpeed	= 	fxTremoloSpeed * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_TREMOLODEPTH ];
	invspeed 		= 	1 - sqrt3(fxTremoloSpeed * temp);
    temp 			= 	temp * (0.7f + invspeed * 0.3f) ;
	fxTremoloDepth	= 	fxTremoloDepth * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ GLOBALFX_TREMOLOENVFOLLOW] * 0.99f;
    temp 			= 	temp * (0.8f + invspeed * 0.2f);
    tremoloEnvFollow= 	tremoloEnvFollow * 0.9f + temp * 0.1f;
    tremoloEnvFollowAbs = fabsf(tremoloEnvFollow);

    // ------

    fxTone 			= 	fxTime * 0.15f + 0.11f - lfo2 * 0.02f;

    // ------ env follow

	envBlocknn += BLOCK_SIZE;

	if(envBlocknn > envDetectSize) {
		if( blocksum > envThreshold) {
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

    // ----------- allpass params -----------

	const float bipolarf = (vcfFreq - 0.5f);
	const float folded = fold(sigmoid(bipolarf * 19 * vcfDiffusion)) * 4; // -1 < folded < 1
	const float lrDelta = 0.0095f * folded;
	const float range = 0.47f + lfo2 * 0.1f;

    f4L = clamp(((vcfFreq + lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
    coef4L = (1.0f - f4L) / (1.0f + f4L);

    f4R = clamp(((vcfFreq - lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
	coef4R = (1.0f - f4R) / (1.0f + f4R);

    // -----------
}

/**
 * add timbre block to bus mix
 */
void FxBus::mixAdd(float *inStereo, int timbreNum) {
	const float level = synthState_->mixerState.instrumentState_[timbreNum].send * 0.02f; // divide for more headroom

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
    const float sampleMultipler = 20.5f * (float) 0x7fffff; // fx level

	float feedbackAttn 	= 0.898f;
	float tremoloEnvFollowMod, tremoloEnvFollowModAttn;
	float tremoloModL, tremoloModR;

	for (int s = 0; s < BLOCK_SIZE; s++) {

    	feedbackReadPosL = feedbackWritePos - fabsf(feedbackDelayLen + timeCvControl);

    	while( feedbackReadPosL < 0 )
    		feedbackReadPosL += feedbackBufferSize;
    	while( feedbackReadPosL >= feedbackBufferSize )
    		feedbackReadPosL -= feedbackBufferSize;

    	feedbackReadPosIntL = (int) feedbackReadPosL;
    	feedbackReadPosIntL &= 0xfffffffe; // keep it even for stereo separation

    	feedbackReadPosR = feedbackWritePos - fabsf(feedbackDelayLen + timeCvControl + 75);

    	while( feedbackReadPosR < 0 )
    		feedbackReadPosR += feedbackBufferSize;
    	while( feedbackReadPosR >= feedbackBufferSize )
    		feedbackReadPosR -= feedbackBufferSize;

    	feedbackReadPosIntR = ((int) feedbackReadPosR)&0xfffffffe;
    	feedbackReadPosIntR += 1;

        // --- audio in

    	inR = *(sample);
    	inL = *(sample + 1);

        // --- cut high

        v6L += inLpF * v7L;
        v7L += inLpF * (inR - v6L - v7L);
        v6L += inLpF * v7L;
        v7L += inLpF * (inR - v6L - v7L);

        v6R += inLpF * v7R;
        v7R += inLpF * (inL - v6R - v7R);
        v6R += inLpF * v7R;
        v7R += inLpF * (inL - v6R - v7R);

    	// --- audio in > harmonic tremolo

        tremoloEnvFollowMod 	= 1 + envelope * tremoloEnvFollow;
        tremoloEnvFollowModAttn = (tremoloEnvFollowMod + (1 - tremoloEnvFollowAbs));

        _lx3L += harmTremoloCutF * _ly3L; // low pass filter
        _ly3L += harmTremoloCutF * ((v6L) - _lx3L - _ly3L);

        lpL = _lx3L;
        hpL = inL - lpL;

        _lx3R += harmTremoloCutF * _ly3R;
        _ly3R += harmTremoloCutF * (v6R - _lx3R - _ly3R);

        lpR = _lx3R;
        hpR = inR - lpR;

    	tremoloModL = ((lfoTremoloSin 		* fxTremoloDepth + 1 - fxTremoloDepth) * 2 - 1) * tremoloEnvFollowModAttn;
    	tremoloModR = (((1 - lfoTremoloSin) * fxTremoloDepth + 1 - fxTremoloDepth) * 2 - 1) * tremoloEnvFollowModAttn;

        combInL = lpL * tremoloModL + hpL * (1 - tremoloModL);
        combInR = lpR * tremoloModR + hpR * (1 - tremoloModR);

    	// --- enveloppe calculation

        blocksum 	+= fabsf(combInL);
        envelope 	= (envelope * envM1 + envDest) * envM2;
        envMod 		= envModDepth * envelope;

    	// feedback & interpolation

        fbL = feedbackHermiteInterpolation(feedbackReadPosIntL);
        fbR = feedbackHermiteInterpolation(feedbackReadPosIntR);

    	// --- input tremolo

    	vcaL = (combInL) * ((tremoloModL * tremoloPanDepth + (1-tremoloPanDepth)) * 2 - 1);
    	vcaR = (combInR) * ((tremoloModR * tremoloPanDepth + (1-tremoloPanDepth)) * 2 - 1);

    	// --- mix node

    	nodeL =  (vcaL +	fbL	* fxFeedback) * feedbackAttn;
    	nodeR =  (vcaR +	fbR	* fxFeedback) * feedbackAttn;

        // --- low pass

        v2L += feedbackLp * v3L;						// lowpass
        v3L += feedbackLp * ( nodeL - v2L - v3L);
        v2L += feedbackLp * v3L;						// lowpass
        v3L += feedbackLp * ( nodeL - v2L - v3L);

        nodeL = v2L;

        v2R += feedbackLp * v3R;						// lowpass
        v3R += feedbackLp * (nodeR - v2R - v3R);
        v2R += feedbackLp * v3R;						// lowpass
        v3R += feedbackLp * (nodeR - v2R - v3R);

        nodeR = v2R;

    	// --- inject in feedback buffer

    	fxCrossover = 0.5f +  (1 + lfo1) * 0.25f;

    	feedbackBuffer[ feedbackWritePos ] 		= (1-fxCrossover) * nodeL + fxCrossover * nodeR;
    	feedbackBuffer[ feedbackWritePos + 1 ] 	= (1-fxCrossover) * nodeR + fxCrossover * nodeL;

        //

        outL = fxInputLevel * fbL + nodeL;
        outR = fxInputLevel * fbR + nodeR;

    	// --- final allpass

    	_ly1L = coef4L * (_ly1L + outL) - _lx1L; // allpass
    	_lx1L = outL;

    	_ly1R = coef4R * (_ly1R + outR) - _lx1R; // allpass
    	_lx1R = outR;

    	_ly1L = coef4L * (_ly1L + outL) - _lx1L; // allpass
    	_lx1L = outL;

    	_ly1R = coef4R * (_ly1R + outR) - _lx1R; // allpass
    	_lx1R = outR;

    	// --- mix out

    	*(outBuff++) += (int32_t) ( _ly1L * sampleMultipler);
    	*(outBuff++) += (int32_t) ( _ly1R * sampleMultipler);

    	sample += 2;

    	// --- write index increment

    	feedbackWritePos 	+= 2;
    	if( feedbackWritePos >= feedbackBufferSize )
    		feedbackWritePos = 0;

    	// --- lfo increment

    	lfo1tri -= lfo1Inc * fxSpeed;// * (1 + envMod);
    	if(lfo1tri >= 1) {
    		lfo1tri = 1;
    		lfo1Inc = -lfo1Inc;
    	}
    	if(lfo1tri <= -1) {
    		lfo1tri = -1;
    		lfo1Inc = -lfo1Inc;
    	}
    	lfo1 = sigmoid(lfo1tri);

    	lfo2 += lfo2Inc;
    	if(lfo2 >= 1) {
    		lfo2 = 1;
    		lfo2Inc = -lfo2Inc;
    	}
    	if(lfo2 <= -1) {
    		lfo2 = -1;
    		lfo2Inc = -lfo2Inc;
    	}

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
    	timeCv 		= (lfo1 * lfoDepth + envMod) * feedbackFxTarget;
    	cvDelta = -(prevTimeCv - timeCv) * 0.07f;
    	timeCvSpeed += cvDelta;
    	timeCvSpeed *= 0.9f;
    	bouncingCv 	+= timeCvSpeed;

    	timeCvControl = bounceLevel * bouncingCv + (1 - bounceLevel) * timeCv;
    }
}

float FxBus::feedbackHermiteInterpolation(int readPos) {

	float y0 = feedbackBuffer[(readPos + feedbackBufferSize - 2) & (feedbackBufferSize - 1) ];
	float y1 = feedbackBuffer[readPos];
	float y2 = feedbackBuffer[(readPos + 2)& (feedbackBufferSize - 1)];
	float y3 = feedbackBuffer[(readPos + 4)& (feedbackBufferSize - 1)];

	float x = y2 - y1;
	x -= floorf(x);

    return hermite4(x, y0, y1, y2, y3);
}
