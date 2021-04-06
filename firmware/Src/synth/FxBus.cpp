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
inline
float tanh3(float x)
{
    return 1.5f * x / (1.7f + fabsf(0.34f * x * x));
}
inline
float tanh4(float x)
{
    return x / sqrt3(x * x + 1);
}
inline
float sat66(float x)
{
    return x * (1 - (x * x * 0.055f));
}
inline float CubicHermite( float t, float A, float B, float C, float D)
{
    float a = -A * 0.5f + (3.0f * B) * 0.5f - 1.5f * C + D * 0.5f;
    float b = A - (5.0f * B) * 0.5f + 2.0f*C - D * 0.5f;
    float c = -A * 0.5f + C * 0.5f;
    float d = B;

    return a*t*t*t + b*t*t + c*t + d;
}
inline float hermite1(float x, float y0, float y1, float y2, float y3)
{
    // 4-point, 3rd-order Hermite (x-form)
    float c0 = y1;
    float c1 = 0.5f * (y2 - y0);
    float c2 = y0 - 2.5f * y1 + 2.f * y2 - 0.5f * y3;
    float c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);

    return ((c3 * x + c2) * x + c1) * x + c0;
}
//interpolates between L0 and H0 taking the previous (L1) and next (H1) points into account
inline float ThirdInterp(const float x,const float L1,const float L0,const float H0,const float H1)
{
    return
    L0 +
    .5f*
    x*(H0-L1 +
       x*(H0 + L0*(-2) + L1 +
          x*( (H0 - L0)*9 + (L1 - H1)*3 +
             x*((L0 - H0)*15 + (H1 -  L1)*5 +
                x*((H0 - L0)*6 + (L1 - H1)*2 )))));
}
inline
float sigmoid(float x)
{
    return x * (1.5f - 0.5f * x * x);
}
inline
float sigmoidPos(float x)
{
    //x : 0 -> 1
    return (sigmoid((x * 2) - 1) + 1) * 0.5f;
}
inline
float fold(float x4) {
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
float tri2sin(float x) {
	//int index = (int) (x * 2048);
    //return sinTable[index];
	// 0 < x < 1

	return x * x * (3 - x - x);
}
//***------------***------------***------------***------------***----- FxBus -------***------------***------------***------------

float FxBus::forwardBuffer[FxBus::forwardBufferSize] __attribute__((section(".ram_d1")));
float FxBus::feedbackBuffer[FxBus::feedbackBufferSize] __attribute__((section(".ram_d1")));

FxBus::FxBus() {
	lfo1 = 0;
}

void FxBus::init(SynthState *synthState) {
    this->synthState_ = synthState;

    for (int s = 0; s < FxBus::forwardBufferSize; s++) {
    	forwardBuffer[ s ] = 0;
    }

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

	fxLp = 0.573f;
	fxCrossover = 0.337f;
}

/**
 * init before sum timbres
 */
void FxBus::mixSumInit() {
    float temp;
	sample = getSampleBlock();

    for (int s = 0; s < BLOCK_SIZE; s++) {
    	*(sample++) = 0;
    	*(sample++) = 0;
    }

    temp 			=  	clamp( synthState_->fullState.masterfxConfig[ MASTERFX_TIME ], 0.0003f, 0.9997f);
    temp			*= 	temp * temp;
    fxTime 			= 	fxTime * 0.9f + temp * 0.1f;

    fxFeedforward	= 	synthState_->fullState.masterfxConfig[ MASTERFX_FFORWARD ];
    fxFeedback 		= 	synthState_->fullState.masterfxConfig[ MASTERFX_FBACK ];
    fxInputLevel 	= 	synthState_->fullState.masterfxConfig[ MASTERFX_INPUTLEVEL ];

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_TREMOLOENVFOLLOW] * 0.99f;
    tremoloEnvFollow= 	tremoloEnvFollow * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_TREMOLOSPEED];
	fxTremoloSpeed	= 	fxTremoloSpeed * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_TREMOLODEPTH ];
	fxTremoloDepth	= 	fxTremoloDepth * 0.9f + temp * 0.1f;

    fxTone 			= 	fxTime * 0.15f + 0.11f - lfo2 * 0.02f;
    fxDiffusion 	= 	clamp(0.5f + lfo2 * 0.45f, 0.1f, 1);

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_SPEED ];
    temp 			*=	temp * temp;
    fxSpeed 		= 	fxSpeed * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_MOD ];
    temp 			= 	temp * (0.1f + (1 - sqrt3(synthState_->fullState.masterfxConfig[ MASTERFX_SPEED ])) * 0.25f);
    fxMod 			= 	fxMod * 0.9f + temp * 0.1f;

    forwardFxTarget  = 	clamp(getQuantizedTime(fxTime, forwardBufferSize), 2, forwardBufferSize - 8);
    forwardDelayLen  = 	forwardDelayLen 	+ ( forwardFxTarget -  forwardDelayLen)	* 0.01f;

    feedbackFxTarget = 	clamp(getQuantizedTime(fxTime, feedbackBufferSize), 2, feedbackBufferSize - 8);
    feedbackDelayLen = 	feedbackDelayLen 	+ (feedbackFxTarget - feedbackDelayLen)	* 0.01f;

    // ----------- VCF -----------

	float OffsetTmp = fxDiffusion;
	vcfDiffusion = clamp((OffsetTmp + 9.0f * vcfDiffusion) * .1f, filterWindowMin, filterWindowMax);

	const float bipolarf = (vcfFreq - 0.5f);
	const float folded = fold(sigmoid(bipolarf * 19 * vcfDiffusion)) * 4; // -1 < folded < 1
	const float offset = vcfDiffusion * vcfDiffusion * 0.17f;
	const float lrDelta = 0.0095f * folded;
	const float range = 0.47f + lfo2 * 0.1f;

	f1L = clamp(((vcfFreq - offset - lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
    f2L = clamp(((vcfFreq + offset + lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
	f3L = clamp(((vcfFreq - (offset * 2) - lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
    f4L = clamp(((vcfFreq + (offset * 2) + lrDelta) * range) * 2, filterWindowMin, filterWindowMax);

	coef1L = (1.0f - f1L) / (1.0f + f1L);
    coef2L = (1.0f - f2L) / (1.0f + f2L);
	coef3L = (1.0f - f3L) / (1.0f + f3L);
    coef4L = (1.0f - f4L) / (1.0f + f4L);

	f1R = clamp(((vcfFreq - offset + lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
    f2R = clamp(((vcfFreq + offset - lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
	f3R = clamp(((vcfFreq - (offset * 2) + lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
    f4R = clamp(((vcfFreq + (offset * 2) - lrDelta) * range) * 2, filterWindowMin, filterWindowMax);

	coef1R = (1.0f - f1R) / (1.0f + f1R);
    coef2R = (1.0f - f2R) / (1.0f + f2R);
	coef3R = (1.0f - f3R) / (1.0f + f3R);
	coef4R = (1.0f - f4R) / (1.0f + f4R);

    // ----------- /vcf -----------

	inHpF = 0.0125f;
	inLpF = 0.56f;
}

/**
 * add timbre block to bus mix
 */
void FxBus::mixSum(float *inStereo, int timbreNum) {
	const float level = synthState_->mixerState.instrumentState_[timbreNum].send;

	sample = getSampleBlock();
	for (int s = 0; s < (BLOCK_SIZE); s++) {
    	*(sample++) += *inStereo++ * level;
    	*(sample++) += *inStereo++ * level;
    }
}

/**
 * process fx on bus mix
 */
void FxBus::processBlock(int32_t *outBuff) {
	sample = getSampleBlock();
    const float sampleMultipler = (float) 0x34ffff; // fx level , max = 0x7fffff

	float fwL, fwR;
	float fbL, fbR;
	float lpR, lpL;
	float hpR, hpL;
	float inR, inL;
	float feedfwAttn 	= 0.83f;
	float feedbackAttn 	= 0.73f;
	float attn = 0.985f;

    for (int s = 0; s < BLOCK_SIZE; s++) {

    	forwardReadPos = forwardWritePos - forwardDelayLen + feedMod();

    	while( forwardReadPos < 0 )
    		forwardReadPos += forwardBufferSize;
    	while( forwardReadPos >= forwardBufferSize )
    		forwardReadPos -= forwardBufferSize;

        forwardReadPosInt = (int) forwardReadPos;
        forwardReadPosInt &= 0xfffffffe; // keep it even for stereo separation

    	feedbackReadPosL = feedbackWritePos - feedbackDelayLen + fdbckMod();

    	while( feedbackReadPosL < 0 )
    		feedbackReadPosL += feedbackBufferSize;
    	while( feedbackReadPosL >= feedbackBufferSize )
    		feedbackReadPosL -= feedbackBufferSize;

    	feedbackReadPosIntL = (int) feedbackReadPosL;
    	feedbackReadPosIntL &= 0xfffffffe; // keep it even for stereo separation

    	feedbackReadPosR = feedbackWritePos - feedbackDelayLen - fdbckMod();

    	while( feedbackReadPosR < 0 )
    		feedbackReadPosR += feedbackBufferSize;
    	while( feedbackReadPosR >= feedbackBufferSize )
    		feedbackReadPosR -= feedbackBufferSize;

    	feedbackReadPosIntR = ((int) feedbackReadPosR)&0xfffffffe;
    	feedbackReadPosIntR += 1;

    	// --- cut low / high

    	inR = *(sample);
    	inL = *(sample + 1);

        v4L = v4L + inHpF * v5L;
        hpL = inR - v4L - v5L;
        v5L = inHpF * hpL + v5L;

        v6L += inLpF * v7L;
        v7L += inLpF * (hpL - v6L - v7L);
        v6L += inLpF * v7L;
        v7L += inLpF * (hpL - v6L - v7L);

        lpL = (v6L * attn);


        v4R = v4R + inHpF * v5R;
        hpR = inL - v4R - v5R;
        v5R = inHpF * hpR + v5R;

        v6R += inLpF * v7R;
        v7R += inLpF * (hpR - v6R - v7R);
        v6R += inLpF * v7R;
        v7R += inLpF * (hpR - v6R - v7R);

        lpR = (v6R * attn);

    	// --- enveloppe follower

    	float tmp = fabsf(inR);
    	if(tmp > envelope)
    	    envelope = attack_coef 	* (envelope - tmp) + tmp;
    	else
    	    envelope = release_coef * (envelope - tmp) + tmp;

    	// --- vcf 1

    	vcf1(forwardReadPosInt);

    	// --- feed forward

    	fwL = forwardHermiteInterpolation(forwardReadPosInt);
    	fwR = forwardHermiteInterpolation(forwardReadPosInt + 1);

        forwardBuffer[ forwardWritePos ] 		= clamp( (lpL 	+ fwL * fxFeedforward) * feedfwAttn, -1, 1);
    	forwardBuffer[ forwardWritePos + 1 ] 	= clamp( (lpR	+ fwR * fxFeedforward) * feedfwAttn, -1, 1);

    	// --- vcf 2

    	vcf2L(feedbackReadPosIntL);
    	vcf2R(feedbackReadPosIntR);

    	//

    	fbL = feedbackHermiteInterpolation(feedbackReadPosIntL) + fxCrossover * feedbackHermiteInterpolation(feedbackReadPosIntR);
    	fbR = feedbackHermiteInterpolation(feedbackReadPosIntR) + fxCrossover * feedbackHermiteInterpolation(feedbackReadPosIntL);

    	// --- feedback

    	feedbackBuffer[ feedbackWritePos ] 		=
    			clamp( (lpL * fxInputLevel +	fwL + fbL		* fxFeedback) * feedbackAttn  , -1, 1);

    	feedbackBuffer[ feedbackWritePos + 1 ] 	=
    			clamp( (lpR * fxInputLevel +	fwR + fbR	 	* fxFeedback) * feedbackAttn  , -1, 1);

    	// --- final allpass

    	_ly4L = coef4L * (_ly4L + feedbackBuffer[ feedbackWritePos ]) - _lx4L; // allpass
    	_lx4L = feedbackBuffer[ feedbackWritePos ];

    	_ly4R = coef4R * (_ly4R + feedbackBuffer[ feedbackWritePos+1 ]) - _lx4R; // allpass
    	_lx4R = feedbackBuffer[ feedbackWritePos+1 ];

    	// --- mix out

    	float tremoloModL = ((	tri2sin(	lfoTremolo) 	* fxTremoloDepth + 1 - fxTremoloDepth) * 2 - 1);
    	float tremoloModR = ((	tri2sin(1 - lfoTremolo) 	* fxTremoloDepth + 1 - fxTremoloDepth) * 2 - 1);

    	*(outBuff++) += (int32_t) ( _ly4L * sampleMultipler * tremoloModL);
    	*(outBuff++) += (int32_t) ( _ly4R * sampleMultipler * tremoloModR);

    	sample += 2;

    	// --- write index increment

    	forwardWritePos 	+= 2;
    	if( forwardWritePos >= forwardBufferSize )
    		forwardWritePos = 0;

    	feedbackWritePos 	+= 2;
    	if( feedbackWritePos >= feedbackBufferSize )
    		feedbackWritePos = 0;

    	// --- lfo increment

    	lfo1 += lfo1Inc * fxSpeed;
    	if(lfo1 >= 1) {
    		lfo1 = 1;
    		lfo1Inc = -lfo1Inc;
    	}
    	if(lfo1 <= -1) {
    		lfo1 = -1;
    		lfo1Inc = -lfo1Inc;
    	}

    	lfo2 += lfo2Inc;
    	if(lfo2 >= 1) {
    		lfo2 = 1;
    		lfo2Inc = -lfo2Inc;
    	}
    	if(lfo2 <= -1) {
    		lfo2 = -1;
    		lfo2Inc = -lfo2Inc;
    	}

    	float envFollowMod = 1 + tremoloEnvFollow * envelope;
    	lfoTremolo += lfoTremoloInc * fxTremoloSpeed * envFollowMod;
    	if(lfoTremolo >= 1 ) {
    		lfoTremolo = 1;
    		lfoTremoloInc = -lfoTremoloInc;
    	}
    	if(lfoTremolo <= 0) {
    		lfoTremolo = 0;
    		lfoTremoloInc = -lfoTremoloInc;
    	}
    }

}

float FxBus::feedMod() {
	return ( lfo1 * fxMod * forwardFxTarget * 0.01f);
}
float FxBus::fdbckMod() {
	return ( lfo1 * fxMod * feedbackFxTarget );
}
float FxBus::forwardHermiteInterpolation(int readPos) {
	float y0 = forwardBuffer[(readPos + forwardBufferSize - 2) % forwardBufferSize];
	float y1 = forwardBuffer[readPos];
	float y2 = forwardBuffer[(readPos + forwardBufferSize + 2) % forwardBufferSize];
	float y3 = forwardBuffer[(readPos + forwardBufferSize + 4) % forwardBufferSize];

	float x = y2 - y1;
	x -= floor(x);

    return hermite1(x, y0, y1, y2, y3);
}
float FxBus::feedbackHermiteInterpolation(int readPos) {

	float y0 = feedbackBuffer[(readPos + feedbackBufferSize - 2) % feedbackBufferSize];
	float y1 = feedbackBuffer[readPos];
	float y2 = feedbackBuffer[(readPos + feedbackBufferSize + 2) % feedbackBufferSize];
	float y3 = feedbackBuffer[(readPos + feedbackBufferSize + 4) % feedbackBufferSize];

	float x = y2 - y1;
	x -= floor(x);

    return hermite1(x, y0, y1, y2, y3);
}
void FxBus::vcf1(int readPos) {

	float inmix, hpL, hpR;
    const float f = fxLp + 0.05f;

	// Left voice
	inmix = forwardBuffer[readPos];

    v0L = v0L + f * v1L;
    hpL = inmix - v0L - v1L;
    v1L = f * hpL + v1L;

    forwardBuffer[readPos] = v0L;

	// Right voice
	inmix = forwardBuffer[readPos + 1];

    v0R = v0R + f * v1R;
    hpR = inmix - v0R - v1R;
    v1R = f * hpR + v1R;

    forwardBuffer[readPos + 1] = v0R;
}

void FxBus::vcf2L(int readPos) {
	// Left voice
	float inmix = feedbackBuffer[readPos];

    v2L += fxLp * v3L;						// lowpass
    v3L += fxLp * ( (inmix) - v2L - v3L);

    v2L += fxLp * v3L;						// lowpass
    v3L += fxLp * ( (inmix) - v2L - v3L);

    feedbackBuffer[readPos] = (v2L);
}
void FxBus::vcf2R(int readPos) {

	// Right voice
    float inmix = feedbackBuffer[readPos];

    v2R += fxLp * v3R;						// lowpass
    v3R += fxLp * ((inmix) - v2R - v3R);

    v2R += fxLp * v3R;						// lowpass
    v3R += fxLp * ((inmix) - v2R - v3R);


    feedbackBuffer[readPos] = (v2R);
}
