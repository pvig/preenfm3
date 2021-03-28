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

	vcfFreq = 0.46f;
	vcfDiffusion = 0.13f;
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

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_LP ] * 0.8f;
    temp			*= 	temp * temp;
    fxLp			= 	fxLp * 0.9f + temp* 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_TIMESHIFT ];
	fxTimeShift		= 	fxTimeShift * 0.9f + temp * 0.1f;

    fxTone 			= 	fxTime * 0.15f + 0.11f - lfo2 * 0.02f;
    fxDiffusion 	= 	clamp(0.5f + lfo2 * 0.45f, 0.1f, 1);

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_SPEED ];
    temp 			*=	temp * temp;
    fxSpeed 		= 	fxSpeed * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_MOD ];
    temp 			= 	temp * (0.1f + (1 - sqrt3(synthState_->fullState.masterfxConfig[ MASTERFX_SPEED ])) * 0.5f);
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
	const float lrDelta = 0.005f * folded;
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

    for (int s = 0; s < BLOCK_SIZE; s++) {

    	forwardReadPos = forwardWritePos - forwardDelayLen + feedMod();

    	while( forwardReadPos < 0 )
    		forwardReadPos += forwardBufferSize;
    	while( forwardReadPos >= forwardBufferSize )
    		forwardReadPos -= forwardBufferSize;

        forwardReadPosInt = (int) forwardReadPos;
        forwardReadPosInt &= 0xfffffffe; // keep it even for stereo separation


    	feedbackReadPos = feedbackWritePos - feedbackDelayLen + fdbckMod();

    	while( feedbackReadPos < 0 )
    		feedbackReadPos += feedbackBufferSize;
    	while( feedbackReadPos >= feedbackBufferSize )
    		feedbackReadPos -= feedbackBufferSize;

    	feedbackReadPosInt = (int) feedbackReadPos;
    	feedbackReadPosInt &= 0xfffffffe; // keep it even for stereo separation

    	// --- vcf 1

    	vcf1(forwardReadPosInt);

    	// --- feed forward

    	fwL = forwardHermiteInterpolation(forwardReadPosInt);
    	fwR = forwardHermiteInterpolation(forwardReadPosInt + 1);

        forwardBuffer[ forwardWritePos ] 		= clamp(*(sample) 		+ tanh4( fwL * fxFeedforward) * 0.87f, -1, 1);
    	forwardBuffer[ forwardWritePos + 1 ] 	= clamp(*(sample+1) 	+ tanh4( fwR * fxFeedforward) * 0.87f, -1, 1);

    	// --- vcf 2

    	vcf2(feedbackReadPosInt);

        // --- feedback

    	feedbackBuffer[ feedbackWritePos ] 		=
    			tanh4( *(sample++) * fxInputLevel +	fwL + clamp(feedbackHermiteInterpolation(feedbackReadPosInt) 		* fxFeedback, -1, 1)) * 0.88f;

    	feedbackBuffer[ feedbackWritePos + 1 ] 	=
    			tanh4( *(sample++) * fxInputLevel +	fwR + clamp(feedbackHermiteInterpolation(feedbackReadPosInt + 1) 	* fxFeedback, -1, 1)) * 0.88f;

    	// mix out

    	*(outBuff++) += (int32_t) ( feedbackBuffer[ feedbackWritePos ] 		* sampleMultipler);
    	*(outBuff++) += (int32_t) ( feedbackBuffer[ feedbackWritePos + 1 ] 	* sampleMultipler);

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
    }

}

float FxBus::feedMod() {
	return ( lfo1 * fxMod * forwardFxTarget);
}
float FxBus::fdbckMod() {
	return ( lfo1 * fxMod * feedbackFxTarget  * (1 + fxTimeShift) );
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

	float inmix;
    const float f = 0.3f;//fxLp;

	float lowL = v0L, bandL = v1L;
	float lowR = v0R, bandR = v1R;

	// Left voice
	inmix = forwardBuffer[readPos];

    lowL += f * bandL;
    bandL += f * ( (inmix) - lowL - bandL);

    lowL += f * bandL;
    bandL += f * ( (inmix) - lowL - bandL);

    _ly1L = coef1L * (_ly1L + (inmix - lowL)) - _lx1L; // allpass
	_lx1L = (inmix - lowL);

	_ly2L = coef2L * (_ly2L + _ly1L) - _lx2L; // allpass
	_lx2L = _ly1L;

    forwardBuffer[readPos] = _ly2L;

	// Right voice
	inmix = forwardBuffer[readPos + 1];

    lowR += f * bandR;
    bandR += f * ((inmix) - lowR - bandR);

    lowR += f * bandR;
    bandR += f * ((inmix) - lowR - bandR);

	_ly1R = coef1R * (_ly1R + (inmix - lowR)) - _lx1R; // allpass
	_lx1R = (inmix - lowR);

	_ly2R = coef2R * (_ly2R + _ly1R) - _lx2R; // allpass
	_lx2R = _ly1R;

    forwardBuffer[readPos + 1] = _ly2R;

    v0L = lowL;
    v1L = bandL;
    v0R = lowR;
    v1R = bandR;
}

void FxBus::vcf2(int readPos) {
	float inmix;
	const float f = fxLp;

	float lowL = v2L, bandL = v3L;
	float lowR = v2R, bandR = v3R;

	// Left voice
	inmix = feedbackBuffer[readPos];

    lowL += f * bandL;
    bandL += f * ( (inmix) - lowL - bandL);

    lowL += f * bandL;
    bandL += f * ( (inmix) - lowL - bandL);

	_ly3L = coef3L * (_ly3L + lowL) - _lx3L; // allpass
	_lx3L = lowL;

	_ly4L = coef4L * (_ly4L + _ly3L) - _lx4L; // allpass
	_lx4L = _ly3L;

    feedbackBuffer[readPos] = (_ly4L);

	// Right voice
    inmix = feedbackBuffer[readPos + 1];

    lowR += f * bandR;
    bandR += f * ((inmix) - lowR - bandR);

    lowR += f * bandR;
    bandR += f * ((inmix) - lowR - bandR);

	_ly3R = coef3R * (_ly3R + lowR) - _lx3R; // allpass
	_lx3R = lowR;

	_ly4R = coef4R * (_ly4R + _ly3R) - _lx4R; // allpass
	_lx4R = _ly3R;

    feedbackBuffer[readPos + 1] = (_ly4R);

    v2L = lowL;
    v3L = bandL;
    v2R = lowR;
    v3R = bandR;
}
