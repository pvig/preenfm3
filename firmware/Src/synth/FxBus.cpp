/*
 * FxBus.cpp
 *
 *  Created on: Feb 14, 2021
 *      Author: patvig
 */

#include "FxBus.h"

#define filterWindowMin 0.01f
#define filterWindowMax 0.99f

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
    v0L = v1L = v2L = v3L = v4L = v5L = v6L = v7L = v8L = v0R = v1R = v2R = v3R = v4R = v5R = v6R = v7R = v8R = v8R = 0.0f;
	vcfFreq = 0.33f;
	vcfDiffusion = 0.47f;
}

/**
 * init before sum timbres
 */
void FxBus::mixSumInit() {

	sample = getSampleBlock();
    for (int s = 0; s < BLOCK_SIZE; s++) {
    	*(sample++) = 0;
    	*(sample++) = 0;
    }
    float temp;
    temp 		=  clamp( synthState_->fullState.masterfxConfig[ MASTERFX_TIME ], 0.0001f, 0.9999f);
    temp		*= temp * temp;
    fxTime 		= fxTime * 0.9f + temp * 0.1f;

    fxFeedforward	= synthState_->fullState.masterfxConfig[ MASTERFX_FFORWARD ] * 0.85f;
    fxFeedback 		= synthState_->fullState.masterfxConfig[ MASTERFX_FBACK ];
    fxInputLevel 	= synthState_->fullState.masterfxConfig[ MASTERFX_INPUTLEVEL ];

    fxTone 		= fxTime * 0.15f + 0.11f;
    fxDiffusion = fxTime * 0.10f + 0.11f;

    temp 	= 	synthState_->fullState.masterfxConfig[ MASTERFX_MOD ];
    temp 	*= 	temp;
    fxMod 	= 	fxMod * 0.9f + temp * 0.1f;

    fxSpeed 	= synthState_->fullState.masterfxConfig[ MASTERFX_SPEED ];
    fxSpeed		= fxSpeed * fxSpeed * fxSpeed * 6;

    forwardDelayLen 	= forwardDelayLen 	* 0.97f 	+ fxTime 	* forwardBufferSize 	* 0.03f;
    feedbackDelayLen 	= feedbackDelayLen 	* 0.97f 	+ fxTime 	* feedbackBufferSize 	* 0.03f;

    if (forwardDelayLen  != prevFxForward) {
        forwardReadPos = forwardWritePos - forwardDelayLen;
    	if( forwardReadPos < 0 )
    		forwardReadPos += forwardBufferSize;
        //forwardReadPos = ((int)forwardReadPos)&0xfffffffe; // make it even
    }

    if (feedbackDelayLen  != prevFxFeedback) {
    	feedbackReadPos = feedbackWritePos - feedbackDelayLen;
    	if( feedbackReadPos < 0 )
    		feedbackReadPos += feedbackBufferSize;
    	//feedbackReadPos = ((int)feedbackReadPos )&0xfffffffe; // make it even
    }

    prevFxForward = forwardDelayLen;
    prevFxFeedback = feedbackDelayLen;
}

/**
 * add timbre to bus mix
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
    const float sampleMultipler = (float) 0x37ffff; // fx level , max = 0x7fffff

    // ----------- VCF -----------

	float fxParamTmp = fxTone;
	fxParamTmp *= fxParamTmp;
	vcfFreq = clamp((fxParamTmp + 9.0f * vcfFreq) * .1f, filterWindowMin, 0.9f);

	float OffsetTmp = fxDiffusion;
	vcfDiffusion = clamp((OffsetTmp + 9.0f * vcfDiffusion) * .1f, filterWindowMin, filterWindowMax);

	const float bipolarf = (vcfFreq - 0.5f);
	const float folded = fold(sigmoid(bipolarf * 19 * vcfDiffusion)) * 4; // -1 < folded < 1
	const float offset = vcfDiffusion * vcfDiffusion * 0.17f;
	const float lrDelta = 0.005f * folded;
	const float range = 0.47f;

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

	_ly1L = v0L;
	_ly1R = v0R;
	_ly2L = v1L;
	_ly2R = v1R;
	_ly3L = v2L;
	_ly3R = v2R;
	_ly4L = v3L;
	_ly4R = v3R;
	_lx1L = v4L;
	_lx1R = v4R;
	_lx2L = v5L;
	_lx2R = v5R;
	_lx3L = v6L;
	_lx3R = v6R;
	_lx4L = v7L;
	_lx4R = v7R;

    // ----------- /vcf -----------

	float fwL, fwR;

    for (int s = 0; s < BLOCK_SIZE; s++) {

    	// --- feed forward

    	fwL = forwardBufferInterpolation(forwardReadPos, false);
    	fwR = forwardBufferInterpolation(forwardReadPos, true);

        forwardBuffer[ forwardWritePos ] 		= *(sample) 	+ tanh4(  	( fwL	* fxFeedforward) ) * 0.85f;
    	forwardBuffer[ forwardWritePos + 1 ] 	= *(sample+1) 	+ tanh4(  	( fwR 	* fxFeedforward) ) * 0.85f;

    	// --- vcf 1

    	vcf1(forwardReadPosInt);

        // --- feedback

    	feedbackBuffer[ feedbackWritePos ] =
    			tanh4( *(sample++) * fxInputLevel +	fwL + (feedbackBufferInterpolation(feedbackReadPos, false) 	* fxFeedback)) * 0.975f;

    	feedbackBuffer[ feedbackWritePos + 1 ] =
    			tanh4( *(sample++) * fxInputLevel +	fwR + (feedbackBufferInterpolation(feedbackReadPos, true) 	* fxFeedback)) * 0.975f;

    	// --- vcf 2

    	vcf2(feedbackWritePos&0xfffffffe);

    	// mix out

    	*(outBuff++) += (int32_t) ( feedbackBuffer[ feedbackWritePos ] 		* sampleMultipler);
    	*(outBuff++) += (int32_t) ( feedbackBuffer[ feedbackWritePos + 1 ] 	* sampleMultipler);


    	// --- feed forward

    	forwardReadPos 		+= 2  + ( lfo1 * fxMod * 0.125f );
    	forwardWritePos 	+= 2;

    	if( forwardWritePos >= forwardBufferSize )
    		forwardWritePos -= forwardBufferSize;

    	while( forwardReadPos < 0 )
    		forwardReadPos += forwardBufferSize;

    	while( forwardReadPos >= forwardBufferSize )
    		forwardReadPos -= forwardBufferSize;

        forwardReadPosInt = (int) forwardReadPos;
        forwardReadPosInt &= 0xfffffffe; // make it even

        // --- feedback

    	feedbackReadPos 	+= 2 + ( lfo1 * fxMod );
    	feedbackWritePos 	+= 2;

    	if( feedbackWritePos >= feedbackBufferSize )
    		feedbackWritePos -= feedbackBufferSize;

    	while( feedbackReadPos < 0 )
    		feedbackReadPos += feedbackBufferSize;

    	while( feedbackReadPos >= feedbackBufferSize )
    		feedbackReadPos -= feedbackBufferSize;

    	// --- lfo

    	lfo1 += lfo1Inc * fxSpeed;
    	if(lfo1 >= 1 || lfo1 <= -1) {
    		lfo1Inc = -lfo1Inc;
    	}

    	/*lfo2 += lfo2Inc;
    	if(lfo2 >= 1 || lfo2 <= -1) {
    		lfo2Inc = -lfo2Inc;
    	}

    	lfo3 += lfo3Inc;
    	if(lfo3 >= 1 || lfo3 <= -1) {
    		lfo3Inc = -lfo3Inc;
    	}

    	lfo4 += lfo4Inc;
    	if(lfo4 >= 1 || lfo4 <= -1) {
    		lfo4Inc = -lfo4Inc;
    	}*/

    }

	v0L = _ly1L;
	v0R = _ly1R;
	v1L = _ly2L;
	v1R = _ly2R;
	v2L = _ly3L;
	v2R = _ly3R;
	v3L = _ly4L;
	v3R = _ly4R;
	v4L = _lx1L;
	v4R = _lx1R;
	v5L = _lx2L;
	v5R = _lx2R;
	v6L = _lx3L;
	v6R = _lx3R;
	v7L = _lx4L;
	v7R = _lx4R;

}

float FxBus::forwardBufferInterpolation(float readPos, bool isRight) {
	float x = readPos - floor(readPos);
	int posInt = ((int) readPos)&0xfffffffe;

	if(isRight) {
		posInt += 1;
	}

	float y0 = forwardBuffer[posInt];
	float y1 = forwardBuffer[(posInt + forwardBufferSize - 2) % forwardBufferSize];

	return y0 * (1 - x) + y1 * x;
}

float FxBus::feedbackBufferInterpolation(float readPos, bool isRight) {
	float x = readPos - floor(readPos);
	int posInt = ((int) readPos)&0xfffffffe;

	if(isRight) {
		posInt += 1;
	}

	float y0 = feedbackBuffer[posInt];
	float y1 = feedbackBuffer[(posInt + feedbackBufferSize - 2) % feedbackBufferSize];

	return y0 * (1 - x) + y1 * x;
}

void FxBus::vcf1(int readPos) {

	const float _feedback = 0.43f;
	const float _crossFeedback = 0.15f;
	float inmix;
	float sp = forwardBuffer[readPos];

	// Left voice
	inmix = sp - _feedback * _ly2L + _crossFeedback * _ly2R;

	_ly1L = coef1L * (_ly1L + inmix) - _lx1L; // do 1st filter
	_lx1L = inmix;
	_ly2L = coef2L * (_ly2L + _ly1L) - _lx2L; // do 2nd filter
	_lx2L = _ly1L;

    forwardBuffer[readPos] = sp + _ly2L;

	sp = forwardBuffer[readPos + 1];

	// Right voice
	inmix = sp - _feedback * _ly2R + _crossFeedback * _ly2L;

	_ly1R = coef1R * (_ly1R + inmix) - _lx1R; // do 1st filter
	_lx1R = inmix;
	_ly2R = coef2R * (_ly2R + _ly1R) - _lx2R; // do 2nd filter
	_lx2R = _ly1R;

    forwardBuffer[readPos + 1] = sp + _ly2R;
}


void FxBus::vcf2(int readPos) {

	const float _feedback = 0.43f;
	const float _crossFeedback = 0.15f;
	float inmix;
	float sp = feedbackBuffer[readPos];

	// Left voice
	inmix = sp - _feedback * _ly4L + _crossFeedback * _ly4R;

	_ly3L = coef3L * (_ly3L + inmix) - _lx3L; // do 3rd filter
	_lx3L = _ly2L;
    _ly4L = coef4L * (_ly4L + _ly3L) - _lx4L; // do 4nth filter
    _lx4L = _ly3L;

    feedbackBuffer[readPos] = sp + _ly2L;

	sp = feedbackBuffer[readPos + 1];

	// Right voice
	inmix = sp - _feedback * _ly4R + _crossFeedback * _ly4L;

	_ly3R = coef3R * (_ly3R + inmix) - _lx3R; // do 3rd filter
	_lx3R = _ly2R;
    _ly4R = coef4R * (_ly4R + _ly3R) - _lx4R; // do 4nth filter
    _lx4R = _ly3R;

    feedbackBuffer[readPos + 1] = sp + _ly2R;
}
