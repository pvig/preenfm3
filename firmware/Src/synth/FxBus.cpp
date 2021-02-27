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
float tanh3(float x)
{
    return 1.5f * x / (1.7f + fabsf(0.34f * x * x));
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
float expf_fast(float a) {
  //https://github.com/ekmett/approximate/blob/master/cbits/fast.c
  union { float f; int x; } u;
  u.x = (int) (12102203 * a + 1064866805);
  return u.f;
}

//***------------***------------***------------***------------***----- FxBus -------***------------***------------***------------

float FxBus::earlyEchoBuffer[FxBus::earlyEchoBufferSize] __attribute__((section(".ram_d1")));
float FxBus::recirculatingBuffer[FxBus::recirculatingBufferSize] __attribute__((section(".ram_d1")));

FxBus::FxBus() {
	lfo1 = 0;
	// https://www.musicdsp.org/en/latest/Effects/44-delay-time-calculation-for-reverberation.html
	/*float t1 = 40.0;        // d0 time
	float g1 = 0.75f;        // d0 gain+
	float gsign = -1;
	float invearlyEchoNum = 1 / earlyEchoNum;*/

	for (int n = 0; n < earlyEchoNum; ++n)
	{
	  /*int dt = 4 * (t1 / powf (2, (float (n) * invearlyEchoNum)));
	  float g = g1 - sqrt3(float (n)/(3+earlyEchoNum));
	  //gsign = - gsign;
	  earlyEchoTimeList[n] = dt &0xfffffffe;
	  earlyEchoGainList[n] = g;*/
		//init delay time : ms to sample
		//earlyEchoTimeList[n] *= 4.8f;
	}

}

void FxBus::init(SynthState *synthState) {
    this->synthState_ = synthState;

    //delayReadPos = earlyEchoSampleCount - delayWritePos;
    for (int s = 0; s < FxBus::earlyEchoBufferSize; s++) {
    	earlyEchoBuffer[ s ] = 0;
    }

	//delayCircReadPos = recirculatingMaxSampleCount - delayCircWritePos;
    for (int s = 0; s < FxBus::recirculatingBufferSize; s++) {
    	recirculatingBuffer[ s ] = 0;
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

    fxType =  synthState_->fullState.masterfxConfig[MASTERFX_TYPE];
    fxTime =  clamp( synthState_->fullState.masterfxConfig[MASTERFX_TIME], 0.003f, 0.999f);
    fxFeedback =  synthState_->fullState.masterfxConfig[MASTERFX_SPACE] * 0.495f;
    fxTone =  synthState_->fullState.masterfxConfig[MASTERFX_TONE];
    fxDiffusion =  synthState_->fullState.masterfxConfig[MASTERFX_DIFFUSION];
    fxWidth =  synthState_->fullState.masterfxConfig[MASTERFX_WIDTH];

    earlyEchoSampleCount = earlyEchoSampleCount * 0.995f + fxTime * FxBus::earlyEchoBufferSize * 0.005f;//  + lfo1 * 0.03f;
    recirculatingMaxSampleCount = recirculatingMaxSampleCount * 0.995f + fxTime * FxBus::recirculatingBufferSize * 0.005f;

    delayReadPos = delayWritePos - earlyEchoSampleCount;
	if( delayReadPos < 0 )
		delayReadPos += FxBus::earlyEchoBufferSize;

	delayCircReadPos = delayCircWritePos - recirculatingMaxSampleCount;
	if( delayCircReadPos < 0 )
		delayCircReadPos += FxBus::recirculatingBufferSize;
}

/**
 * add timbre to bus mix
 */
void FxBus::mixSum(float *inStereo, int timbreNum) {
	float level = synthState_->mixerState.instrumentState_[timbreNum].send;

	sample = getSampleBlock();

	for (int s = 0; s < BLOCK_SIZE; s++) {
    	*(sample++) += *inStereo++ * level;
    	*(sample++) += *inStereo++ * level;
    }
}

/**
 * process fx on bus mix
 */
void FxBus::processBlock(int32_t *outBuff) {
	sample = getSampleBlock();
    float sampleMultipler = (float) 0x7fffff;

    // ----------- VCF -----------

	float fxParamTmp = FxBus::fxTone;
	fxParamTmp *= fxParamTmp;
	vcfFreq = clamp((fxParamTmp + 9.0f * vcfFreq) * .1f, filterWindowMin, 0.9f);

	float OffsetTmp = FxBus::fxDiffusion;
	vcfDiffusion = clamp((OffsetTmp + 9.0f * vcfDiffusion) * .1f, filterWindowMin, filterWindowMax);

	const float bipolarf = (vcfFreq - 0.5f);
	const float folded = fold(sigmoid(bipolarf * 19 * vcfDiffusion)) * 4; // -1 < folded < 1

	const float offset = vcfDiffusion * vcfDiffusion * 0.17f;
	const float lrDelta = 0.005f * folded;

	const float range = 0.47f;

	const float f1L = clamp(fold((vcfFreq - offset - lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
	const float f2L = clamp(fold((vcfFreq + offset + lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
	const float f3L = clamp(fold((vcfFreq - (offset * 2) - lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
	const float f4L = clamp(fold((vcfFreq + (offset * 2) + lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
	coef1L = (1.0f - f1L) / (1.0f + f1L);
	coef2L = (1.0f - f2L) / (1.0f + f2L);
	coef3L = (1.0f - f3L) / (1.0f + f3L);
	coef4L = (1.0f - f4L) / (1.0f + f4L);
	const float f1R = clamp(fold((vcfFreq - offset + lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
	const float f2R = clamp(fold((vcfFreq + offset - lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
	const float f3R = clamp(fold((vcfFreq - (offset * 2) + lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
	const float f4R = clamp(fold((vcfFreq + (offset * 2) - lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
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

	float dL, dR;
	float rL, rR;
	int lineDelay;
	float lineGain, lineFeedback;
	int crossSign, crossDest1, crossDest2;
	float feed1, feed2;

    for (int s = 0; s < BLOCK_SIZE; s++) {

    	// set limits

    	delayReadPos += earlySpeed;
    	delayWritePos += earlySpeed;

    	if( delayWritePos >= earlyEchoBufferSize )
    		delayWritePos -= earlyEchoBufferSize;

    	if( delayReadPos >= earlyEchoBufferSize )
    		delayReadPos -= earlyEchoBufferSize;

        delayReadPosInt = (int) delayReadPos;
        delayReadPosInt &= 0xfffffffe;//make it even



    	// early echo in

        earlyEchoBuffer[ delayWritePos ] = *(sample++);
    	earlyEchoBuffer[ delayWritePos + 1 ] = *(sample++);

        // early echo generators

    	dL = 0;
		dR = 0;

    	for (int n = 0; n < earlyEchoNum; ++n)
    	{
    		lineGain = earlyEchoGainList[n];
    		lineFeedback = earlyEchoFeebackList[n];
    		lineDelay = delayReadPosInt - 56 +  earlyEchoTimeList[n] + (lfo1 * 10);

        	if( lineDelay < 0 )
        		lineDelay += earlyEchoBufferSize;
        	if( lineDelay >= earlyEchoBufferSize )
        		lineDelay -= earlyEchoBufferSize;

        	lineDelay &= 0xfffffffe;//make it even

    		earlyEchoBuffer[delayWritePos] += earlyEchoBuffer[lineDelay] * lineFeedback;
    		earlyEchoBuffer[delayWritePos + 1] += earlyEchoBuffer[lineDelay + 1] * lineFeedback;

    		dL += earlyEchoBuffer[lineDelay] * lineGain;
    		dR += earlyEchoBuffer[lineDelay + 1] * lineGain;


    	}

    	// filter

    	vcf(delayCircReadPosInt);

    	// recirculating delays

        delayCircReadPos += 2;
        delayCircWritePos += 2;

    	if( delayCircReadPos >= recirculatingBufferSize )
    		delayCircReadPos -= recirculatingBufferSize;

    	if( delayCircWritePos >= recirculatingBufferSize )
    		delayCircWritePos -= recirculatingBufferSize;

    	delayCircReadPosInt = (int) delayCircReadPos;
        delayCircReadPosInt &= 0xfffffffe;//make it even



    	rL = dL;
		rR = dR;

    	for (int n = 0; n < recirculatingEchoNum; ++n)
    	{
    		crossSign = (float) recirculatingCrossList[n][0];
    		crossDest1 =recirculatingCrossList[n][1];
    		crossDest2 =recirculatingCrossList[n][2];

    		lineGain = recirculatingGainList[n] * crossSign;
    		lineDelay = delayCircReadPosInt - 96 + recirculatingTimeList[n] + (lfo2 * 10);

        	if( lineDelay < 0 )
        		lineDelay += recirculatingBufferSize;
        	if( lineDelay >= recirculatingBufferSize )
        		lineDelay -= recirculatingBufferSize;

    		lineDelay &= 0xfffffffe;//make it even

    		feed1 = recirculatingBuffer[lineDelay] * lineGain;
    		feed2 = recirculatingBuffer[lineDelay + 1] * lineGain;

    		recirculatingBuffer[ crossDest1 ] += feed1;
    		recirculatingBuffer[ crossDest2 ] += feed2;

    		rL += feed1;
    		rR += feed2;
    	}

    	recirculatingBuffer[delayCircWritePos] = - rL * fxFeedback ;
    	recirculatingBuffer[delayCircWritePos + 1] = - rR * fxFeedback ;


    	// mix out

    	*(outBuff++) += (int32_t) ((dL + rL ) * sampleMultipler);
    	*(outBuff++) += (int32_t) ((dR + rR ) * sampleMultipler);
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


	// --- lfo

	lfo1 += lfo1Inc;
	if(lfo1 >= 1 || lfo1 <= -1) {
		lfo1Inc = -lfo1Inc;
	}

	lfo2 += lfo2Inc;
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
	}

}

float FxBus::interpolation(int delayReadPosInt) {

	float y0 = earlyEchoBuffer[(delayReadPosInt + earlyEchoBufferSize - 2) % earlyEchoBufferSize];
	float y1 = earlyEchoBuffer[(delayReadPosInt + 0)];
	float y2 = earlyEchoBuffer[(delayReadPosInt + 2)];
	float y3 = earlyEchoBuffer[(delayReadPosInt + 4)];

	float x = y2 - y1;
	x -= floor(x);

    return hermite1(x, y0, y1, y2, y3);
}

void FxBus::vcf(int delayReadPosInt) {

	const float _feedback = 0.43f;
	const float _crossFeedback = 0.15f;
	float inmix;
	float sp = recirculatingBuffer[delayReadPosInt];

	// Left voice
	inmix = sp - _feedback * _ly3L + _crossFeedback * _ly3R;

	_ly1L = coef1L * (_ly1L + inmix) - _lx1L; // do 1st filter
	_lx1L = inmix;
	_ly2L = coef2L * (_ly2L + _ly1L) - _lx2L; // do 2nd filter
	_lx2L = _ly1L;
	_ly3L = coef3L * (_ly3L + _ly2L) - _lx3L; // do 3rd filter
	_lx3L = _ly2L;
	_ly4L = coef4L * (_ly4L + _ly3L) - _lx4L; // do 4nth filter
	_lx4L = _ly3L;

	recirculatingBuffer[delayReadPosInt + 1] = sp + _ly4L;

	sp = recirculatingBuffer[delayReadPosInt + 1];

	// Right voice
	inmix = sp - _feedback * _ly3R + _crossFeedback * _ly3L;

	_ly1R = coef1R * (_ly1R + inmix) - _lx1R; // do 1st filter
	_lx1R = inmix;
	_ly2R = coef2R * (_ly2R + _ly1R) - _lx2R; // do 2nd filter
	_lx2R = _ly1R;
	_ly3R = coef3R * (_ly3R + _ly2R) - _lx3R; // do 3rd filter
	_lx3R = _ly2R;
	_ly4R = coef4R * (_ly4R + _ly3R) - _lx4R; // do 4nth filter
	_lx4R = _ly3R;

	recirculatingBuffer[delayReadPosInt + 1] = sp + _ly4R;
}
