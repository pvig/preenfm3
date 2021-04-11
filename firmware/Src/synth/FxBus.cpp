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
float expf_fast(float a) {
  //https://github.com/ekmett/approximate/blob/master/cbits/fast.c
  union { float f; int x; } u;
  u.x = (int) (12102203 * a + 1064866805);
  return u.f;
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
inline float hermite1(float x, float y0, float y1, float y2, float y3)
{
    // 4-point, 3rd-order Hermite (x-form)
    float c0 = y1;
    float c1 = 0.5f * (y2 - y0);
    float c2 = y0 - 2.5f * y1 + 2.f * y2 - 0.5f * y3;
    float c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);

    return ((c3 * x + c2) * x + c1) * x + c0;
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
float tri2sin(float x) {
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

	feedbackLp = 0.573f;
	fxCrossover = 0.337f;
	harmTremoloCutF = 0.4f;

	envFollowLpF = 0.1f;

	inLpF = 0.56f;
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

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_ENVMOD] * 0.99f;
    envModDepth 	= 	envModDepth * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_SPEED ];
    temp 			*=	temp * temp;
    fxSpeed 		= 	fxSpeed * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_MOD ];

    temp 			= 	temp * (0.1f + (1 - sqrt3(synthState_->fullState.masterfxConfig[ MASTERFX_SPEED ])) * 0.35f);
    fxMod 			= 	fxMod * 0.9f + temp * 0.1f;

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_ENVATTACK] * 0.99f;
    envAttackA	 	= 	envAttackA * 0.9f + temp * 0.1f;
    if(envAttackA != envAttackB) {
    	attack_coef 	= expf_fast( -1.0f / (48 * (1 + envAttackA) ));
    	envAttackB = envAttackA;
    }

    temp 			= 	synthState_->fullState.masterfxConfig[ MASTERFX_ENVRELEASE] * 0.99f;
    envReleaseA	 	= 	envReleaseA * 0.9f + temp * 0.1f;
    if(envReleaseA != envReleaseB) {
    	release_coef 	= expf_fast( -1.0f / (48 * (1 + envReleaseA)));
    	envReleaseB = envReleaseA;
    }

    forwardFxTarget  = 	getQuantizedTime(fxTime, forwardBufferSizeReadable);
    forwardDelayLen  = 	forwardDelayLen 	+ ( forwardFxTarget -  forwardDelayLen)	* 0.01f;

    feedbackFxTarget = 	getQuantizedTime(fxTime, feedbackBufferSizeReadable);
    feedbackDelayLen = 	feedbackDelayLen 	+ (feedbackFxTarget - feedbackDelayLen)	* 0.01f;

    // ----------- VCF -----------

	float OffsetTmp = fxDiffusion;
	vcfDiffusion = clamp((OffsetTmp + 9.0f * vcfDiffusion) * .1f, filterWindowMin, filterWindowMax);

	const float bipolarf = (vcfFreq - 0.5f);
	const float folded = fold(sigmoid(bipolarf * 19 * vcfDiffusion)) * 4; // -1 < folded < 1
	const float offset = vcfDiffusion * vcfDiffusion * 0.34f;
	const float lrDelta = 0.0095f * folded;
	const float range = 0.47f + lfo2 * 0.1f;


    f4L = clamp(((vcfFreq + offset + lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
    coef4L = (1.0f - f4L) / (1.0f + f4L);

    f4R = clamp(((vcfFreq + offset - lrDelta) * range) * 2, filterWindowMin, filterWindowMax);
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


	float feedfwAttn 	= 0.93f;
	float feedbackAttn 	= 0.78f;
	float attn = 0.985f;
	float tremoloEnvFollowMod, tremoloEnvFollowModAttn;
	float tremoloModL, tremoloModR;
	float fdbckModVal;
	float fbLInterpol, fbRInterpol;

    for (int s = 0; s < BLOCK_SIZE; s++) {

    	forwardReadPos = forwardWritePos - forwardDelayLen + feedMod();

    	while( forwardReadPos < 0 )
    		forwardReadPos += forwardBufferSizeReadable;
    	while( forwardReadPos >= forwardBufferSizeReadable )
    		forwardReadPos -= forwardBufferSizeReadable;

        forwardReadPosInt = (int) forwardReadPos;
        forwardReadPosInt &= 0xfffffffe; // keep it even for stereo separation

        fdbckModVal = fdbckMod();

    	feedbackReadPosL = feedbackWritePos - feedbackDelayLen + fdbckModVal;

    	while( feedbackReadPosL < 0 )
    		feedbackReadPosL += feedbackBufferSizeReadable;
    	while( feedbackReadPosL >= feedbackBufferSizeReadable )
    		feedbackReadPosL -= feedbackBufferSizeReadable;

    	feedbackReadPosIntL = (int) feedbackReadPosL;
    	feedbackReadPosIntL &= 0xfffffffe; // keep it even for stereo separation

    	feedbackReadPosR = feedbackWritePos - feedbackDelayLen - fdbckModVal;

    	while( feedbackReadPosR < 0 )
    		feedbackReadPosR += feedbackBufferSizeReadable;
    	while( feedbackReadPosR >= feedbackBufferSizeReadable )
    		feedbackReadPosR -= feedbackBufferSizeReadable;

    	feedbackReadPosIntR = ((int) feedbackReadPosR)&0xfffffffe;
    	feedbackReadPosIntR += 1;

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


    	// --- enveloppe follower

        // 		first lowpass input

        v4L = v4L + envFollowLpF * v5L;
        v5L = envFollowLpF * (v6L - v4L - v5L) + v5L;

        float envIn = clamp( fabsf(v4L * 16) , 0, 1);

        if( envIn > envelope  )
        	envelope = envIn + attack_coef * (envelope - envIn);
        else
        	envelope = envIn + release_coef * (envelope - envIn);

        envMod = envelope * envModDepth ;

        tremoloEnvFollowMod = envelope * tremoloEnvFollow;
        tremoloEnvFollowModAttn = (tremoloEnvFollowMod + (1 - fabsf(tremoloEnvFollow)));


    	// --- harmonic tremolo

        _lx3L += harmTremoloCutF * _ly3L;
        _ly3L += harmTremoloCutF * (v6L - _lx3L - _ly3L);

        lpL = (_lx3L * attn);
        hpL = inL - _lx3L;


        _lx3R += harmTremoloCutF * _ly3R;
        _ly3R += harmTremoloCutF * (v6R - _lx3R - _ly3R);

        lpR = (_lx3R * attn);
        hpR = inR - _lx3R;


    	tremoloModL = ((lfoTremoloSin 		* fxTremoloDepth + 1 - fxTremoloDepth) * 2 - 1) * tremoloEnvFollowModAttn;
    	tremoloModR = (((1 - lfoTremoloSin) * fxTremoloDepth + 1 - fxTremoloDepth) * 2 - 1) * tremoloEnvFollowModAttn;

        combInL = lpL * tremoloModL + hpL * (1 - tremoloModL);
        combInR = lpR * tremoloModR + hpR * (1 - tremoloModR);



    	// --- vcf 1

    	vcf1(forwardReadPosInt);

    	// --- feed forward

    	fwL = forwardHermiteInterpolation(forwardReadPosInt);
    	fwR = forwardHermiteInterpolation(forwardReadPosInt + 1);

    	float frwrd = clamp(fxFeedforward , -feedfwAttn, feedfwAttn );

        forwardBuffer[ forwardWritePos ] 		= clamp( (combInL 	+ fwL * frwrd) * feedfwAttn, -1, 1);
    	forwardBuffer[ forwardWritePos + 1 ] 	= clamp( (combInR	+ fwR * frwrd) * feedfwAttn, -1, 1);

    	// --- vcf 2

    	vcf2L(feedbackReadPosIntL);
    	vcf2R(feedbackReadPosIntR);

    	//
    	fbLInterpol = feedbackHermiteInterpolation(feedbackReadPosIntL);
    	fbRInterpol = feedbackHermiteInterpolation(feedbackReadPosIntR);

    	fbL = fbLInterpol + fxCrossover * fbRInterpol;
    	fbR = fbRInterpol + fxCrossover * fbLInterpol;

    	// --- feedback

    	float fdbck = clamp(fxFeedback , -feedbackAttn, feedbackAttn );

    	feedbackBuffer[ feedbackWritePos ] 		=
    			clamp( (lpL * fxInputLevel +	fwL + fbL		* fdbck ) * feedbackAttn  , -1, 1);

    	feedbackBuffer[ feedbackWritePos + 1 ] 	=
    			clamp( (lpR * fxInputLevel +	fwR + fbR	 	* fdbck ) * feedbackAttn  , -1, 1);

    	// --- final allpass

    	_ly1L = coef4L * (_ly1L + feedbackBuffer[  feedbackWritePos  ]) - _lx1L; // allpass
    	_lx1L = feedbackBuffer[  feedbackWritePos  ];

    	_ly1R = coef4R * (_ly1R + feedbackBuffer[ feedbackWritePos+1 ]) - _lx1R; // allpass
    	_lx1R = feedbackBuffer[ feedbackWritePos+1 ];

    	_ly1L = coef4L * (_ly1L + feedbackBuffer[  feedbackWritePos  ]) - _lx1L; // allpass
    	_lx1L = feedbackBuffer[  feedbackWritePos  ];

    	_ly1R = coef4R * (_ly1R + feedbackBuffer[ feedbackWritePos+1 ]) - _lx1R; // allpass
    	_lx1R = feedbackBuffer[ feedbackWritePos+1 ];

    	// --- tremolo

    	vcaL = _ly1L ;//* ((tremoloModL * tremoloPanDepth + (1-tremoloPanDepth)) * 2 - 1);
    	vcaR = _ly1R ;//* ((tremoloModR * tremoloPanDepth + (1-tremoloPanDepth)) * 2 - 1);

    	// --- mix out

    	*(outBuff++) += (int32_t) ( vcaL * sampleMultipler);
    	*(outBuff++) += (int32_t) ( vcaR * sampleMultipler);

    	sample += 2;

    	// --- write index increment

    	forwardWritePos 	+= 2;
    	if( forwardWritePos >= forwardBufferSize )
    		forwardWritePos = 0;

    	feedbackWritePos 	+= 2;
    	if( feedbackWritePos >= feedbackBufferSize )
    		feedbackWritePos = 0;

    	// --- lfo increment

    	lfo1 += lfo1Inc * clamp(fxSpeed + envMod, 0, 1);
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


    	lfoTremolo += lfoTremoloInc * clamp(fxTremoloSpeed + tremoloEnvFollowMod, 0, 1);
    	lfoTremoloSin = tri2sin(lfoTremolo);

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
	float y2 = forwardBuffer[(readPos + 2) ];
	float y3 = forwardBuffer[(readPos + 4) ];

	float x = y2 - y1;
	x -= floor(x);

    return hermite4(x, y0, y1, y2, y3);
}
float FxBus::feedbackHermiteInterpolation(int readPos) {

	float y0 = feedbackBuffer[(readPos + feedbackBufferSize - 2) % feedbackBufferSize];
	float y1 = feedbackBuffer[readPos];
	float y2 = feedbackBuffer[(readPos + 2) ];
	float y3 = feedbackBuffer[(readPos + 4) ];

	float x = y2 - y1;
	x -= floor(x);

    return hermite4(x, y0, y1, y2, y3);
}
void FxBus::vcf1(int readPos) {

	float inmix, hpL, hpR;
    const float f = feedbackLp + 0.05f;

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

    v2L += feedbackLp * v3L;						// lowpass
    v3L += feedbackLp * ( (inmix) - v2L - v3L);

    v2L += feedbackLp * v3L;						// lowpass
    v3L += feedbackLp * ( (inmix) - v2L - v3L);

    feedbackBuffer[readPos] = (v2L);
}
void FxBus::vcf2R(int readPos) {

	// Right voice
    float inmix = feedbackBuffer[readPos];

    v2R += feedbackLp * v3R;						// lowpass
    v3R += feedbackLp * ((inmix) - v2R - v3R);

    v2R += feedbackLp * v3R;						// lowpass
    v3R += feedbackLp * ((inmix) - v2R - v3R);


    feedbackBuffer[readPos] = (v2R);
}
