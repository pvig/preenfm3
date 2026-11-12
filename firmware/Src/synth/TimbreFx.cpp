#include "Timbre.h"
#include "Voice.h"

#if defined(__GNUC__) && !defined(DEBUG)
#pragma GCC push_options
#pragma GCC optimize ("Ofast", "fast-math")
#endif

#define INV_BLOCK_SIZE (1.0f / BLOCK_SIZE)

extern float noise[32];

inline float modulo2(float readPos, int bufferLen)
{
    return unlikely(readPos < 0) ? readPos + bufferLen : readPos;
}
inline float modulo(float d, float max)
{
    return unlikely(d >= max) ? d - max : d;
}
inline float clamp(float d, float min, float max)
{
    const float t = unlikely(d < min) ? min : d;
    return unlikely(t > max) ? max : t;
}
inline float sqrt3(const float x)
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
inline float foldAbs(float x)
{
    x *= 0.5f;
    float f = (fabsf(x - roundf(x)));
    return f + f;
}
inline float fold(float x4)
{
    return fabsf(x4 + 0.25f - roundf(x4 + 0.25f)) - 0.25f;
}
inline float sigmoid(float x)
{
    return x * (1.5f - 0.5f * x * x);
}
inline float sigmoidPos(float x)
{
    // x : 0 -> 1
    return (sigmoid((x * 2) - 1) + 1) * 0.5f;
}
inline float tanh4(float x)
{
    return x / sqrt3(x * x + 1);
}
inline float fastSin(float x)
{
    return 3.9961f * x * (1 - fabsf(x));
}
inline float hann(float x)
{
    float s = sqrt3(x * (1 - x));
    return s + s;
}
inline float max(float x1, float x2)
{
    return (x1 > x2) ? x1 : x2;
}
static inline float fast_log2f(float x)
{
    union
    {
        float f;
        uint32_t i;
    } vx = {x};
    int exp = (int)((vx.i >> 23) & 0xFF) - 127;        // exposant IEEE754
    float mant = (vx.i & 0x7FFFFF) / (float)(1 << 23); // mantisse normalisée
    return exp + mant;                                 // approx log2
}
inline float fast_expf(float x)
{
    union
    {
        float f;
        int32_t i;
    } u;
    u.i = (int32_t)(12102203.0f * x + 127 * (1 << 23));
    return u.f;
}
inline float fast_pow2(float p)
{
    // approx 2^p
    union
    {
        uint32_t i;
        float f;
    } v;
    v.i = (uint32_t)((1 << 23) * (p + 126.94269504f));
    return v.f;
}
inline float fastExpNeg(float x)
{
    // x > 0
    return fast_expf(-x);
}
inline float fast_cos_2pi(float x)
{
    // x entre 0 et 0.5 environ (fréquence normalisée)
    float x2 = x * x;
    return 1.0f - 19.7392088f * x2 + 64.3776844f * x2 * x2;
}
static inline float soft_clip(float x) {
    return x / (1.0f + fabsf(x));
}
static inline float blend_screen_smooth(float x1, float x2) {
    float a = 0.5f * (x1 + 1.0f);
    float b = 0.5f * (x2 + 1.0f);
    float y = a + b - a * b; // smooth polynomial, no branches
    return 2.0f * y - 1.0f;
}
static inline float blend_softlight(float x1, float x2) {
    float a = 0.5f * (x1 + 1.0f);
    float b = 0.5f * (x2 + 1.0f);
    float y = a - (1.0f - 2.0f * b) * a * (1.0f - a); // continuous polynomial
    return 2.0f * y - 1.0f;
}
static inline float mix_soft(float x1, float x2, float drive) {
    float m = 0.5f * (x1 + x2);
    float d = 0.5f * (x1 - x2);
    float y = m + d * (1.0f - drive * fabsf(m));
    return soft_clip(y);
};
// all pass params
const float f1 = 0.0156f;
const float apcoef1 = (1.0f - f1) / (1.0f + f1);
const float f2 = (0.17f + f1);
const float apcoef2 = (1.0f - f2) / (1.0f + f2);
const float f3 = (0.17f + f2);
const float apcoef3 = (1.0f - f3) / (1.0f + f3);
const float f4 = (0.17f + f3);
const float apcoef4 = (1.0f - f4) / (1.0f + f4);

// delay sizes
const float delayBufferSizeF = delayBufferSize;
const float delayBufferSize90 = delayBufferSize * 0.25f;
const float delayBufferSize180 = delayBufferSize * 0.5f;
const int delayBufferSizeM1 = delayBufferSize - 1;
const int delayBufferSizeM4 = delayBufferSize - 4;
const float delayBufferSizeInv = 1.0f / delayBufferSize;
const int delayBufStereoSize = delayBufferSize * 0.5f;
const int delayBufStereoSizeM1 = delayBufStereoSize - 1;
const float delayBufStereoDiv4 = delayBufStereoSize * 0.25f;
const float delayBufStereoSizeInv = 1.0f / delayBufStereoSize;

void Timbre::initFx()
{
    hb_x1[0] = &hb1_x1;
    hb_x1[1] = &hb2_x1;
    hb_x1[2] = &hb3_x1;
    hb_x1[3] = &hb4_x1;
    hb_x1[4] = &hb5_x1;
    hb_x1[5] = &hb6_x1;
    hb_x1[6] = &hb7_x1;
    hb_x1[7] = &hb8_x1;

    hb_x2[0] = &hb1_x2;
    hb_x2[1] = &hb2_x2;
    hb_x2[2] = &hb3_x2;
    hb_x2[3] = &hb4_x2;
    hb_x2[4] = &hb5_x2;
    hb_x2[5] = &hb6_x2;
    hb_x2[6] = &hb7_x2;
    hb_x2[7] = &hb8_x2;

    hb_y1[0] = &hb1_y1;
    hb_y1[1] = &hb2_y1;
    hb_y1[2] = &hb3_y1;
    hb_y1[3] = &hb4_y1;
    hb_y1[4] = &hb5_y1;
    hb_y1[5] = &hb6_y1;
    hb_y1[6] = &hb7_y1;
    hb_y1[7] = &hb8_y1;

    hb_y2[0] = &hb1_y2;
    hb_y2[1] = &hb2_y2;
    hb_y2[2] = &hb3_y2;
    hb_y2[3] = &hb4_y2;
    hb_y2[4] = &hb5_y2;
    hb_y2[5] = &hb6_y2;
    hb_y2[6] = &hb7_y2;
    hb_y2[7] = &hb8_y2;
}

void Timbre::fxAfterBlock()
{
    int fx2Type = params_.effect2.type;
    if (fx2Type == FILTER2_OFF)
    {
        prevFx2Type = fx2Type;
        this->newNotePlayed = false;
        return;
    }

    if (!voices_[lastPlayedNote_]->isPlaying())
    {
        // this voice is not playing but still need to calculate lfo
        voices_[lastPlayedNote_]->matrix.computeAllDestinations();
    }

    float matrixFilterFrequency = voices_[lastPlayedNote_]->matrix.getDestination(FILTER2_PARAM1);
    float matrixFilterParam2 = voices_[lastPlayedNote_]->matrix.getDestination(FILTER2_PARAM2);
    float matrixFilterAmp = voices_[lastPlayedNote_]->matrix.getDestination(FILTER2_AMP);
    float matrixFilterPan = clamp(voices_[lastPlayedNote_]->matrix.getDestination(ALL_PAN), -1, 1);
    float gainTmp = clamp(this->params_.effect2.param3 + matrixFilterAmp, 0, 16);

    if (prevFx2Type != fx2Type)
    {
        // anti click on fx change
        mixerGain_ = 0;
        feedbackInput = 0;
        feedback = 0;
        for (int s = 0; s < delayBufferSize; s++)
        {
            delayBuffer_[s] = 0;
        }
        low1 = low2 = low3 = low4 = low5 = low6 = 0;
        band1 = band2 = band3 = band4 = band5 = band6 = 0;
        hb1_x1 = hb1_x2 = hb1_y1 = hb1_y2 = 0;
        hb2_x1 = hb2_x2 = hb2_y1 = hb2_y2 = 0;
        hb3_x1 = hb3_x2 = hb3_y1 = hb3_y2 = 0;
        hb4_x1 = hb4_x2 = hb4_y1 = hb4_y2 = 0;
        hb5_x1 = hb5_x2 = hb5_y1 = hb5_y2 = 0;
        hb6_x1 = hb6_x2 = hb6_y1 = hb6_y2 = 0;
        hb7_x1 = hb7_x2 = hb7_y1 = hb7_y2 = 0;
        hb8_x1 = hb8_x2 = hb8_y1 = hb8_y2 = 0;
        inputWritePos1 = inputWritePos2 = inputWritePos3 = inputWritePos4 = inputWritePos5 = 0;
    }
    prevFx2Type = fx2Type;

    switch (fx2Type)
    {
    case FILTER2_FLANGE:
    {
        mixerGain_ = 0.02f * gainTmp + .98f * mixerGain_;
        float mixerGain_01 = clamp(mixerGain_, 0, 1);
        int mixerGain255 = mixerGain_01 * 255;
        float dry = panTable[255 - mixerGain255];
        float wet = panTable[mixerGain255] * 0.5f;
        float extraAmp = clamp(mixerGain_ - 1, 0, 1);
        wet += extraAmp;

        param1S = 0.02f * this->params_.effect2.param1 + .98f * param1S;
        float speed = clamp(param1S * param1S * param1S, 0, 1);

        float lfo1Inc = clamp(hb8_y1 * speed * 0.015f, -1, 1);
        hb8_x1 += lfo1Inc;
        if (hb8_x1 >= 1)
        {
            hb8_x1 = 1;
            hb8_y1 = -1;
        }
        if (hb8_x1 <= 0)
        {
            hb8_x1 = 0;
            hb8_y1 = 1;
        }

        float lfoRaw = foldAbs(hb8_x1 + matrixFilterFrequency) * 0.5f;
        float fxParamTmp = sigmoidPos(lfoRaw); 
        delayReadFrac = (fxParamTmp + 99 * delayReadFrac) * 0.01f; // smooth change

        float currentDelaySize1 = clamp(delaySize1, 2.0f, (float)(delayBufStereoSize - 3));
        delaySize1 = clamp(delayBufStereoDiv4 * delayReadFrac, 2.0f, (float)(delayBufStereoSize - 3));
        float delaySizeInc1 = (delaySize1 - currentDelaySize1) * INV_BLOCK_SIZE;

        float currentFeedback = feedback;
        feedback = clamp(sqrt3(fabsf(this->params_.effect2.param2 + matrixFilterParam2)), -1, 1) * 0.95f;
        float feedbackInc = (feedback - currentFeedback) * INV_BLOCK_SIZE;

        float * __restrict__ sp = sampleBlock_;
        float * __restrict__ delayBuf = delayBuffer_;

        float delayReadPos90;

        float filterB2 = 0.1f;
        float filterB = (filterB2 * filterB2 * 0.5f);

        float _in3_b1 = (1 - filterB);
        float _in3_a0 = (1 + _in3_b1 * _in3_b1 * _in3_b1) * 0.5f;

        const float f = 0.8f;
        const float f2 = 0.85f;
        const float fnotch = 1.03f;
        const float driftOffset = 37.f + fabsf(noise[6]) * 0.8f;

        int dWritePos = delayWritePos;
        float dReadPos = delayReadPos;
        float l1 = low1, b1 = band1, l2 = low2, b2 = band2;
        float l3 = low3, b3 = band3, l4 = low4, b4 = band4;
        float l5 = low5, l6 = low6;
        float ly1 = _ly1, lx1 = _lx1, ly2 = _ly2, lx2 = _lx2;
        float h4_x1 = hb4_x1, h4_y1 = hb4_y1, h4_x2 = hb4_x2, h4_y2 = hb4_y2;
        float h1_x1 = hb1_x1, h1_y1 = hb1_y1, h1_x2 = hb1_x2, h1_y2 = hb1_y2;
        float hp1_x1 = hp_in_x1, hp1_y1 = hp_in_y1, hp1_y0 = hp_in_y0;
        float hp2_x1 = hp_in2_x1, hp2_y1 = hp_in2_y1, hp2_y0 = hp_in2_y0;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {
            l3 += f * b3;
            b3 += f * (*sp - l3 - b3);

            l4 += f * b4;
            b4 += f * (*(sp + 1) - l4 - b4);

            // feedback
            float feedL = l5 * currentFeedback;
            float feedR = l6 * currentFeedback;

            // L
            ly1 = apcoef2 * (ly1 + feedL) - lx1; // allpass
            lx1 = feedL;

            h4_y1 = apcoef4 * (h4_y1 + ly1) - h4_x1; // allpass
            h4_x1 = ly1;

            // R
            ly2 = apcoef2 * (ly2 + feedR) - lx2; // allpass
            lx2 = feedR;

            h4_y2 = apcoef4 * (h4_y2 + ly2) - h4_x2; // allpass
            h4_x2 = ly2;

            // audio in hp
            float hp_in_x0 = tanh4(l3 + l3 - h4_y1);
            hp1_y0 = _in3_a0 * (hp_in_x0 - hp1_x1) + _in3_b1 * hp1_y1;
            hp1_y1 = hp1_y0;
            hp1_x1 = hp_in_x0;

            float hp_in2_x0 = tanh4(l4 + l4 - h4_y2);
            hp2_y0 = _in3_a0 * (hp_in2_x0 - hp2_x1) + _in3_b1 * hp2_y1;
            hp2_y1 = hp2_y0;
            hp2_x1 = hp_in2_x0;

            dWritePos = (dWritePos + 1) & delayBufStereoSizeM1;
            delayBuf[dWritePos] = hp1_y0;
            delayBuf[dWritePos + delayBufStereoSize] = hp2_y0;

            dReadPos = modulo2((float)dWritePos - currentDelaySize1, delayBufStereoSize);
            delayReadPos90 = modulo2(dReadPos - driftOffset, delayBufStereoSize);

            l5 = delayInterpolationHermite(dReadPos, delayBuf, delayBufStereoSizeM1, 0);
            l6 = delayInterpolationHermite(delayReadPos90, delayBuf, delayBufStereoSizeM1, delayBufStereoSize);

            l1 += f2 * b1;
            b1 += f2 * (l5 - l1 - b1);

            l2 += f2 * b2;
            b2 += f2 * (l6 - l2 - b2);

            // notch L
            h1_x1 += fnotch * h1_y1;
            float high7 = l1 - h1_x1 - h1_y1;
            h1_y1 += fnotch * high7;
            float notchL = (high7 + h1_x1);

            // notch R
            h1_x2 += fnotch * h1_y2;
            float high8 = l2 - h1_x2 - h1_y2;
            h1_y2 += fnotch * high8;
            float notchR = (high8 + h1_x2);

            *sp = *sp * dry + notchL * wet;
            sp++;
            *sp = *sp * dry + notchR * wet;
            sp++;

            currentDelaySize1 += delaySizeInc1;
            currentFeedback += feedbackInc;
        }

        delayWritePos = dWritePos;
        delayReadPos = dReadPos;
        low1 = l1; band1 = b1; low2 = l2; band2 = b2;
        low3 = l3; band3 = b3; low4 = l4; band4 = b4;
        low5 = l5; low6 = l6;
        _ly1 = ly1; _lx1 = lx1; _ly2 = ly2; _lx2 = lx2;
        hb4_x1 = h4_x1; hb4_y1 = h4_y1; hb4_x2 = h4_x2; hb4_y2 = h4_y2;
        hb1_x1 = h1_x1; hb1_y1 = h1_y1; hb1_x2 = h1_x2; hb1_y2 = h1_y2;
        hp_in_x1 = hp1_x1; hp_in_y1 = hp1_y1; hp_in_y0 = hp1_y0;
        hp_in2_x1 = hp2_x1; hp_in2_y1 = hp2_y1; hp_in2_y0 = hp2_y0;
    }
    break;
    case FILTER2_CHORUS:
    {
        mixerGain_ = 0.02f * gainTmp + .98f * mixerGain_;
        float mixerGain_01 = clamp(mixerGain_, 0, 1);
        int mixerGain255 = mixerGain_01 * 255;
        float dry = panTable[255 - mixerGain255];
        float wet = panTable[mixerGain255] * 0.5f;
        float extraAmp = clamp(mixerGain_ - 1, 0, 1);
        wet += extraAmp;

        float wetL = wet * (1 + matrixFilterPan);
        float wetR = wet * (1 - matrixFilterPan);

        param1S = 0.02f * this->params_.effect2.param1 + .98f * param1S;
        float speed = clamp(param1S * param1S + matrixFilterFrequency, 0, 1);
        float speedMod = speed * 4;

        param2S = 0.01f * (this->params_.effect2.param2) + .99f * param2S;
        float width = param2S * 0.6f;

        float lfo1Inc = clamp(hb8_y1 * 0.0007f * speedMod, -1.f, 1.f);
        hb8_x1 += lfo1Inc;
        if (hb8_x1 >= 1)
        {
            hb8_x1 = 1;
            hb8_y1 = -1;
        }
        if (hb8_x1 <= 0)
        {
            hb8_x1 = 0;
            hb8_y1 = 1;
        }

        float lfo2Inc = clamp(hb8_y2 * 0.008f * speedMod, -1.f, 1.f);
        hb8_x2 += lfo2Inc;
        if (hb8_x2 >= 1)
        {
            hb8_x2 = 1;
            hb8_y2 = -1;
        }
        if (hb8_x2 <= 0)
        {
            hb8_x2 = 0;
            hb8_y2 = 1;
        }

        float lfo = (hb8_x1 * 0.3f + sigmoidPos(hb8_x2) * width) * 0.5f;

        float fxParamTmp = (lfo + matrixFilterParam2);
        delayReadFrac = (fxParamTmp + 99 * delayReadFrac) * 0.01f; // smooth change

        float readPos1 = foldAbs(delayReadFrac);
        float readPos2 = foldAbs(readPos1 + 0.3333f);
        float readPos3 = foldAbs(readPos1 + 0.6666f);

        float currentDelaySize1 = clamp(delaySize1, 2.0f, (float)(delayBufferSize - 3));
        float currentDelaySize2 = clamp(delaySize2, 2.0f, (float)(delayBufferSize - 3));
        float currentDelaySize3 = clamp(delaySize3, 2.0f, (float)(delayBufferSize - 3));
        float chorusSize = delayBufStereoSize; // 21ms
        delaySize1 = clamp(1 + chorusSize * readPos1, 2.0f, (float)(delayBufferSize - 3));
        delaySize2 = clamp(1 + chorusSize * readPos2, 2.0f, (float)(delayBufferSize - 3));
        delaySize3 = clamp(1 + chorusSize * readPos3, 2.0f, (float)(delayBufferSize - 3));
        float delaySizeInc1 = (delaySize1 - currentDelaySize1) * INV_BLOCK_SIZE;
        float delaySizeInc2 = (delaySize2 - currentDelaySize2) * INV_BLOCK_SIZE;
        float delaySizeInc3 = (delaySize3 - currentDelaySize3) * INV_BLOCK_SIZE;

        float * __restrict__ sp = sampleBlock_;
        float * __restrict__ delayBuf = delayBuffer_;

        float delReadPos1, delReadPos2, delReadPos3, monoIn;
        float delayOut1 = 0, delayOut2 = 0, delayOut3 = 0;

        const float f = 0.8f;
        const float f2 = 0.87f;

        // hi pass params
        float filterB2 = 0.25f;
        float filterB = (filterB2 * filterB2 * 0.5f);

        float _in_b1 = (1 - filterB);
        float _in_a0 = (1 + _in_b1 * _in_b1 * _in_b1) * 0.5f;

        int dWritePos = delayWritePos;
        float l1 = low1, b1 = band1, l2 = low2, b2 = band2;
        float l3 = low3, b3 = band3, l5 = low5, b5 = band5, l6 = low6, b6 = band6;
        float hp1_x1 = hp_in_x1, hp1_y1 = hp_in_y1, hp1_y0 = hp_in_y0;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {
            monoIn = (*sp + *(sp + 1)) * 0.5f;

            // input lp
            l3 += f * b3;
            b3 += f * (monoIn - l3 - b3);

            // audio in hp
            float hp_in_x0 = l3;
            hp1_y0 = _in_a0 * (hp_in_x0 - hp1_x1) + _in_b1 * hp1_y1;
            hp1_y1 = hp1_y0;
            hp1_x1 = hp_in_x0;

            float bass = l3 - hp1_y0;

            dWritePos = (dWritePos + 1) & delayBufferSizeM1;
            delayBuf[dWritePos] = hp1_y0;

            float dWritePosF = (float)dWritePos;

            delReadPos1 = modulo2(dWritePosF - currentDelaySize1, delayBufferSize);
            delayOut1 = delayInterpolationHermite(delReadPos1, delayBuf, delayBufferSizeM1, 0);

            delReadPos2 = modulo2(dWritePosF - currentDelaySize2, delayBufferSize);
            delayOut2 = delayInterpolationHermite(delReadPos2, delayBuf, delayBufferSizeM1, 0);

            delReadPos3 = modulo2(dWritePosF - currentDelaySize3, delayBufferSize);
            delayOut3 = delayInterpolationHermite(delReadPos3, delayBuf, delayBufferSizeM1, 0);

            float delaySumOut = bass + (delayOut1 - delayOut3 + delayOut2);
            float delaySumOut2 = bass + (delayOut3 - delayOut1 + delayOut2);

            l1 += f2 * b1;
            b1 += f2 * (delaySumOut - l1 - b1);

            l2 += f2 * b2;
            b2 += f2 * (l1 - l2 - b2);

            l5 += f2 * b5;
            b5 += f2 * (delaySumOut2 - l5 - b5);

            l6 += f2 * b6;
            b6 += f2 * (l5 - l6 - b6);

            *sp = *sp * dry + l2 * wetL;
            sp++;
            *sp = *sp * dry + l6 * wetR;
            sp++;

            currentDelaySize1 += delaySizeInc1;
            currentDelaySize2 += delaySizeInc2;
            currentDelaySize3 += delaySizeInc3;
        }

        delayWritePos = dWritePos;
        delayWritePosF = (float)dWritePos;
        low1 = l1; band1 = b1; low2 = l2; band2 = b2;
        low3 = l3; band3 = b3; low5 = l5; band5 = b5; low6 = l6; band6 = b6;
        hp_in_x1 = hp1_x1; hp_in_y1 = hp1_y1; hp_in_y0 = hp1_y0;
    }
    break;
    case FILTER2_DIMENSION:
    {
        mixerGain_ = 0.02f * gainTmp + .98f * mixerGain_;
        float mixerGain_01 = clamp(mixerGain_, 0, 1);
        int mixerGain255 = mixerGain_01 * 255;
        float dry = panTable[255 - mixerGain255];
        float wet = panTable[mixerGain255];
        float extraAmp = clamp(mixerGain_ - 1, 0, 1);
        wet += extraAmp;

        float param1 = this->params_.effect2.param1;
        float speed = clamp(param1 * param1 * 0.5f, 0, 1);

        float lfo1Inc = clamp(hb8_y1 * speed * 0.01f, -1, 1);
        hb8_x1 += lfo1Inc;
        if (hb8_x1 >= 1)
        {
            hb8_x1 = 1;
            hb8_y1 = -1;
        }
        if (hb8_x1 <= 0)
        {
            hb8_x1 = 0;
            hb8_y1 = 1;
        }

        float lfo2Inc = clamp(hb8_y2 * speed * 0.0097f, -1, 1); // ~3% plus lent
        hb8_x2 += lfo2Inc;
        if (hb8_x2 >= 1) { hb8_x2 = 1; hb8_y2 = -1; }
        if (hb8_x2 <= 0) { hb8_x2 = 0; hb8_y2 = 1; }
        
        float lfo = hb8_x1 * 0.5f;
        float lfo2 = hb8_x2 * 0.5f;

        param1S = 0.02f * matrixFilterFrequency + .98f * param1S;

        matrixFilterFrequency *= 0.5f;

        float fxParamTmp = sigmoidPos(foldAbs(0.25f + (0.125f * ((lfo + param1S)))));
        delayReadFrac = (fxParamTmp + 99 * delayReadFrac) * 0.01f; // smooth change       
        float fxParamTmp2 = sigmoidPos(foldAbs(0.25f + (0.125f * (lfo2 + param1S))));
        delayReadFrac2 = (fxParamTmp2 + 99 * delayReadFrac2) * 0.01f; // smooth change

        float currentDelaySize1 = clamp(delaySize1, 2.0f, (float)(delayBufStereoSize - 3));
        float currentDelaySize2 = clamp(delaySize2, 2.0f, (float)(delayBufStereoSize - 3));
        delaySize1 = clamp(1 + delayBufStereoSize * delayReadFrac, 2.0f, (float)(delayBufStereoSize - 3));
        delaySize2 = clamp(1 + delayBufStereoSize * delayReadFrac2, 2.0f, (float)(delayBufStereoSize - 3));
        float delaySizeInc1 = (delaySize1 - currentDelaySize1) * INV_BLOCK_SIZE;
        float delaySizeInc2 = (delaySize2 - currentDelaySize2) * INV_BLOCK_SIZE;

        param2S = 0.1f * (this->params_.effect2.param2 + matrixFilterParam2) + .9f * param2S;

        const float f = 0.82f;
        const float f2 = 0.2f;

        // hi pass params
        float filterB = (f2 * f2 * 0.5f);

        float _in2_b1 = (1 - filterB);
        float _in2_a0 = (1 + _in2_b1 * _in2_b1 * _in2_b1) * 0.5f;

        float * __restrict__ sp = sampleBlock_;
        float * __restrict__ delayBuf = delayBuffer_;

        float delReadPos, delReadPos2;
        float delayOut1, delayOut3 = 0;

        // mid / side
        float outL, outR, mid, side;
        float width = clamp(param2S * 2, 0, 2);
        float tmp = 1.f / (1.f + width);
        float coef_M = tmp * wet;
        float coef_S = width * tmp * wet;

        int dWritePos = delayWritePos;
        float l1 = low1, b1 = band1, l2 = low2, b2 = band2;
        float hp1_x1 = hp_in_x1, hp1_y1 = hp_in_y1, hp1_y0 = hp_in_y0;
        float hp2_x1 = hp_in2_x1, hp2_y1 = hp_in2_y1, hp2_y0 = hp_in2_y0;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {
            l1 += f * b1;
            b1 += f * (*sp - l1 - b1);

            l2 += f * b2;
            b2 += f * (*(sp + 1) - l2 - b2);

            // audio in hp
            float hp_in_x0 = l1;
            hp1_y0 = _in2_a0 * (hp_in_x0 - hp1_x1) + _in2_b1 * hp1_y1;
            hp1_y1 = hp1_y0;
            hp1_x1 = hp_in_x0;

            float hp_in2_x0 = l2;
            hp2_y0 = _in2_a0 * (hp_in2_x0 - hp2_x1) + _in2_b1 * hp2_y1;
            hp2_y1 = hp2_y0;
            hp2_x1 = hp_in2_x0;

            float lpc1 = (l1 - hp1_y0);
            float lpc2 = (l2 - hp2_y0);

            dWritePos = (dWritePos + 1) & delayBufStereoSizeM1;

            delayBuf[dWritePos] = hp1_y0 + hp2_y0 * 0.05f; // cross feed for more dimension
            delayBuf[dWritePos + delayBufStereoSize] = hp2_y0 + hp1_y0 * 0.05f;

            float dWritePosF = (float)dWritePos;

            delReadPos = modulo2(dWritePosF - currentDelaySize1, delayBufStereoSize);
            delReadPos2 = modulo2(dWritePosF - currentDelaySize2, delayBufStereoSize);

            delayOut1 = delayInterpolationHermite(delReadPos, delayBuf, delayBufStereoSizeM1, 0);
            delayOut3 = delayInterpolationHermite(delReadPos2, delayBuf, delayBufStereoSizeM1, delayBufStereoSize);

            outL = lpc1 + delayOut3;
            outR = lpc2 + delayOut1;

            mid = coef_M * (outL + outR);
            side = coef_S * (outL - outR);

            *sp = (*sp) * dry + (mid + side);
            sp++;
            *sp = (*sp) * dry + (mid - side);
            sp++;

            currentDelaySize1 += delaySizeInc1;
            currentDelaySize2 += delaySizeInc2;
        }

        delayWritePos = dWritePos;
        delayWritePosF = (float)dWritePos;
        low1 = l1; band1 = b1; low2 = l2; band2 = b2;
        hp_in_x1 = hp1_x1; hp_in_y1 = hp1_y1; hp_in_y0 = hp1_y0;
        hp_in2_x1 = hp2_x1; hp_in2_y1 = hp2_y1; hp_in2_y0 = hp2_y0;
    }
    break;
    case FILTER2_DOUBLER:
    {
        mixerGain_ = 0.02f * gainTmp + .98f * mixerGain_;
        float mixerGain_01 = clamp(mixerGain_, 0, 1);
        int mixerGain255 = mixerGain_01 * 255;
        float dry = panTable[255 - mixerGain255];
        float wet = panTable[mixerGain255] * 0.75f;
        float extraAmp = clamp(mixerGain_ - 1, 0, 1);
        wet += extraAmp;

        float wetL = wet * (1 + matrixFilterPan);
        float wetR = wet * (1 - matrixFilterPan);

        param1S = 0.1f * this->params_.effect2.param1 + .9f * param1S;
        param2S = 0.1f * (this->params_.effect2.param2 + matrixFilterParam2) + .9f * param2S;

        float currentShift = shift;
        shift = clamp(fabsf(param1S * 2 + matrixFilterFrequency * 0.5f), 0, 16);
        float shiftInc = (shift - currentShift) * INV_BLOCK_SIZE;

        float quadrant = fabsf((clamp(shift, 0, 2) - 1)); // 2 quadrant for up & down shift
        float feedbackZeroZone = clamp(0.82f + quadrant * 40, 0, 0.999f);

        float feed = clamp(param2S * feedbackZeroZone, 0, 1);
        float currentFeedback = feedback;
        feedback = clamp(feed, -1, 1) * 0.44f;
        float feedbackInc = (feedback - currentFeedback) * INV_BLOCK_SIZE;

        float lpZeroZone = clamp(quadrant * 50, 0, 1);

        float filterA2 = 0.5f + lpZeroZone * 0.2f;
        float filterA = (filterA2 * filterA2 * 0.5f);
        float _in_lp_b = 1 - filterA;
        float _in_lp_a = 1 - _in_lp_b;

        const float f = 0.7f;
        const float f2 = 0.38f;
        const float f3 = 0.72f;

        float * __restrict__ sp = sampleBlock_;
        float * __restrict__ delayBuf = delayBuffer_;

        float delayReadPos180, level1, level2;
        float delayOut1, delayOut2;

        // hi pass params
        float hpZeroZone = 1 - clamp(sqrt3(quadrant) * 4.5f, 0, 1);

        float filterB2 = 0.1f + (0.32f - (param1S * 0.2f) + hpZeroZone * 0.3f) * clamp(param2S * 1.2f, 0, 1);
        float filterB = (filterB2 * filterB2 * 0.5f);

        float _in2_b1 = (1 - filterB);
        float _in2_a0 = (1 + _in2_b1 * _in2_b1 * _in2_b1) * 0.5f;

        int dWritePos = delayWritePos;
        float dReadPos = delayReadPos;
        float fbInput = feedbackInput;
        float h1_x1 = hb1_x1, h1_x2 = hb1_x2;
        float hp1_x1 = hp_in_x1, hp1_y1 = hp_in_y1, hp1_y0 = hp_in_y0;
        float l1 = low1, b1 = band1, l2 = low2, b2 = band2;
        float l3 = low3, b3 = band3, l4 = low4, b4 = band4;
        float l5 = low5, b5 = band5, l6 = low6, b6 = band6;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {
            float monoIn = (*sp + *(sp + 1)) * 0.5f;

            // feedback lp
            h1_x1 = _in_lp_a * fbInput + h1_x1 * _in_lp_b;
            h1_x2 = _in_lp_b * h1_x1 + h1_x2 * _in_lp_b;

            // delay in hp
            float hp_in_x0 = clamp(tanh4((monoIn + h1_x2 * currentFeedback)), -1, 1);
            hp1_y0 = _in2_a0 * (hp_in_x0 - hp1_x1) + _in2_b1 * hp1_y1;
            hp1_y1 = hp1_y0;
            hp1_x1 = hp_in_x0;

            dWritePos = (dWritePos + 1) & delayBufferSizeM1;

            l1 += f * b1;
            b1 += f * (hp1_y0 - l1 - b1);
            l2 += f * b2;
            b2 += f * (l1 - l2 - b2);

            delayBuf[dWritePos] = l2;

            dReadPos = modulo(dReadPos + currentShift, delayBufferSize);
            delayReadPos180 = modulo2(dReadPos - delayBufferSize180, delayBufferSize);

            delayOut1 = delayInterpolation(dReadPos, delayBuf, delayBufferSizeM1);
            delayOut2 = delayInterpolation(delayReadPos180, delayBuf, delayBufferSizeM1);

            float dWritePosF = (float)dWritePos;
            float rwp1 = modulo2(dWritePosF - dReadPos, delayBufferSize);
            float rwp3 = modulo2(dWritePosF - delayReadPos180, delayBufferSize);

            level1 = hann(rwp1 * delayBufferSizeInv);
            level2 = hann(rwp3 * delayBufferSizeInv);

            float out1 = delayOut1 * (level1);
            float out2 = delayOut2 * (level2);

            float delaySumOut = out1 + out2;
            float delaySumOut2 = out1 - out2;

            fbInput = delaySumOut;

            l5 += f3 * b5;
            b5 += f3 * (delaySumOut - l5 - b5);

            l6 += f3 * b6;
            b6 += f3 * (delaySumOut2 - l6 - b6);

            // bass boost
            l3 += f2 * b3;
            b3 += f2 * (l5 - l3 - b3);
            l4 += f2 * b4;
            b4 += f2 * (l6 - l4 - b4);

            *sp = (*sp * dry + ((l3 + l5) * wetL));
            sp++;
            *sp = (*sp * dry + ((l4 + l6) * wetR));
            sp++;

            currentFeedback += feedbackInc;
            currentShift += shiftInc;
        }

        delayWritePos = dWritePos;
        delayWritePosF = (float)dWritePos;
        delayReadPos = dReadPos;
        feedbackInput = fbInput;
        hb1_x1 = h1_x1; hb1_x2 = h1_x2;
        hp_in_x1 = hp1_x1; hp_in_y1 = hp1_y1; hp_in_y0 = hp1_y0;
        low1 = l1; band1 = b1; low2 = l2; band2 = b2;
        low3 = l3; band3 = b3; low4 = l4; band4 = b4;
        low5 = l5; band5 = b5; low6 = l6; band6 = b6;
    }
    break;
    case FILTER2_TRIPLER:
    {
        mixerGain_ = 0.02f * gainTmp + .98f * mixerGain_;
        float mixerGain_01 = clamp(mixerGain_, 0, 1);
        int mixerGain255 = mixerGain_01 * 255;
        float dry = panTable[255 - mixerGain255];
        float wet = panTable[mixerGain255] * 0.75f;
        float extraAmp = clamp(mixerGain_ - 1, 0, 1);
        wet += extraAmp;

        float wetL = wet * (1 + matrixFilterPan);
        float wetR = wet * (1 - matrixFilterPan);

        param1S = 0.005f * this->params_.effect2.param1 + .995f * param1S;
        param2S = 0.005f * this->params_.effect2.param2 + .995f * param2S;

        float currentShift = shift;
        shift = clamp(fabsf(param1S * 2 + matrixFilterFrequency * 0.5f), 0, 16);
        float shiftInc = (shift - currentShift) * INV_BLOCK_SIZE;

        float currentShift2 = shift2;
        shift2 = clamp(fabsf(param2S * 2 + matrixFilterParam2 * 0.5f), 0, 16);
        float shiftInc2 = (shift2 - currentShift2) * INV_BLOCK_SIZE;

        const float f = 0.75f;
        const float f2 = 0.72f;

        float * __restrict__ sp = sampleBlock_;
        float * __restrict__ delayBuf = delayBuffer_;

        float level1, level2, level3, level4;
        float delayOut1 = 0, delayOut2 = 0, delayOut3 = 0, delayOut4 = 0;
        float delayReadPos180;

        int dWritePos = delayWritePos;
        float dReadPos = delayReadPos;
        float dReadPos2 = delayReadPos2;
        float l1 = low1, b1 = band1, l3 = low3, b3 = band3;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {
            float monoIn = (*sp + *(sp + 1)) * 0.5f;

            dWritePos = (dWritePos + 1) & delayBufferSizeM1;

            l1 += f * b1;
            b1 += f * (monoIn - l1 - b1);

            delayBuf[dWritePos] = l1;

            //--------------- shifter 1

            dReadPos = modulo(dReadPos + currentShift, delayBufferSize);
            delayReadPos180 = modulo(dReadPos + delayBufferSize180, delayBufferSize);

            delayOut1 = delayInterpolation(dReadPos, delayBuf, delayBufferSizeM1);
            delayOut2 = delayInterpolation(delayReadPos180, delayBuf, delayBufferSizeM1);

            float dWritePosF = (float)dWritePos;
            float rwp1 = modulo2(dWritePosF - dReadPos, delayBufferSize);
            float rwp2 = modulo2(dWritePosF - delayReadPos180, delayBufferSize);

            level1 = hann(rwp1 * delayBufferSizeInv);
            level2 = hann(rwp2 * delayBufferSizeInv);

            float out1 = delayOut1 * (level1);
            float out2 = delayOut2 * (level2);

            float delaySumOut = out1 + out2;

            //--------------- shifter 2

            dReadPos2 = modulo(dReadPos2 + currentShift2, delayBufferSize);
            delayReadPos180 = modulo(dReadPos2 + delayBufferSize90, delayBufferSize);

            delayOut3 = delayInterpolation(dReadPos2, delayBuf, delayBufferSizeM1);
            delayOut4 = delayInterpolation(delayReadPos180, delayBuf, delayBufferSizeM1);

            float rwp3 = modulo2(dWritePosF - dReadPos2, delayBufferSize);
            float rwp4 = modulo2(dWritePosF - delayReadPos180, delayBufferSize);

            level3 = hann(rwp3 * delayBufferSizeInv);
            level4 = hann(rwp4 * delayBufferSizeInv);

            float out3 = delayOut3 * (level3);
            float out4 = delayOut4 * (level4);

            delaySumOut -= out3 + out4;

            // lp output
            l3 += f2 * b3;
            b3 += f2 * ((delaySumOut)-l3 - b3);

            *sp = *sp * dry + l3 * wetL;
            sp++;
            *sp = *sp * dry + l3 * wetR;
            sp++;

            currentShift += shiftInc;
            currentShift2 += shiftInc2;
        }

        delayWritePos = dWritePos;
        delayReadPos = dReadPos;
        delayReadPos2 = dReadPos2;
        low1 = l1; band1 = b1; low3 = l3; band3 = b3;
    }
    break;
    case FILTER2_BODE:
    {
        // dry wet
        mixerGain_ = 0.02f * gainTmp + .98f * mixerGain_;
        float mixerGain_01 = clamp(mixerGain_, 0, 1);
        int mixerGain255 = mixerGain_01 * 255;
        float dry = panTable[255 - mixerGain255];
        float wet = panTable[mixerGain255] * 0.5f;
        float extraAmp = clamp(mixerGain_ - 1, 0, 1);
        wet += extraAmp;

        float wetL = wet * (1 + matrixFilterPan);
        float wetR = wet * (1 - matrixFilterPan);

        float param1 = this->params_.effect2.param1;

        // mix between freq shift + and - :
        int param255 = 0;
        if (param1 <= 0.4f)
        {
            param255 = 0;
        }
        else if (param1 >= 0.6f)
        {
            param255 = 255;
        }
        else
        {
            param255 = 255 * (5.f * (param1 - 0.4f));
        }
        float shiftMinus = panTable[255 - param255];
        float shiftPlus = panTable[param255];

        // shift val
        param1S = 0.1f * param1 + 0.9f * param1S;
        float quadrant = fabsf(param1S - 0.5f); // 2 quadrant for up & down shift
        float quadrantSq = sqrt3(quadrant);

        // shift increment :
        shift = clamp(shift, 0, 0.1f);
        float currentShift = shift;
        float shiftval = clamp(fabsf(quadrant * 1.34f + matrixFilterFrequency * 0.5f), 0, 0.9999f);
        shiftval *= shiftval * 0.05f;
        shift = shift * 0.96f + 0.04f * shiftval;
        float shiftInc = (shift - currentShift) * INV_BLOCK_SIZE;

        // feedback
        float feedbackZeroZone = clamp(0.6f + (quadrant * 6), 0, 1);

        float feedbackParam = clamp(this->params_.effect2.param2 + matrixFilterParam2, 0, 1) * 0.8f;

        feedback = clamp(feedback, 0, 1);

        float currentFeedback = feedback;
        feedback = feedbackParam * feedbackZeroZone;
        float feedbackInc = (feedback - currentFeedback) * INV_BLOCK_SIZE;

        float currentDelaySize1 = clamp(delaySize1, 0, delayBufferSize);
        delaySize1 = clamp(430 + 70 * quadrant, 0, delayBufferSize);
        float delaySizeInc1 = (delaySize1 - currentDelaySize1) * INV_BLOCK_SIZE;

        float * __restrict__ sp = sampleBlock_;
        float * __restrict__ const delayBuf = delayBuffer_;

        float iirFilter1, iirFilter2, iirFilter3, iirFilter4, iirFilter5, iirFilter6, iirFilter7, iirFilter8;
        float cos, sin;
        float phase2;
        float shifterIn;
        float shifterOutR = 0, shifterOutI = 0, shifterOut = 0, shifterOut2 = 0;

        const float f = 0.73f;
        const float f2 = 0.75f;
        const float f3 = 0.25f;

        // hi pass params
        float filterB2 = 0.2f;
        float filterB = (filterB2 * filterB2 * 0.5f);

        float _in2_b1 = (1 - filterB);
        float _in2_a0 = (1 + _in2_b1 * _in2_b1 * _in2_b1) * 0.5f;

        float hpZeroZone = clamp((1.3f - param1S) * sqrt3(param1 * quadrantSq * 16), 0, 1);
        filterB2 = 0.32f - 0.3f * hpZeroZone * (1 - (feedbackZeroZone * fabsf(feedbackParam) * 0.125f));
        filterB = (filterB2 * filterB2 * 0.5f);

        float _in3_b1 = (1 - filterB);
        float _in3_a0 = (1 + _in3_b1 * _in3_b1 * _in3_b1) * 0.5f;

        float delayOut1;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {
            float monoIn = (*sp + *(sp + 1)) * 0.5f;

            // monoIn lp
            low1 += f * band1;
            band1 += f * (monoIn - low1 - band1);

            // feedback hp
            float feedbackIn = tanh4(feedbackInput * 1.35f) * currentFeedback;

            _ly1 = apcoef3 * (_ly1 + feedbackIn) - _lx1; // allpass
            _lx1 = feedbackIn;

            _ly2 = apcoef4 * (_ly2 + _ly1) - _lx2; // allpass
            _lx2 = _ly1;

            float hp_in_x0 = (feedbackIn + _ly2) * 0.5f;
            hp_in_y0 = _in3_a0 * (hp_in_x0 - hp_in_x1) + _in3_b1 * hp_in_y1;
            hp_in_y1 = hp_in_y0;
            hp_in_x1 = hp_in_x0;

            // feedback lp
            low2 += f2 * band2;
            band2 += f2 * (hp_in_y0 - low2 - band2);
            low3 += f2 * band3;
            band3 += f2 * (low2 - low3 - band3);

            delayWritePos = (delayWritePos + 1) & delayBufferSizeM1;
            delayBuf[delayWritePos] = low3;

            delayReadPos = modulo2(delayWritePos - currentDelaySize1, delayBufferSize);
            delayOut1 = delayInterpolation(delayReadPos, delayBuf, delayBufferSizeM1);

            shifterIn = clamp(low1 - delayOut1, -1.f, 1.f);

            float hp_in2_x0 = shifterIn;
            hp_in2_y0 = _in2_a0 * (hp_in2_x0 - hp_in2_x1) + _in2_b1 * hp_in2_y1;
            hp_in2_y1 = hp_in2_y0;
            hp_in2_x1 = hp_in2_x0;

            low4 += f * band4;
            band4 += f * (hp_in2_y0 - low4 - band4);

            // Frequency shifter

            //     Phase reference path
            iirFilter1 = iirFilter(low4, 0.48645677879491144857f, &hb1_x1, &hb1_x2, &hb1_y1, &hb1_y2);
            iirFilter2 = iirFilter(iirFilter1, 0.88068726735639790704f, &hb2_x1, &hb2_x2, &hb2_y1, &hb2_y2);
            iirFilter3 = iirFilter(iirFilter2, 0.97790456293916316888f, &hb3_x1, &hb3_x2, &hb3_y1, &hb3_y2);
            iirFilter4 = iirFilter(iirFilter3, 0.99767037906310385154f, &hb4_x1, &hb4_x2, &hb4_y1, &hb4_y2);

            //     +90 deg path
            iirFilter5 = iirFilter(low4, 0.16507919004304125177f, &hb5_x1, &hb5_x2, &hb5_y1, &hb5_y2);
            iirFilter6 = iirFilter(iirFilter5, 0.73969068299070206418f, &hb6_x1, &hb6_x2, &hb6_y1, &hb6_y2);
            iirFilter7 = iirFilter(iirFilter6, 0.94788883423814862539f, &hb7_x1, &hb7_x2, &hb7_y1, &hb7_y2);
            iirFilter8 = iirFilter(iirFilter7, 0.99119752093109647628f, &hb8_x1, &hb8_x2, &hb8_y1, &hb8_y2);

            //     sin
            phase1 = phase1 + currentShift;
            if (phase1 >= 1)
            {
                phase1 -= 2.f;
            }
            //     cos = sin( x + 90°)
            phase2 = phase1 + 0.5f;
            if (phase2 >= 1)
            {
                phase2 -= 2.f;
            }

            sin = fastSin(phase1);
            cos = fastSin(phase2);

            shifterOutR = sin * iirFilter4;
            shifterOutI = cos * iirFilter8;
            shifterOut = shifterOutI + shifterOutR;
            shifterOut2 = shifterOutI - shifterOutR;
            float shifterOutMixA = shifterOut2 * shiftPlus;
            float shifterOutMixB = shifterOut * shiftMinus;
            float shifterOutMix = shifterOutMixA + shifterOutMixB;
            float shifterOutMix2 = shifterOutMixA - shifterOutMixB;

            feedbackInput = shifterOutMix;

            // bass boost
            low5 += f3 * band5;
            band5 += f3 * (shifterOutMix - low5 - band5);

            *sp = (*sp * dry) - tanh4((shifterOutMix + low5) * 1.5f) * wetL;
            sp++;
            *sp = (*sp * dry) - tanh4((shifterOutMix2 + low5) * 1.5f) * wetR;
            sp++;

            currentFeedback += feedbackInc;
            currentShift += shiftInc;
            currentDelaySize1 += delaySizeInc1;
        }
    }
    break;
    case FILTER2_WIDE:
    {
        mixerGain_ = 0.02f * gainTmp + .98f * mixerGain_;
        float mixerGain_01 = clamp(mixerGain_, 0, 1);
        int mixerGain255 = mixerGain_01 * 255;
        float dry = panTable[255 - mixerGain255];
        float wet = panTable[mixerGain255] * 0.3333f;
        float extraAmp = clamp(mixerGain_ - 1, 0, 1);
        wet += extraAmp;

        float wetL = wet * (1 + matrixFilterPan);
        float wetR = wet * (1 - matrixFilterPan);

        param1S = 0.02f * fabsf(this->params_.effect2.param1 + matrixFilterFrequency) + .98f * param1S;
        param2S = 0.02f * clamp(this->params_.effect2.param2 + matrixFilterParam2, 0.f, 1.f) + .98f * param2S;

        float detune = param1S * param1S * 0.03125f;

        float currentShift = shift;
        shift = clamp(1 + detune, 0, 16);
        float shiftInc = (shift - currentShift) * INV_BLOCK_SIZE;

        float currentShift2 = shift2;
        shift2 = clamp(1 - detune, 0, 16);
        float shiftInc2 = (shift2 - currentShift2) * INV_BLOCK_SIZE;

        const int delaySizeInt = 500;

        const float f = 0.72f;
        const float f2 = 0.7f + f * 0.1f;
        const float f3 = 0.32f;

        // hi pass params
        float filterB2 = clamp(param2S * param2S * 0.93f, 0, 1);
        float filterB = (filterB2 * filterB2 * 0.5f);

        const float _hp_b1 = (1 - filterB);
        const float _hp_a0 = (1 + _hp_b1 * _hp_b1 * _hp_b1) * 0.5f;

        float * __restrict__ sp = sampleBlock_;
        float * __restrict__ const delayBuf = delayBuffer_;

        float level1, level2, level3, level4;
        float delayReadPos180;
        float delayOut1 = 0, delayOut2 = 0, delayOut3 = 0, delayOut4 = 0;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {
            float monoIn = (*sp + *(sp + 1)) * 0.5f;

            low1 += f * band1;
            band1 += f * (monoIn - low1 - band1);

            delayWritePos = (delayWritePos + 1) & delayBufStereoSizeM1;
            delayWritePosF = (float)delayWritePos;

            // predelay
            delayBuf[delayBufStereoSize + delayWritePos] = low1;
            int wp = (delayWritePos - delaySizeInt) & delayBufStereoSizeM1;

            // hp
            float hp_in_x0 = delayBuf[delayBufStereoSize + wp];
            hp_in_y0 = _hp_a0 * (hp_in_x0 - hp_in_x1) + _hp_b1 * hp_in_y1;
            hp_in_y1 = hp_in_y0;
            hp_in_x1 = hp_in_x0;

            delayBuf[delayWritePos] = hp_in_y0;

            float hpComplement = hp_in_x0 - hp_in_y0;

            //--------------- shifter 1

            delayReadPos = modulo(delayReadPos + currentShift, delayBufStereoSize);
            delayReadPos180 = modulo(delayReadPos + delayBufferSize90, delayBufStereoSize);

            delayOut1 = delayInterpolation(delayReadPos, delayBuf, delayBufStereoSizeM1);
            delayOut2 = delayInterpolation(delayReadPos180, delayBuf, delayBufStereoSizeM1);

            float rwp1 = modulo2(delayWritePosF - delayReadPos, delayBufStereoSize);
            float rwp2 = modulo2(delayWritePosF - delayReadPos180, delayBufStereoSize);

            level1 = hann(rwp1 * delayBufStereoSizeInv);
            level2 = hann(rwp2 * delayBufStereoSizeInv);

            float out1 = delayOut1 * (level1);
            float out2 = delayOut2 * (level2);

            //--------------- shifter 2

            delayReadPos2 = modulo(delayReadPos2 + currentShift2, delayBufStereoSize);
            delayReadPos180 = modulo(delayReadPos2 + delayBufferSize90, delayBufStereoSize);

            delayOut3 = delayInterpolation(delayReadPos2, delayBuf, delayBufStereoSizeM1);
            delayOut4 = delayInterpolation(delayReadPos180, delayBuf, delayBufStereoSizeM1);

            float rwp3 = modulo2(delayWritePosF - delayReadPos2, delayBufStereoSize);
            float rwp4 = modulo2(delayWritePosF - delayReadPos180, delayBufStereoSize);

            level3 = hann(rwp3 * delayBufStereoSizeInv);
            level4 = hann(rwp4 * delayBufStereoSizeInv);

            float out3 = delayOut3 * (level3);
            float out4 = delayOut4 * (level4);

            // lp output
            low3 += f2 * band3;
            band3 += f2 * ((hpComplement + out1 + out2) - low3 - band3);

            low4 += f2 * band4;
            band4 += f2 * ((hpComplement + out3 + out4) - low4 - band4);

            // bass boost
            low5 += f3 * band5;
            band5 += f3 * (low3 - low5 - band5);
            low6 += f3 * band6;
            band6 += f3 * (low4 - low6 - band6);

            *sp = *sp * dry + (low3 + low3 + low5) * wetL;
            sp++;
            *sp = *sp * dry + (low4 + low4 + low6) * wetR;
            sp++;

            currentShift += shiftInc;
            currentShift2 += shiftInc2;
        }
    }
    break;
    case FILTER2_DELAYCRUNCH:
    {
        mixerGain_ = 0.02f * gainTmp + .98f * mixerGain_;
        float mixerGain_01 = clamp(mixerGain_, 0, 1);
        int mixerGain255 = mixerGain_01 * 255;
        float dry = panTable[255 - mixerGain255];
        float wet = panTable[mixerGain255] * 1.2f;
        float extraAmp = clamp(mixerGain_ - 1, 0, 1);
        wet += extraAmp;

        float wetL = wet * (1 + matrixFilterPan);
        float wetR = wet * (1 - matrixFilterPan);

        param1S = 0.005f * (this->params_.effect2.param1) + .995f * param1S;
        matrixFilterFrequencyS = 0.02f * (matrixFilterFrequency * matrixFilterFrequency) + .98f * matrixFilterFrequencyS;
        param2S = 0.05f * (this->params_.effect2.param2 + matrixFilterParam2) + .95f * param2S;

        param2S = clamp(param2S, 0, 1.f);
        float param2Square = param2S * param2S;

        feedback = param2S * 1.2f;

        const float sampleRateDivide = 4;
        const float sampleRateDivideInv = 1 / sampleRateDivide;
        float inputIncCount = 0;

        float currentDelaySize1 = clamp(delaySize1, 0, delayBufferSize);
        delaySize1 = 1.f + (delayBufferSize - 16) * clamp(param1S + (matrixFilterFrequencyS * 0.0625f), 0.f, 1.f);
        float delaySizeInc1 = (delaySize1 - currentDelaySize1) * sampleRateDivideInv * INV_BLOCK_SIZE;

        // hp input
        const float filterB2 = 0.1f + param2Square * 0.4f;
        const float filterB = (filterB2 * filterB2 * 0.5f);
        const float _in2_b1 = (1 - filterB);
        const float _in2_a0 = (1 + _in2_b1 * _in2_b1 * _in2_b1) * 0.5f;

        // hp feedback
        float filterC2 = 0.1f + param2Square * 0.3f;
        float filterC = (filterC2 * filterC2 * 0.5f);

        float _in3_b1 = (1 - filterC);
        float _in3_a0 = (1 + _in3_b1 * _in3_b1 * _in3_b1) * 0.5f;

        const float f = 0.8f + param2Square * 0.198f;
        const float f2 = 0.9f - param2Square * 0.15f;

        const float fnotch = 1.03f;

        float * __restrict__ sp = sampleBlock_;
        float * __restrict__ const delayBuf = delayBuffer_;

        float delayOut1 = 0;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {

            if (++inputIncCount >= sampleRateDivide)
            {
                float monoIn = (*sp + *(sp + 1)) * 0.5f;

                inputIncCount = 0;
                delayWritePos = (delayWritePos + 1) & delayBufferSizeM1;
                delayWritePosF = (float)delayWritePos;

                // hp
                hb4_x1 = monoIn;
                hb4_y1 = _in2_a0 * (hb4_x1 - hb4_x2) + _in2_b1 * hb4_y2;
                hb4_y2 = hb4_y1;
                hb4_x2 = hb4_x1;

                hb6_x1 = hb4_y1;
                hb6_y1 = _in2_a0 * (hb6_x1 - hb6_x2) + _in2_b1 * hb6_y2;
                hb6_y2 = hb6_y1;
                hb6_x2 = hb6_x1;

                float hp_in_x0 = hb6_y1;
                hp_in_y0 = _in2_a0 * (hp_in_x0 - hp_in_x1) + _in2_b1 * hp_in_y1;
                hp_in_y1 = hp_in_y0;
                hp_in_x1 = hp_in_x0;

                // lp & hp feedback
                low1 += f * band1;
                band1 += f * ((delayOut1 * feedback) - low1 - band1);

                float hp_in2_x0 = low1;
                hp_in2_y0 = _in3_a0 * (hp_in2_x0 - hp_in2_x1) + _in3_b1 * hp_in2_y1;
                hp_in2_y1 = hp_in2_y0;
                hp_in2_x1 = hp_in2_x0;

                hb5_x1 = hp_in2_y0 + hp_in_y0;
                hb5_y1 = _in3_a0 * (hb5_x1 - hb5_x2) + _in3_b1 * hb5_y2;
                hb5_y2 = hb5_y1;
                hb5_x2 = hb5_x1;

                delayBuf[delayWritePos] = hb5_y1;
            }

            delayReadPos = modulo2(delayWritePosF - currentDelaySize1, delayBufferSize);
            delayOut1 = delayInterpolation(delayReadPos, delayBuf, delayBufferSizeM1);
            delayWritePosF += sampleRateDivideInv;

            // lp
            low3 += f2 * band3;
            band3 += f2 * (delayOut1 - low3 - band3);
            low4 += f2 * band4;
            band4 += f2 * (low3 - low4 - band4);
            low5 += f2 * band5;
            band5 += f2 * (low4 - low5 - band5);
            low6 += f2 * band6;
            band6 += f2 * (low5 - low6 - band6);

            // notch
            hb1_x1 += fnotch * hb1_y1;
            float high7 = low6 - hb1_x1 - hb1_y1;
            hb1_y1 += fnotch * high7;
            float notch = (high7 + hb1_x1);

            *sp = *sp * dry + notch * wetL;
            sp++;
            *sp = *sp * dry + notch * wetR;
            sp++;

            currentDelaySize1 += delaySizeInc1;
        }
    }
    break;
    case FILTER2_PINGPONG:
    {
        mixerGain_ = 0.02f * gainTmp + .98f * mixerGain_;
        float mixerGain_01 = clamp(mixerGain_, 0, 1);
        int mixerGain255 = mixerGain_01 * 255;
        float dry = panTable[255 - mixerGain255];
        float wet = panTable[mixerGain255] * 1.25f;
        float extraAmp = clamp(mixerGain_ - 1, 0, 1);
        wet += extraAmp;

        float wetL = wet * (1 + matrixFilterPan);
        float wetR = wet * (1 - matrixFilterPan);

        param1S = 0.005f * (this->params_.effect2.param1) + .995f * param1S;
        matrixFilterFrequencyS = 0.02f * (matrixFilterFrequency * matrixFilterFrequency) + .98f * matrixFilterFrequencyS;

        // Interpolation linéaire de param2S sur le bloc pour éviter l'aliasing
        // lors d'une modulation rapide (enveloppe). La valeur cible est calculée
        // une fois, puis param2S est incrémenté sample par sample dans la boucle.
        float param2Target = clamp(0.05f * (this->params_.effect2.param2 + matrixFilterParam2) + .95f * param2S, 0.f, 1.f);
        float param2SInc = (param2Target - param2S) * INV_BLOCK_SIZE;

        // Calcul des paramètres du filtre à partir de la valeur courante de param2S
        // (sera mis à jour chaque sample dans la boucle)
        float param2Square = param2S * param2S;

        const float sampleRateDivide = 4;
        const float sampleRateDivideInv = 1 / sampleRateDivide;
        float inputIncCount = 0;

        float currentDelaySize1 = clamp(delaySize1, 0, delayBufferSizeF);
        delaySize1 = 1.f + (delayBufferSize - 16) * clamp(param1S + (matrixFilterFrequencyS * 0.0625f), 0.f, 1.f);
        float delaySizeInc1 = (delaySize1 - currentDelaySize1) * sampleRateDivideInv * INV_BLOCK_SIZE;

        float currentDelaySize2 = clamp(delaySize2, 0, delayBufferSizeF);
        delaySize2 = delaySize1 * 0.5f;
        float delaySizeInc2 = (delaySize2 - currentDelaySize2) * sampleRateDivideInv * INV_BLOCK_SIZE;

        // hp 1 — coefficients calculés avant la boucle, stables sur le bloc
        float filterB2 = 0.1f + param2Square * 0.43f;
        float filterB = (filterB2 * filterB2 * 0.5f);
        float _in2_b1 = (1 - filterB);
        float _in2_a0 = (1 + _in2_b1 * _in2_b1 * _in2_b1) * 0.5f;

        // hp 2
        float filterC2 = 0.1f + param2Square * 0.33f;
        float filterC = (filterC2 * filterC2 * 0.5f);

        float _in3_b1 = (1 - filterC);
        float _in3_a0 = (1 + _in3_b1 * _in3_b1 * _in3_b1) * 0.5f;

        float f = 0.8f + param2Square * 0.198f;
        float airAttn = (param2Square + param2Square + param1S) * 0.03f;
        float f2 = 0.85f - airAttn;
        float f3 = f2 - 0.03f;
        const float fnotch = 1.03f;

        float delayOut1, delayOut2;
        float * __restrict__ sp = sampleBlock_;
        float * __restrict__ const delayBuf = delayBuffer_;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {
            // Lissage sample par sample de param2S pour éviter l'aliasing
            param2S += param2SInc;
            param2Square = param2S * param2S;

            // f2/f3 servent à 48 kHz : mis à jour à chaque sample (coût négligeable)
            airAttn = (param2Square + param2Square + param1S) * 0.03f;
            f2 = 0.85f - airAttn;
            f3 = f2 - 0.03f;

            if (++inputIncCount >= sampleRateDivide)
            {
                // Coefficients utilisés uniquement à SR/4 : calculés ici (4x moins souvent)
                feedback = param2S * 1.22f;
                filterB2 = 0.1f + param2Square * 0.43f;
                filterB = (filterB2 * filterB2 * 0.5f);
                _in2_b1 = (1 - filterB);
                _in2_a0 = (1 + _in2_b1 * _in2_b1 * _in2_b1) * 0.5f;
                filterC2 = 0.1f + param2Square * 0.33f;
                filterC = (filterC2 * filterC2 * 0.5f);
                _in3_b1 = (1 - filterC);
                _in3_a0 = (1 + _in3_b1 * _in3_b1 * _in3_b1) * 0.5f;
                f = 0.8f + param2Square * 0.198f;

                float monoIn = (*sp + *(sp + 1)) * 0.5f;

                inputIncCount = 0;
                delayWritePos = (delayWritePos + 1) & delayBufferSizeM1;
                delayWritePosF = (float)delayWritePos;

                low1 += f * band1;
                band1 += f * ((feedbackInput * feedback) - low1 - band1);

                // hp
                hb4_x1 = monoIn;
                hb4_y1 = _in2_a0 * (hb4_x1 - hb4_x2) + _in2_b1 * hb4_y2;
                hb4_y2 = hb4_y1;
                hb4_x2 = hb4_x1;

                hb6_x1 = hb4_y1;
                hb6_y1 = _in2_a0 * (hb6_x1 - hb6_x2) + _in2_b1 * hb6_y2;
                hb6_y2 = hb6_y1;
                hb6_x2 = hb6_x1;

                float hp_in_x0 = hb6_y1;
                hp_in_y0 = _in2_a0 * (hp_in_x0 - hp_in_x1) + _in2_b1 * hp_in_y1;
                hp_in_y1 = hp_in_y0;
                hp_in_x1 = hp_in_x0;

                float hp_in2_x0 = low1;
                hp_in2_y0 = _in3_a0 * (hp_in2_x0 - hp_in2_x1) + _in3_b1 * hp_in2_y1;
                hp_in2_y1 = hp_in2_y0;
                hp_in2_x1 = hp_in2_x0;

                hb5_x1 = hp_in2_y0 + hp_in_y0;
                hb5_y1 = _in3_a0 * (hb5_x1 - hb5_x2) + _in3_b1 * hb5_y2;
                hb5_y2 = hb5_y1;
                hb5_x2 = hb5_x1;

                delayBuf[delayWritePos] = hb5_y1;
            }

            delayReadPos = modulo2(delayWritePosF - currentDelaySize1, delayBufferSize);
            delayOut1 = delayInterpolation(delayReadPos, delayBuf, delayBufferSizeM1);
            delayReadPos = modulo2(delayWritePosF - currentDelaySize2, delayBufferSize);
            delayOut2 = delayInterpolation(delayReadPos, delayBuf, delayBufferSizeM1);

            feedbackInput = delayOut1;

            delayWritePosF += sampleRateDivideInv;

            // lp L
            low3 += f2 * band3;
            band3 += f2 * (delayOut1 - low3 - band3);
            low4 += f2 * band4;
            band4 += f2 * (low3 - low4 - band4);
            low5 += f2 * band5;
            band5 += f2 * (low4 - low5 - band5);
            low6 += f2 * band6;
            band6 += f2 * (low5 - low6 - band6);

            // lp R
            hb2_x1 += f3 * hb2_y1;
            hb2_y1 += f3 * (delayOut2 - hb2_x1 - hb2_y1);
            hb2_x2 += f3 * hb2_y2;
            hb2_y2 += f3 * (hb2_x1 - hb2_x2 - hb2_y2);
            hb3_x1 += f3 * hb3_y1;
            hb3_y1 += f3 * (hb2_x2 - hb3_x1 - hb3_y1);
            hb3_x2 += f3 * hb3_y2;
            hb3_y2 += f3 * (hb3_x1 - hb3_x2 - hb3_y2);

            // notch L
            hb1_x1 += fnotch * hb1_y1;
            float high7 = low6 - hb1_x1 - hb1_y1;
            hb1_y1 += fnotch * high7;
            float notchL = (high7 + hb1_x1);

            // notch R
            hb1_x2 += fnotch * hb1_y2;
            float high8 = hb3_x2 - hb1_x2 - hb1_y2;
            hb1_y2 += fnotch * high8;
            float notchR = (high8 + hb1_x2);

            *sp = *sp * dry + notchL * wetL;
            sp++;
            *sp = *sp * dry + notchR * wetR;
            sp++;

            currentDelaySize1 += delaySizeInc1;
            currentDelaySize2 += delaySizeInc2;
        }
    }
    break;
    case FILTER2_DIFFUSER:
    {
        mixerGain_ = 0.02f * gainTmp + .98f * mixerGain_;
        float mixerGain_01 = clamp(mixerGain_, 0, 1);
        int mixerGain255 = mixerGain_01 * 255;
        float dry = panTable[255 - mixerGain255];
        float wet = panTable[mixerGain255];
        float extraAmp = clamp(mixerGain_ - 1, 0, 1);
        wet += extraAmp;

        float wetL = wet * (1 + matrixFilterPan);
        float wetR = wet * (1 - matrixFilterPan);

        param1S = 0.005f * (this->params_.effect2.param1) + .995f * param1S;
        matrixFilterFrequencyS = 0.01f * (matrixFilterFrequency) + .99f * matrixFilterFrequencyS;
        param2S = 0.05f * (this->params_.effect2.param2 + matrixFilterParam2) + .95f * param2S;

        param1S = clamp(param1S, 0, 1);
        param2S = clamp(param2S, 0, 1);

        feedback = param2S * 0.4999f;

        const float sampleRateDivide = 4;
        const float sampleRateDivideInv = 1 / sampleRateDivide;
        float inputIncCount = 0;

        float currentDelaySize1 = clamp(delaySize1, 0, delayBufStereoSize);
        delaySize1 = 1.f + (delayBufStereoSize - 120) * clamp(param1S + (matrixFilterFrequencyS * 0.125f), 0.f, 1.f) * 0.5f;
        float delaySizeInc1 = (delaySize1 - currentDelaySize1) * sampleRateDivideInv * INV_BLOCK_SIZE;

        float filterB2 = 0.15f + param2S * 0.1f;
        float filterB = (filterB2 * filterB2 * 0.5f);

        float _in3_b1 = (1 - filterB);
        float _in3_a0 = (1 + _in3_b1 * _in3_b1 * _in3_b1) * 0.5f;

        float f = 0.7f - clamp(feedback - (param1S * 0.3f), 0, 1) * 0.6f;
        float f2 = 0.75f - param1S * 0.1f;
        const float fnotch = 1.03f;

        const float inputCoef1 = 0.7f;
        const float inputCoef2 = 0.625f;

        float diff1Out = 0, diff2Out = 0, diff3Out = 0, diff4Out = 0, diff5Out = 0;

        float sizeParamInpt = 0.2f + param1S * 0.2f;

        // Modulation lente et indépendante de chaque allpass — brise les résonances métalliques.
        // shift/shift2 : marcheurs aléatoires bornés (± ~2 samples à SR/4), signes alternés
        // pour éviter un pitch-shift global uniforme.
        shift  = shift  * 0.9985f + noise[4] * 0.04f;
        shift2 = shift2 * 0.9990f + noise[5] * 0.03f;
        const float modAmt = 4.5f;

        float inputBuffer1ReadLen = inputBufferLen1 * sizeParamInpt + shift  * modAmt;
        float inputBuffer2ReadLen = inputBufferLen2 * sizeParamInpt - shift2 * modAmt;
        float inputBuffer3ReadLen = inputBufferLen3 * sizeParamInpt + shift2 * modAmt * 0.7f;
        float inputBuffer4ReadLen = inputBufferLen4 * sizeParamInpt - shift  * modAmt * 0.5f;
        float inputBuffer5ReadLen = inputBufferLen5 * sizeParamInpt + (shift - shift2) * modAmt * 0.4f;

        // Offsets de lecture entiers, constants sur le bloc : (int)modulo2(wp - len) == wrap(wp - ceil(len))
        const int readOff1 = (int)ceilf(inputBuffer1ReadLen);
        const int readOff2 = (int)ceilf(inputBuffer2ReadLen);
        const int readOff3 = (int)ceilf(inputBuffer3ReadLen);
        const int readOff4 = (int)ceilf(inputBuffer4ReadLen);
        const int readOff5 = (int)ceilf(inputBuffer5ReadLen);

        const int bufferStart1 = delayBufStereoSize;
        const int bufferStart2 = bufferStart1 + inputBufferLen1;
        const int bufferStart3 = bufferStart2 + inputBufferLen2;
        const int bufferStart4 = bufferStart3 + inputBufferLen3;
        const int bufferStart5 = bufferStart4 + inputBufferLen4;

        float * __restrict__ sp = sampleBlock_;
        float * __restrict__ const delayBuf = delayBuffer_;
        float delayOut1, delayOut2;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {

            if (++inputIncCount >= sampleRateDivide)
            {
                float monoIn = (*sp + *(sp + 1)) * 0.5f;

                inputIncCount = 0;
                delayWritePos = (delayWritePos + 1) & delayBufStereoSizeM1;
                delayWritePosF = (float)delayWritePos;
                if (++inputWritePos1 >= inputBufferLen1) inputWritePos1 = 0;
                if (++inputWritePos2 >= inputBufferLen2) inputWritePos2 = 0;
                if (++inputWritePos3 >= inputBufferLen3) inputWritePos3 = 0;
                if (++inputWritePos4 >= inputBufferLen4) inputWritePos4 = 0;
                if (++inputWritePos5 >= inputBufferLen5) inputWritePos5 = 0;

                // ---- feedback lp
                low1 += f * band1;
                band1 += f * ((feedbackInput * feedback) - low1 - band1);

                // ---- diffuser 1
                int rp1 = inputWritePos1 - readOff1;
                if (rp1 < 0) rp1 += inputBufferLen1;
                int inputReadPos1 = bufferStart1 + rp1;
                float in_apSum1 = (monoIn - low1) + delayBuf[inputReadPos1] * inputCoef1;
                diff1Out = delayBuf[inputReadPos1] - in_apSum1 * inputCoef1;
                delayBuf[bufferStart1 + inputWritePos1] = in_apSum1;

                // ---- diffuser 2
                int rp2 = inputWritePos2 - readOff2;
                if (rp2 < 0) rp2 += inputBufferLen2;
                int inputReadPos2 = bufferStart2 + rp2;
                float in_apSum2 = diff1Out + delayBuf[inputReadPos2] * inputCoef2;
                diff2Out = delayBuf[inputReadPos2] - in_apSum2 * inputCoef2;
                delayBuf[bufferStart2 + inputWritePos2] = in_apSum2;

                // ---- diffuser 3
                int rp3 = inputWritePos3 - readOff3;
                if (rp3 < 0) rp3 += inputBufferLen3;
                int inputReadPos3 = bufferStart3 + rp3;
                float in_apSum3 = -diff2Out + delayBuf[inputReadPos3] * inputCoef2;
                diff3Out = delayBuf[inputReadPos3] - in_apSum3 * inputCoef2;
                delayBuf[bufferStart3 + inputWritePos3] = in_apSum3;

                // ---- diffuser 4
                int rp4 = inputWritePos4 - readOff4;
                if (rp4 < 0) rp4 += inputBufferLen4;
                int inputReadPos4 = bufferStart4 + rp4;
                float in_apSum4 = diff3Out + delayBuf[inputReadPos4] * inputCoef2;
                diff4Out = delayBuf[inputReadPos4] - in_apSum4 * inputCoef2;
                delayBuf[bufferStart4 + inputWritePos4] = in_apSum4;

                // ---- diffuser 5
                int rp5 = inputWritePos5 - readOff5;
                if (rp5 < 0) rp5 += inputBufferLen5;
                int inputReadPos5 = bufferStart5 + rp5;
                float in_apSum5 = -diff4Out + delayBuf[inputReadPos5] * inputCoef2;
                diff5Out = delayBuf[inputReadPos5] - in_apSum5 * inputCoef2;
                delayBuf[bufferStart5 + inputWritePos5] = in_apSum5;

                low1 += f2 * band1;
                band1 += f2 * (diff5Out - low1 - band1);

                // hp
                float hp_in_x0 = diff5Out;
                hp_in_y0 = _in3_a0 * (hp_in_x0 - hp_in_x1) + _in3_b1 * hp_in_y1;
                hp_in_y1 = hp_in_y0;
                hp_in_x1 = hp_in_x0;

                float hp_in2_x0 = hp_in_y0;
                hp_in2_y0 = _in3_a0 * (hp_in2_x0 - hp_in2_x1) + _in3_b1 * hp_in2_y1;
                hp_in2_y1 = hp_in2_y0;
                hp_in2_x1 = hp_in2_x0;

                hb5_x1 = hp_in2_y0;
                hb5_y1 = _in3_a0 * (hb5_x1 - hb5_x2) + _in3_b1 * hb5_y2;
                hb5_y2 = hb5_y1;
                hb5_x2 = hb5_x1;

                low2 += f * band2;
                band2 += f * (hb5_y1 - low2 - band2);

                delayBuf[delayWritePos] = low1 + low2;
            }

            delayReadPos = modulo2(delayWritePosF - currentDelaySize1, delayBufStereoSize);
            delayOut1 = delayInterpolation(delayReadPos, delayBuf, delayBufStereoSizeM1);

            delayReadPos = modulo2(delayReadPos - 512, delayBufStereoSize);
            delayOut2 = delayInterpolation(delayReadPos, delayBuf, delayBufStereoSizeM1);

            feedbackInput = delayOut1;

            delayWritePosF += sampleRateDivideInv;

            // lp L
            low3 += f2 * band3;
            band3 += f2 * (delayOut1 - low3 - band3);
            low4 += f2 * band4;
            band4 += f2 * (low3 - low4 - band4);
            low5 += f2 * band5;
            band5 += f2 * (low4 - low5 - band5);
            low6 += f2 * band6;
            band6 += f2 * (low5 - low6 - band6);

            // lp R
            hb2_x1 += f2 * hb2_y1;
            hb2_y1 += f2 * (delayOut2 - hb2_x1 - hb2_y1);
            hb2_x2 += f2 * hb2_y2;
            hb2_y2 += f2 * (hb2_x1 - hb2_x2 - hb2_y2);
            hb3_x1 += f2 * hb3_y1;
            hb3_y1 += f2 * (hb2_x2 - hb3_x1 - hb3_y1);
            hb3_x2 += f2 * hb3_y2;
            hb3_y2 += f2 * (hb3_x1 - hb3_x2 - hb3_y2);

            // notch L
            hb1_x1 += fnotch * hb1_y1;
            float high7 = low6 - hb1_x1 - hb1_y1;
            hb1_y1 += fnotch * high7;
            float notchL = (high7 + hb1_x1);

            // notch R
            hb1_x2 += fnotch * hb1_y2;
            float high8 = hb3_x2 - hb1_x2 - hb1_y2;
            hb1_y2 += fnotch * high8;
            float notchR = (high8 + hb1_x2);

            *sp = *sp * dry + notchL * wetL;
            sp++;
            *sp = *sp * dry + notchR * wetR;
            sp++;

            currentDelaySize1 += delaySizeInc1;
        }
    }
    break;
    case FILTER2_GRAIN1:
    {
        mixerGain_ = 0.02f * gainTmp + .98f * mixerGain_;
        float mixerGain_01 = clamp(mixerGain_, 0, 1);
        int mixerGain255 = mixerGain_01 * 255;
        float dry = panTable[255 - mixerGain255];
        float wet = panTable[mixerGain255] * 1.5f;
        float extraAmp = clamp(mixerGain_ - 1, 0, 1);
        wet += extraAmp;

        float wetL = wet * (1 + matrixFilterPan);
        float wetR = wet * (1 - matrixFilterPan);

        float param2 = clamp(fabsf(this->params_.effect2.param2 + matrixFilterParam2), 0, 1);
        param2 *= param2;

        float lock = ((param2 > 0.99f) ? 0.f : 1.f);
        lockA = lockA * 0.98f + lock * 0.02f;
        lockB = (1 - lockA);

        matrixFilterFrequencyS = 0.02f * (matrixFilterFrequency * matrixFilterFrequency * 0.125f) + 0.98f * matrixFilterFrequencyS;
        param1S = 0.02f * fabs(this->params_.effect2.param1 + matrixFilterFrequencyS) + .98f * param1S;

        if (lockA >= 0.9999f)
        {
            param2S = 0.05f * param2 + 0.95f * param2S;
        }

        const float sampleRateDivide = 4;
        const float sampleRateDivideInv = 1 / sampleRateDivide;
        float inputIncCount = 0;

        bool grainProb = 0.985f >= fabs(noise[0]);

        const float filterB2 = 0.25f;
        const float filterB = (filterB2 * filterB2 * 0.5f);

        const float _in3_b1 = (1 - filterB);
        const float _in3_a0 = (1 + _in3_b1 * _in3_b1 * _in3_b1) * 0.5f;

        if (grainProb)
        {
            if (grainTable[grainNext][GRAIN_RAMP] >= 1 && grainTable[grainPrev][GRAIN_RAMP] > (0.33f - param2S * 0.27f))
            {
                // grain done, compute another one
                float jitter = foldAbs(param2S);
                float param2sq = sqrt3(param2S);
                float grainRate = sampleRateDivideInv * (1 + jitter * noise[4] * 0.0025f);
                grainTable[grainNext][GRAIN_RAMP] = 0;
                grainTable[grainNext][GRAIN_SIZE] = clamp((1800 + (noise[2]) * 40 * jitter * jitter) * param1S * param1S, 432, delayBufferSize - 100);
                grainTable[grainNext][GRAIN_POS] = modulo2(delayWritePosF - (400 + param2sq * noise[5] * 1290), delayBufferSize);
                float invGrainSize = 1 / grainTable[grainNext][GRAIN_SIZE];
                grainTable[grainNext][GRAIN_CURRENT_SHIFT] = grainTable[grainNext][GRAIN_NEXT_SHIFT];
                grainTable[grainNext][GRAIN_NEXT_SHIFT] = grainRate * invGrainSize;
                grainTable[grainNext][GRAIN_INC] = clamp((grainTable[grainNext][GRAIN_NEXT_SHIFT] - grainTable[grainNext][GRAIN_CURRENT_SHIFT]) * invGrainSize, 0, 0.5f);
                grainTable[grainNext][GRAIN_VOL] = clamp(0.25f + sqrt3(fabsf(noise[3])), 0, 1);
                grainTable[grainNext][GRAIN_PAN] = clamp(1 + (noise[6]) * param2sq, 0, 2) * 0.5f;
                grainPrev = grainNext;

                if (lockB < 0.02f)
                {
                    loopSize = clamp((1 - param1S) * 1400, 48, 1800);
                }
            }
            if (++grainNext > 2)
            {
                grainNext = 0;
            }
        }

        float env;

        const float f2 = 0.85f;
        const float fnotch = 1.03f;

        float grain1, grain2, grain3;
        float grain1L, grain1R;
        float grain2L, grain2R;
        float grain3L, grain3R;
        float grainSumR, grainSumL;

        float * __restrict__ sp = sampleBlock_;
        float * __restrict__ const delayBuf = delayBuffer_;
        float monoIn;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {

            if (++inputIncCount >= sampleRateDivide)
            {
                inputIncCount = 0;
                delayWritePos = (delayWritePos + 1) & delayBufferSizeM1;
                delayWritePosF = (float)delayWritePos;

                delayReadPos = modulo(delayWritePos + loopSize, delayBufferSize);
                float feedbackIn = delayInterpolation(delayReadPos, delayBuf, delayBufferSizeM1);

                monoIn = (*sp + *(sp + 1)) * 0.5f;

                // hp
                float hp_in_x0 = monoIn;
                hp_in_y0 = _in3_a0 * (hp_in_x0 - hp_in_x1) + _in3_b1 * hp_in_y1;
                hp_in_y1 = hp_in_y0;
                hp_in_x1 = hp_in_x0;

                delayBuf[delayWritePos] = hp_in_y0 * lockA + feedbackIn * lockB;
            }

            ///-------- grain 1
            grain1L = grain1R = 0;
            if (grainTable[0][GRAIN_RAMP] < 1)
            {
                delayReadPos = modulo(grainTable[0][GRAIN_POS] + grainTable[0][GRAIN_RAMP] * grainTable[0][GRAIN_SIZE], delayBufferSize);
                env = hann(grainTable[0][GRAIN_RAMP]) * grainTable[0][GRAIN_VOL];
                grain1 = env * delayInterpolation(delayReadPos, delayBuf, delayBufferSizeM1);
                grain1L = grain1 * grainTable[0][GRAIN_PAN];
                grain1R = grain1 - grain1L;
            }
            grainTable[0][GRAIN_RAMP] += grainTable[0][GRAIN_CURRENT_SHIFT];
            grainTable[0][GRAIN_CURRENT_SHIFT] += grainTable[0][GRAIN_INC];

            ///-------- grain 2
            grain2L = grain2R = 0;
            if (grainTable[1][GRAIN_RAMP] < 1)
            {
                delayReadPos = modulo(grainTable[1][GRAIN_POS] + grainTable[1][GRAIN_RAMP] * grainTable[1][GRAIN_SIZE], delayBufferSize);
                env = hann(grainTable[1][GRAIN_RAMP]) * grainTable[1][GRAIN_VOL];
                grain2 = env * delayInterpolation(delayReadPos, delayBuf, delayBufferSizeM1);
                grain2L = grain2 * grainTable[1][GRAIN_PAN];
                grain2R = grain2 - grain2L;
            }
            grainTable[1][GRAIN_RAMP] += grainTable[1][GRAIN_CURRENT_SHIFT];
            grainTable[1][GRAIN_CURRENT_SHIFT] += grainTable[1][GRAIN_INC];

            ///-------- grain 3
            grain3L = grain3R = 0;
            if (grainTable[2][GRAIN_RAMP] < 1)
            {
                delayReadPos = modulo(grainTable[2][GRAIN_POS] + grainTable[2][GRAIN_RAMP] * grainTable[2][GRAIN_SIZE], delayBufferSize);
                env = hann(grainTable[2][GRAIN_RAMP]) * grainTable[2][GRAIN_VOL];
                grain3 = env * delayInterpolation(delayReadPos, delayBuf, delayBufferSizeM1);
                grain3L = grain3 * grainTable[2][GRAIN_PAN];
                grain3R = grain3 - grain3L;
            }
            grainTable[2][GRAIN_RAMP] += grainTable[2][GRAIN_CURRENT_SHIFT];
            grainTable[2][GRAIN_CURRENT_SHIFT] += grainTable[2][GRAIN_INC];

            grainSumL = grain1L + grain2L + grain3L;
            grainSumR = grain1R + grain2R + grain3R;

            // lp L
            low1 += f2 * band1;
            band1 += f2 * ((grainSumL)-low1 - band1);
            low3 += f2 * band3;
            band3 += f2 * ((low1)-low3 - band3);
            low5 += f2 * band5;
            band5 += f2 * ((low3)-low5 - band5);

            // lp R
            low2 += f2 * band2;
            band2 += f2 * ((grainSumR)-low2 - band2);
            low4 += f2 * band4;
            band4 += f2 * ((low2)-low4 - band4);
            low6 += f2 * band6;
            band6 += f2 * ((low4)-low6 - band6);

            // notch L
            hb1_x1 += fnotch * hb1_y1;
            float high7 = low5 - hb1_x1 - hb1_y1;
            hb1_y1 += fnotch * high7;
            float notchL = (high7 + hb1_x1);

            // notch R
            hb1_x2 += fnotch * hb1_y2;
            float high8 = low6 - hb1_x2 - hb1_y2;
            hb1_y2 += fnotch * high8;
            float notchR = (high8 + hb1_x2);

            *sp = *sp * dry + notchL * wetL;
            sp++;
            *sp = *sp * dry + notchR * wetR;
            sp++;
        }
    }

    break;
    case FILTER2_GRAIN2:
    {
        mixerGain_ = 0.02f * gainTmp + .98f * mixerGain_;
        float mixerGain_01 = clamp(mixerGain_, 0, 1);
        int mixerGain255 = mixerGain_01 * 255;
        float dry = panTable[255 - mixerGain255];
        float wet = panTable[mixerGain255] * 1.5f;
        float extraAmp = clamp(mixerGain_ - 1, 0, 1);
        wet += extraAmp;

        float wetL = wet * (1 + matrixFilterPan);
        float wetR = wet * (1 - matrixFilterPan);

        float param2 = clamp(fabsf(this->params_.effect2.param2 + matrixFilterParam2), 0, 1);
        param2 *= param2;
        matrixFilterFrequencyS = 0.02f * (matrixFilterFrequency * matrixFilterFrequency * 0.125f) + 0.98f * matrixFilterFrequencyS;
        param1S = 0.02f * fabs(this->params_.effect2.param1 + matrixFilterFrequencyS) + .98f * param1S;
   
        float lock = ((param2 > 0.99f) ? 0.f : 1.f);
        lockA = lockA * 0.98f + lock * 0.02f;
        lockB = (1 - lockA);

        if (lockA >= 0.9999f)
        {
            param1S = 0.0005f * fabs(this->params_.effect2.param1 + matrixFilterFrequency) + .9995f * param1S;
            param2S = 0.05f * param2 + 0.95f * param2S;
        }
        else
        {
            param1S = 0.005f * fabs(this->params_.effect2.param1 + matrixFilterFrequency) + .995f * param1S;
        }

        const float sampleRateDivide = 4;
        const float sampleRateDivideInv = 1 / sampleRateDivide;
        float inputIncCount = 0;

        bool grainProb = clamp(param1S * 3, 0.05f, 0.95f) >= fabs(noise[0]);

        const float filterB2 = 0.25f;
        const float filterB = (filterB2 * filterB2 * 0.5f);

        const float _in3_b1 = (1 - filterB);
        const float _in3_a0 = (1 + _in3_b1 * _in3_b1 * _in3_b1) * 0.5f;

        if (grainProb)
        {
            if (grainTable[grainNext][GRAIN_RAMP] >= 1 && grainTable[grainPrev][GRAIN_RAMP] > (0.33f - param2S * 0.27f))
            {
                // grain done, compute another one
                float param2sq = sqrt3(param2S);
                float grainRate = sampleRateDivideInv * clamp(0.5f + param1S + param2S * noise[4] * 0.025f, 0, 2);
                grainTable[grainNext][GRAIN_RAMP] = 0;
                grainTable[grainNext][GRAIN_SIZE] = clamp((1800 + (noise[2]) * 40 * param2S * param2S) * param1S * param1S, 432, delayBufferSize - 100);
                grainTable[grainNext][GRAIN_POS] = modulo2(delayWritePosF - (800 + param2sq * noise[5] * 390 * (1.25f - param1S)), delayBufferSize);
                float invGrainSize = 1 / grainTable[grainNext][GRAIN_SIZE];
                grainTable[grainNext][GRAIN_CURRENT_SHIFT] = grainTable[grainNext][GRAIN_NEXT_SHIFT];
                grainTable[grainNext][GRAIN_NEXT_SHIFT] = grainRate * invGrainSize;
                grainTable[grainNext][GRAIN_INC] = clamp((grainTable[grainNext][GRAIN_NEXT_SHIFT] - grainTable[grainNext][GRAIN_CURRENT_SHIFT]) * invGrainSize, 0, 0.5f);
                grainTable[grainNext][GRAIN_VOL] = clamp(0.25f + sqrt3(fabsf(noise[3])), 0, 1);
                grainTable[grainNext][GRAIN_PAN] = clamp(1 + (noise[6]) * param2sq, 0, 2) * 0.5f;
                grainPrev = grainNext;

                if (lockB < 0.02f)
                {
                    loopSize = clamp((1 - param1S) * 1400, 48, 1800);
                }
            }
            if (++grainNext > 2)
            {
                grainNext = 0;
            }
        }

        float env;

        const float f2 = 0.85f;
        const float fnotch = 1.03f;

        float grain1, grain2, grain3;
        float grain1L, grain1R;
        float grain2L, grain2R;
        float grain3L, grain3R;
        float grainSumR, grainSumL;

        float * __restrict__ sp = sampleBlock_;
        float * __restrict__ const delayBuf = delayBuffer_;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {

            if (++inputIncCount >= sampleRateDivide)
            {
                inputIncCount = 0;
                delayWritePos = (delayWritePos + 1) & delayBufferSizeM1;
                delayWritePosF = (float)delayWritePos;

                delayReadPos = modulo(delayWritePos + loopSize, delayBufferSize);
                float feedbackIn = delayInterpolation(delayReadPos, delayBuf, delayBufferSizeM1);

                // hp
                float hp_in_x0 = (*sp + *(sp + 1)) * 0.5f;
                hp_in_y0 = _in3_a0 * (hp_in_x0 - hp_in_x1) + _in3_b1 * hp_in_y1;
                hp_in_y1 = hp_in_y0;
                hp_in_x1 = hp_in_x0;

                delayBuf[delayWritePos] = hp_in_y0 * lockA + feedbackIn * lockB;
            }

            ///-------- grain 1
            grain1L = grain1R = 0;
            if (grainTable[0][GRAIN_RAMP] < 1)
            {
                delayReadPos = modulo(grainTable[0][GRAIN_POS] + grainTable[0][GRAIN_RAMP] * grainTable[0][GRAIN_SIZE], delayBufferSize);
                env = hann(grainTable[0][GRAIN_RAMP]) * grainTable[0][GRAIN_VOL];
                grain1 = env * delayInterpolation(delayReadPos, delayBuf, delayBufferSizeM1);
                grain1L = grain1 * grainTable[0][GRAIN_PAN];
                grain1R = grain1 - grain1L;
            }
            grainTable[0][GRAIN_RAMP] += grainTable[0][GRAIN_CURRENT_SHIFT];
            grainTable[0][GRAIN_CURRENT_SHIFT] += grainTable[0][GRAIN_INC];

            ///-------- grain 2
            grain2L = grain2R = 0;
            if (grainTable[1][GRAIN_RAMP] < 1)
            {
                delayReadPos = modulo(grainTable[1][GRAIN_POS] + grainTable[1][GRAIN_RAMP] * grainTable[1][GRAIN_SIZE], delayBufferSize);
                env = hann(grainTable[1][GRAIN_RAMP]) * grainTable[1][GRAIN_VOL];
                grain2 = env * delayInterpolation(delayReadPos, delayBuf, delayBufferSizeM1);
                grain2L = grain2 * grainTable[1][GRAIN_PAN];
                grain2R = grain2 - grain2L;
            }
            grainTable[1][GRAIN_RAMP] += grainTable[1][GRAIN_CURRENT_SHIFT];
            grainTable[1][GRAIN_CURRENT_SHIFT] += grainTable[1][GRAIN_INC];

            ///-------- grain 3
            grain3L = grain3R = 0;
            if (grainTable[2][GRAIN_RAMP] < 1)
            {
                delayReadPos = modulo(grainTable[2][GRAIN_POS] + grainTable[2][GRAIN_RAMP] * grainTable[2][GRAIN_SIZE], delayBufferSize);
                env = hann(grainTable[2][GRAIN_RAMP]) * grainTable[2][GRAIN_VOL];
                grain3 = env * delayInterpolation(delayReadPos, delayBuf, delayBufferSizeM1);
                grain3L = grain3 * grainTable[2][GRAIN_PAN];
                grain3R = grain3 - grain3L;
            }
            grainTable[2][GRAIN_RAMP] += grainTable[2][GRAIN_CURRENT_SHIFT];
            grainTable[2][GRAIN_CURRENT_SHIFT] += grainTable[2][GRAIN_INC];

            grainSumL = grain1L + grain2L + grain3L;
            grainSumR = grain1R + grain2R + grain3R;

            // lp L
            low1 += f2 * band1;
            band1 += f2 * ((grainSumL)-low1 - band1);
            low3 += f2 * band3;
            band3 += f2 * ((low1)-low3 - band3);
            low5 += f2 * band5;
            band5 += f2 * ((low3)-low5 - band5);

            // lp R
            low2 += f2 * band2;
            band2 += f2 * ((grainSumR)-low2 - band2);
            low4 += f2 * band4;
            band4 += f2 * ((low2)-low4 - band4);
            low6 += f2 * band6;
            band6 += f2 * ((low4)-low6 - band6);

            // notch L
            hb1_x1 += fnotch * hb1_y1;
            float high7 = low5 - hb1_x1 - hb1_y1;
            hb1_y1 += fnotch * high7;
            float notchL = (high7 + hb1_x1);

            // notch R
            hb1_x2 += fnotch * hb1_y2;
            float high8 = low6 - hb1_x2 - hb1_y2;
            hb1_y2 += fnotch * high8;
            float notchR = (high8 + hb1_x2);

            *sp = *sp * dry + notchL * wetL;
            sp++;
            *sp = *sp * dry + notchR * wetR;
            sp++;
        }
    }

    break;
    case FILTER2_STEREO_BP:
    {

        mixerGain_ = 0.02f * gainTmp + .98f * mixerGain_;
        float mixerGain_01 = clamp(mixerGain_, 0, 1);
        int mixerGain255 = mixerGain_01 * 255;
        float dry = panTable[255 - mixerGain255];
        float wet = panTable[mixerGain255];
        float extraAmp = clamp(mixerGain_ - 1, 0, 1);
        wet += extraAmp;

        param1S = 0.02f * (this->params_.effect2.param1) + .98f * param1S;

        const float f = param1S * param1S * param1S * 0.9f;
        const float matrixFreqAtnn = matrixFilterFrequency * 0.125f;

        float bpf1 = clamp(0.015f + fold((f + matrixFreqAtnn) * 0.25f) * 3.8f, 0.01f, 0.9f);
        float bpf2 = clamp(0.015f + fold((f - matrixFreqAtnn) * 0.25f) * 3.8f, 0.01f, 0.9f);

        float * __restrict__ sp = sampleBlock_;
        float * __restrict__ const delayBuf = delayBuffer_;

        float filterParam2 = clamp(matrixFilterParam2 + this->params_.effect2.param2, 0, 1) * (1 - param1S * param1S * 0.06f);

        const float fb = sqrt3(0.5f - filterParam2 * 0.497f);
        const float scale = sqrt3(fb);
        const float fb2 = fb * 0.982f;
        const float scale2 = sqrt3(fb2);

        const float inputGain = 1.0f;
        const float finalGain =  3 * (1 - filterParam2 * filterParam2 * 0.5f);

        wet *= finalGain;

        float wetL = wet * (1 + matrixFilterPan);
        float wetR = wet * (1 - matrixFilterPan);

        float high1 = 0;
        float high2 = 0;
        float high3 = 0;
        float high4 = 0;
        float high5 = 0;
        float high6 = 0;
        float high7 = 0;
        float high8 = 0;

        const float f1 = clamp(0.15f + f * 0.5f, 0.01f, 0.99f);
        float coef1 = (1.0f - f1) / (1.0f + f1);

        float drift = _ly1;
        float nexDrift = noise[7] * 0.005f;
        float deltaD = (nexDrift - drift) * 0.000625f;
        _ly1 = nexDrift;

        shift = shift * 0.999f + noise[4] * 0.0001f;

        //   h1 = fondamentale, feedback tanh (résonance non-linéaire, singing quality)
        //   h2 = fondamentale légèrement décalée (~1.5%, drift analogique)
        //   h3 = octave  (×2), excité par la sortie de h1
        //   h4 = tierce  (×3), excité par la sortie de h1
        float jitter   = shift * 0.004f;
        float bpf_h1   = bpf1 * (1.0f + jitter);
        float bpf_h2   = bpf1 * (1.015f - jitter * 0.5f);
        float bpf_h3   = clamp(bpf1 * 2.0f * (1.0f + jitter * 0.3f), 0.01f, 0.92f);
        float bpf_h1_r = bpf2 * (1.0f + jitter);
        float bpf_h2_r = bpf2 * (1.015f - jitter * 0.5f);
        float bpf_h3_r = clamp(bpf2 * 2.0f * (1.0f + jitter * 0.3f), 0.01f, 0.92f);
        float bpf_h4   = clamp(bpf1 * 3.0f * (1.0f - jitter * 0.2f), 0.01f, 0.92f);
        float bpf_h4_r = clamp(bpf2 * 3.0f * (1.0f - jitter * 0.2f), 0.01f, 0.92f);

        // input hp coefs calc :
        const float cutoff = 0.05f;
        const float _in_b1 = (1 - cutoff);
        const float _in_a0 = (1 + _in_b1 * _in_b1 * _in_b1) * 0.5f;
        const float _in_a1 = -_in_a0;

        // limiter
        const float threshold = 0.7f;
        const float kneeWidth = 0.2f;
        const float kneeWidthInv = 1 / (2 * kneeWidth);
        const float threshKneeP = threshold + kneeWidth * 0.5f;
        const float threshKneeM = threshold - kneeWidth * 0.5f;
        const float makeup = 1 / threshold;

        const int delaySize = 256;
        const int delaySizeM1 = delaySize - 1;

        const float attackCoeff = 0.5f;
        const float releaseCoeff = 0.996f;
        const float holdTime = 0.03f;
        const int holdSampleCount = static_cast<int>(holdTime * PREENFM_FREQUENCY);
        // hold counter persistant entre blocs (hb4_y1 inutilisé dans cet effet, pas de RAM en plus)
        int holdSamples = (int)hb4_y1;

        hb4_x1 = clamp(hb4_x1, 0, 1);
        hb4_x2 = clamp(hb4_x2, 0, 1);

        const float lpCoef = 0.7f;

        float target_gain = 1.0f;

        // hb6/hb8 : l'ancien code exécutait 2 fois le même HP sur le même état ; la 2e passe
        // valait exactement y = b1 * y_pass1. Forme fusionnée équivalente (même son, 1 seule passe) :
        const float _in2_a0 = _in_b1 * _in_a0;
        const float _in2_b1 = _in_b1 * _in_b1;

        const float wetLMakeup = wetL * makeup;
        const float wetRMakeup = wetR * makeup;

        for (int k = BLOCK_SIZE; k--;)
        {
            float fbM = fb + drift;
            drift += deltaD;

            // hp input L
            hb5_x1 = (*sp) * inputGain;
            hb5_y1 = _in_a0 * hb5_x1 + _in_a1 * hb5_x2 + _in_b1 * hb5_y2;
            hb5_y2 = hb5_y1;
            hb5_x2 = hb5_x1;

            hb6_x1 = hb5_y1;
            hb6_y1 = _in2_a0 * (hb6_x1 - hb6_x2) + _in2_b1 * hb6_y2;
            hb6_y2 = hb6_y1;
            hb6_x2 = hb6_x1;

            // hp input R
            hb7_x1 = *(sp + 1) * inputGain;
            hb7_y1 = _in_a0 * hb7_x1 + _in_a1 * hb7_x2 + _in_b1 * hb7_y2;
            hb7_y2 = hb7_y1;
            hb7_x2 = hb7_x1;

            hb8_x1 = hb7_y1;
            hb8_y1 = _in2_a0 * (hb8_x1 - hb8_x2) + _in2_b1 * hb8_y2;
            hb8_y2 = hb8_y1;
            hb8_x2 = hb8_x1;

            // Left voice

            hb1_y1 = coef1 * (hb1_y1 + hb6_y1) - hb1_x1; // allpass (shared phase smear)
            hb1_x1 = hb6_y1;

            // h1 — fondamentale, feedback non-linéaire (transistor ladder)
            low1 = low1 + bpf_h1 * band1;
            high1 = scale2 * hb1_y1 - low1 - fb2 * tanh4(band1);
            band1 = bpf_h1 * high1 + band1;

            // h2 — fondamentale légèrement décalée (drift analogique)
            low2 = low2 + bpf_h2 * band2;
            high2 = scale * hb1_y1 - low2 - fbM * band2;
            band2 = bpf_h2 * high2 + band2;

            // h3 — 2e harmonique (octave), en cascade de h1
            low3 = low3 + bpf_h3 * band3;
            high3 = scale2 * band1 * 0.7f - low3 - fb2 * band3;
            band3 = bpf_h3 * high3 + band3;

            // h4 — 3e harmonique (3f), en cascade de h1
            hb2_x1 = hb2_x1 + bpf_h4 * hb2_x2;
            high7 = scale2 * band1 * 0.5f - hb2_x1 - fb2 * hb2_x2;
            hb2_x2 = bpf_h4 * high7 + hb2_x2;

            float outL = band1 + band2 * 0.5f + band3 * 0.7f + hb2_x2 * 0.5f;
            // LP post-SVF tracking BP center
            hb3_y1 += bpf_h1 * (outL - hb3_y1);

            // Right voice

            hb1_y2 = coef1 * (hb1_y2 + hb8_y1) - hb1_x2; // allpass (shared phase smear)
            hb1_x2 = hb8_y1;

            // h1 — fondamentale, feedback non-linéaire (transistor ladder)
            low4 = low4 + bpf_h1_r * band4;
            high4 = scale2 * hb1_y2 - low4 - fb2 * tanh4(band4);
            band4 = bpf_h1_r * high4 + band4;

            // h2 — fondamentale légèrement décalée (drift analogique)
            low5 = low5 + bpf_h2_r * band5;
            high5 = scale * hb1_y2 - low5 - fbM * band5;
            band5 = bpf_h2_r * high5 + band5;

            // h3 — 2e harmonique (octave), en cascade de h1
            low6 = low6 + bpf_h3_r * band6;
            high6 = scale2 * band4 * 0.7f - low6 - fb2 * band6;
            band6 = bpf_h3_r * high6 + band6;

            // h4 — 3e harmonique (3f), en cascade de h1
            hb2_y1 = hb2_y1 + bpf_h4_r * hb2_y2;
            high8 = scale2 * band4 * 0.5f - hb2_y1 - fb2 * hb2_y2;
            hb2_y2 = bpf_h4_r * high8 + hb2_y2;

            float outR = band4 + band5 * 0.5f + band6 * 0.7f + hb2_y2 * 0.5f;
            // LP post-SVF tracking BP center
            hb3_y2 += bpf_h1 * (outR - hb3_y2);

            // limiter delay

            delayWritePos = (delayWritePos + 1) & delaySizeM1;

            _lx1 = _lx1 * lpCoef + tanh4(hb3_y1) * (1.0f - lpCoef);
            _lx2 = _lx2 * lpCoef + tanh4(hb3_y2) * (1.0f - lpCoef);

            delayBuf[delayWritePos] = _lx1;
            delayBuf[delayWritePos + delaySize] = _lx2;

            // limiter — detect on _lx1/_lx2 so envelope matches what is in the delay buffer

            float gain = hb4_x1;
            float envelope = hb4_x2;

            int readpos = (delayWritePos - delaySize) & delaySizeM1;

            float fltOut1 = delayBuf[readpos];
            float fltOut2 = delayBuf[delaySize + readpos];

            float absLeft = fabsf(_lx1);
            float absRight = fabsf(_lx2);
            float absSample = (absLeft > absRight) ? absLeft : absRight;

            envelope = max(absSample, envelope * releaseCoeff);

            target_gain = 1.0f;

            if (envelope > threshKneeP)
            {
                target_gain = threshKneeP / envelope;
                holdSamples = holdSampleCount;
            }
            else if (envelope > threshKneeM)
            {
                // soft knee
                float x = envelope - threshKneeM;
                float y = x * x * kneeWidthInv;
                target_gain = (threshKneeM + y) / envelope;
            }

            // Attaque toujours autorisée ; le hold (30 ms, persistant entre blocs) ne bloque que le relâchement.
            if (target_gain < gain || holdSamples <= 0)
            {
                gain = gain * attackCoeff + target_gain * (1.0f - attackCoeff);
            }
            if (holdSamples > 0)
            {
                holdSamples--;
            }

            hb4_x1 = gain;
            hb4_x2 = envelope;

            //  ------------

            *sp = *sp * dry + fltOut1 * wetLMakeup * gain;
            sp++;

            *sp = *sp * dry + fltOut2 * wetRMakeup * gain;
            sp++;
        }
        hb4_y1 = (float)holdSamples;
    }
    break;
    default:
        // NO EFFECT
        break;
    }

    this->newNotePlayed = false;
}

inline float Timbre::iirFilter(float x, float a1, float *xn1, float *xn2, float *yn1, float *yn2)
{
    // https://dsp.stackexchange.com/a/59157
    // 𝑦[𝑘] = 𝑐( 𝑥[𝑘] + 𝑦[𝑘−2] ) − 𝑥[𝑘−2]
    float y = a1 * (x + *yn2) - *xn2;
    *yn2 = *yn1;
    *yn1 = y;
    *xn2 = *xn1;
    *xn1 = x;

    return y;
}

inline float Timbre::delayInterpolation(float readPos, float buffer[], int bufferLenM1)
{
    int readPosInt = readPos;
    float y1 = buffer[readPosInt];
    float y0 = buffer[(readPosInt - 1) & bufferLenM1];
    float x = 1 - (readPos - floorf(readPos));
    return (y0 - y1) * x + y1;
}

inline float Timbre::delayInterpolation2(float readPos, float buffer[], int bufferLenM1, int offset)
{
    int readPosInt = readPos;
    float y1 = buffer[offset + readPosInt];
    float y0 = buffer[offset + ((readPosInt - 1) & bufferLenM1)];
    float x = 1 - (readPos - floorf(readPos));
    return (y0 - y1) * x + y1;
}

inline float Timbre::delayInterpolationHermite(float readPos, const float buffer[], int bufferLenM1, int offset)
{
    int readPosInt = (int)readPos;
    float frac = readPos - (float)readPosInt;
    float xm1 = buffer[offset + ((readPosInt - 1) & bufferLenM1)];
    float x0  = buffer[offset + (readPosInt & bufferLenM1)];
    float x1  = buffer[offset + ((readPosInt + 1) & bufferLenM1)];
    float x2  = buffer[offset + ((readPosInt + 2) & bufferLenM1)];
    return hermiteInterpolation(frac, xm1, x0, x1, x2);
}

inline float Timbre::hermiteInterpolation(float frac, float xm1, float x0, float x1, float x2)
{
    float c0 = x0;
    float c1 = 0.5f * (x1 - xm1);
    float c2 = xm1 - 2.5f * x0 + 2.0f * x1 - 0.5f * x2;
    float c3 = 0.5f * (x2 - xm1) + 1.5f * (x0 - x1);

    return ((c3 * frac + c2) * frac + c1) * frac + c0;
}

#if defined(__GNUC__) && !defined(DEBUG)
#pragma GCC pop_options
#endif
