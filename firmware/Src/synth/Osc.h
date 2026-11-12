/*
 * Copyright 2013 Xavier Hosxe
 *
 * Author: Xavier Hosxe (xavier . hosxe (at) gmail . com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <math.h>

#include "SynthStateAware.h"
#include "Matrix.h"

#ifndef PFM3_FAST_FEEDBACK_INTERP
#define PFM3_FAST_FEEDBACK_INTERP 0
#endif

extern float sinTable[];
struct OscState {
    // Current wavetable phase/index (wrapped to table size each sample/block).
    float index;
    // Per-sample oscillator frequency used by render paths.
    float frequency;
    // Base frequency after harmonic multiplier + matrix pitch modulation.
    float mainFrequencyPlusMatrix;
    // Base note frequency before matrix pitch modulation.
    float mainFrequency;
    float fromFrequency;
    float nextFrequency;

    // Output quantization depth (1..19 in decimation modes).
    uint8_t waveDecimationBits;
    // Enables sample-rate and bit-depth decimation path when non-zero.
    uint8_t waveDecimationEnabled;
    // Enables HQ interpolation path when non-zero (ignored if decimation is enabled).
    uint8_t waveInterpolationEnabled;
    // Keeps D=2 sample-and-hold phase continuity across block boundaries.
    uint8_t waveDecimationStepPhase;
    // Precomputed scale factors for quantization/dequantization.
    float waveDecimationScale;
    float waveDecimationInvScale;
    // Reused as held sample in decimation and previous sample in HQ interpolation.
    float waveDecimationHeldSample;
    // Warp value after matrix modulation and clamping.
    float effectiveWarp;
};


extern struct WaveTable waveTables[];
extern float* waveMipmapTables[][5];
extern float exp2_harm[];


class Osc : public SynthStateAware
{
public:
    Osc() {};
    virtual ~Osc() {};

    void init(SynthState* synthState, struct OscillatorParams *oscParams, struct OperatorPhaseRowParams* phaseParamsBase, DestinationEnum df);

    void newNote(struct OscState* oscState, float newNoteFrequency, float phase);
    float getNoteRealFrequencyEstimation(struct OscState* oscState, float newNoteFrequency);
    void glideToNote(struct OscState* oscState, float newNoteFrequency);
    inline void glideStep(struct OscState* oscState, float phase) {
        oscState->mainFrequency = oscState->fromFrequency + (oscState->nextFrequency - oscState->fromFrequency) * phase;
    }

    inline __attribute__((always_inline)) void updateWarpWithMatrix(struct OscState *oscState, Matrix* matrix) {
        float warp = *phaseWarpParam;
        if (destWarp != DESTINATION_NONE) {
            warp += matrix->getDestination(destWarp);
        }

        if (!(warp == warp)) {
            warp = 0.0f;
        } else if (warp > 4.0f) {
            warp = 4.0f;
        } else if (warp < -4.0f) {
            warp = -4.0f;
        }

        // Tiny deadband keeps the no-warp fast path stable under near-zero modulation.
        if (warp > -0.0001f && warp < 0.0001f) {
            warp = 0.0f;
        }
        oscState->effectiveWarp = warp;
    }

    inline __attribute__((always_inline)) void calculateFrequencyWithMatrix(struct OscState *oscState, Matrix* matrix, float expHarm) {
        oscState->mainFrequencyPlusMatrix = oscState->mainFrequency;
        oscState->mainFrequencyPlusMatrix *= expHarm;
        oscState->mainFrequencyPlusMatrix +=  (oscState->mainFrequency  * (matrix->getDestination(destFreq) + matrix->getDestination(ALL_OSC_FREQ)) * .1f);
    }

    // Quantize a raw wavetable sample to the active decimation bit depth.
    inline __attribute__((always_inline)) float quantizeWaveSample(struct OscState *oscState, float value) {
        float scaled = value * oscState->waveDecimationScale;
        int q = (int)scaled;
        return (float)q * oscState->waveDecimationInvScale;
    }

    inline __attribute__((always_inline)) float quantizeOscOutputBeforeEnvelope(struct OscState *oscState, float outputSample) {
        if (!oscState->waveDecimationEnabled) {
            return outputSample;
        }
        return quantizeWaveSample(oscState, outputSample);
    }

    // Fast piecewise-linear half-cycle warp remap in index space (used by decimation paths).
    inline __attribute__((always_inline)) int getWarpedIndexFast(int iIndex, int max, int halfSize, float slopeFirstHalf, float slopeSecondHalf, float secondHalfOffset) {
        float fi = (float)iIndex;
        float warped = iIndex < halfSize
                ? fi * slopeFirstHalf
                : fi * slopeSecondHalf + secondHalfOffset;
        int indexInteger = (int)warped;
        return indexInteger & max;
    }

    // Linear wavetable lookup for arbitrary (possibly negative) floating-point index.
    inline __attribute__((always_inline)) float sampleLinear(const float* wave, float fIndex, int max) {
        int i0 = (int)fIndex;
        float frac = fIndex - (float)i0;
        if (unlikely(frac < 0.0f)) {
            frac += 1.0f;
            i0 -= 1;
        }
        i0 &= max;
        float s0 = wave[i0];
        float s1 = wave[(i0 + 1) & max];
        return s0 + frac * (s1 - s0);
    }

    // Branchless linear wavetable lookup when fIndex >= 0.0f is guaranteed.
    inline __attribute__((always_inline)) float sampleLinearPositive(const float* wave, float fIndex, int max) {
        int i0 = (int)fIndex;
        float frac = fIndex - (float)i0;
        i0 &= max;
        float s0 = wave[i0];
        float s1 = wave[(i0 + 1) & max];
        return s0 + frac * (s1 - s0);
    }

    // Continuous-phase piecewise-linear half-cycle warp with linear wavetable interpolation.
    inline __attribute__((always_inline)) float sampleWarpedLinear(const float* wave, float fPhaseInTable, int max, float halfSize, float slopeFirstHalf, float slopeSecondHalf, float secondHalfOffset) {
        float warped = fPhaseInTable < halfSize
                ? fPhaseInTable * slopeFirstHalf
                : fPhaseInTable * slopeSecondHalf + secondHalfOffset;
        return sampleLinear(wave, warped, max);
    }

    // Select band-limited mipmap table pointer for the given shape and oscillator frequency (Hz).
    inline __attribute__((always_inline)) float* getWaveTableMipmap(int shape, float freqHz) {
        float absFreq = fabsf(freqHz);
        if (likely(absFreq <= 1200.0f)) {
            return waveTables[shape].table;
        }
        int level;
        if (absFreq <= 2400.0f) {
            level = 1;
        } else if (absFreq <= 4800.0f) {
            level = 2;
        } else if (absFreq <= 9600.0f) {
            level = 3;
        } else {
            level = 4;
        }
        return waveMipmapTables[shape][level];
    }

    // Single-sample renderer used by non-block paths.
    // Mode behavior:
    // - Decimation enabled: D=2 sample-and-hold + bit-depth quantization.
    // - Decimation disabled, interpolation disabled: linear wavetable interpolation (Full).
    // - Decimation disabled, interpolation enabled: linear wavetable interpolation + 2-point smoothing (HQ).
    inline __attribute__((always_inline)) float getNextSample(struct OscState *oscState)  {
        int shape = (int) oscillator->shape;
        struct WaveTable* waveTable = &waveTables[shape];
        int max = waveTable->max;
        float* wave = getWaveTableMipmap(shape, oscState->frequency);
        float phaseIncrement = oscState->frequency * waveTable->precomputedValue + waveTable->floatToAdd;
        float warp = oscState->effectiveWarp;
        bool waveDecimationEnabled = oscState->waveDecimationEnabled != 0;
        bool waveInterpolationEnabled = (oscState->waveInterpolationEnabled != 0) && !waveDecimationEnabled;

        if (likely(warp == 0.0f)) {
            if (unlikely(waveDecimationEnabled)) {
                if (oscState->waveDecimationStepPhase != 0) {
                    oscState->waveDecimationStepPhase = 0;
                    return oscState->waveDecimationHeldSample;
                }

                oscState->waveDecimationStepPhase = 1;
                float phaseIncrement2 = phaseIncrement + phaseIncrement;
                oscState->index += phaseIncrement2;

                int indexInteger = oscState->index;
                oscState->index -= indexInteger;
                indexInteger &= max;
                oscState->index += indexInteger;

                float sample = wave[indexInteger];
                sample = quantizeWaveSample(oscState, sample);
                oscState->waveDecimationHeldSample = sample;
                return sample;
            }

            float fIndexMid = oscState->index + 0.5f * phaseIncrement;

            float fIndex = oscState->index + phaseIncrement;
            int indexInteger = (int)fIndex;
            float frac = fIndex - (float)indexInteger;
            if (unlikely(frac < 0.0f)) {
                frac += 1.0f;
                indexInteger -= 1;
            }
            indexInteger &= max;
            oscState->index = (float)indexInteger + frac;

            float s0 = wave[indexInteger];
            float s1 = wave[(indexInteger + 1) & max];
            float currentValue = s0 + frac * (s1 - s0);

            if (unlikely(waveInterpolationEnabled)) {
                float prevValue = oscState->waveDecimationHeldSample;
                oscState->waveDecimationHeldSample = currentValue;
                float midValue = sampleLinear(wave, fIndexMid, max);
                return 0.5f * midValue + 0.25f * (prevValue + currentValue);
            }

            return currentValue;
        }

        int size = max + 1;
        int halfSize = size >> 1;
        float fHalfSize = (float)halfSize;
        float slopeFirstHalf = 1.0f + warp;
        float slopeSecondHalf = 1.0f - warp;
        float secondHalfOffset = ((float)size) * warp;

        if (unlikely(waveDecimationEnabled)) {
            if (oscState->waveDecimationStepPhase != 0) {
                oscState->waveDecimationStepPhase = 0;
                return oscState->waveDecimationHeldSample;
            }

            oscState->waveDecimationStepPhase = 1;
            float phaseIncrement2 = phaseIncrement + phaseIncrement;
            oscState->index += phaseIncrement2;

            int indexInteger = oscState->index;
            oscState->index -= indexInteger;
            indexInteger &= max;
            oscState->index += indexInteger;

            indexInteger = getWarpedIndexFast(indexInteger, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);

            float sample = wave[indexInteger];
            sample = quantizeWaveSample(oscState, sample);
            oscState->waveDecimationHeldSample = sample;
            return sample;
        }

        float fIndexMid = oscState->index + 0.5f * phaseIncrement;

        float fIndex = oscState->index + phaseIncrement;
        int indexInteger = (int)fIndex;
        float frac = fIndex - (float)indexInteger;
        if (unlikely(frac < 0.0f)) {
            frac += 1.0f;
            indexInteger -= 1;
        }
        indexInteger &= max;
        float fPhaseInTable = (float)indexInteger + frac;
        oscState->index = fPhaseInTable;

        float currentValue = sampleWarpedLinear(wave, fPhaseInTable, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);

        if (unlikely(waveInterpolationEnabled)) {
            int idxMid = (int)fIndexMid;
            float fracMid = fIndexMid - (float)idxMid;
            if (unlikely(fracMid < 0.0f)) {
                fracMid += 1.0f;
                idxMid -= 1;
            }
            float fPhaseMid = (float)(idxMid & max) + fracMid;
            float midValue = sampleWarpedLinear(wave, fPhaseMid, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            float prevValue = oscState->waveDecimationHeldSample;
            oscState->waveDecimationHeldSample = currentValue;
            return 0.5f * midValue + 0.25f * (prevValue + currentValue);
        }

        return currentValue;
    }

    inline __attribute__((always_inline)) float getPhase(struct OscState *oscState)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        return oscState->index * waveTable->phaseMul;
    }

    inline __attribute__((always_inline)) float geIndexFromtPhase(float phase)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        return phase * waveTable->max;
    }

    // Phase increment per sample (wavetable index units) for the current oscState->frequency.
    inline __attribute__((always_inline)) float getPhaseIncrement(struct OscState *oscState)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        return oscState->frequency * waveTable->precomputedValue + waveTable->floatToAdd;
    }

    // Sub-sample accurate hard sync of this (slave) oscillator on a master oscillator wrap.
    // Must be called after the master getNextSample() that wrapped (masterPhase = elapsed phase
    // since the wrap) and after oscState->frequency of this slave is set for the current sample,
    // and before this slave getNextSample().
    // The slave restarts at phase 0 at the exact wrap instant, so it only advances by the
    // fraction of the sample elapsed since the wrap instead of being rounded to a whole sample.
    inline __attribute__((always_inline)) void hardSyncToMaster(struct OscState *oscState, Osc *masterOsc, struct OscState *masterState, float masterPhase)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        float masterIncrement = masterOsc->getPhaseIncrement(masterState) * waveTables[(int) masterOsc->oscillator->shape].phaseMul;

        if (unlikely(oscState->waveDecimationEnabled != 0 || masterState->waveDecimationEnabled != 0 || !(masterIncrement > 0.0f))) {
            // Decimated or reversed phase: keep the legacy sample-rounded reset.
            oscState->index = masterPhase * waveTable->max;
            return;
        }

        // Fraction of the current sample that elapsed since the master wrap (0..1).
        float elapsed = masterPhase / masterIncrement;
        if (elapsed > 1.0f) {
            elapsed = 1.0f;
        }
        // getNextSample() will add one full increment: pre-compensate so that the slave
        // phase after that call is exactly elapsed * increment. Offset by one table length
        // to keep the index positive.
        float slaveIncrement = getPhaseIncrement(oscState);
        oscState->index = (float) (waveTable->max + 1) + (elapsed - 1.0f) * slaveIncrement;
    }

    // Block renderer for Full mode. If HQ interpolation is requested, this function
    // delegates to getNextBlockHQ so the HQ path stays centralized.
    inline __attribute__((always_inline)) float* getNextBlock(struct OscState *oscState)  {

        if (unlikely(oscState->waveInterpolationEnabled != 0)) {
            return getNextBlockHQ(oscState);
        }

        int shape = (int) oscillator->shape;
        int max = waveTables[shape].max;
        float *wave = getWaveTableMipmap(shape, oscState->frequency);
        float freq = oscState->frequency * waveTables[shape].precomputedValue + waveTables[shape].floatToAdd;
        float warp = oscState->effectiveWarp;
        bool warpEnabled = warp != 0.0f;
        bool waveDecimationEnabled = oscState->waveDecimationEnabled != 0;
        float fIndex = oscState->index;
        if (unlikely(freq < 0.0f)) {
            int shift = ((int)(-32.0f * freq) | max) + 1;
            fIndex += (float)shift;
        }
        int iIndex;
        float* oscValuesToFill = oscValues[oscValuesCpt];
        oscValuesCpt++;
        oscValuesCpt &= 0x3;

        if (unlikely(waveDecimationEnabled)) {
            float freq2 = freq + freq;
            int k = 0;
            if (oscState->waveDecimationStepPhase != 0) {
                oscValuesToFill[k++] = oscState->waveDecimationHeldSample;
                oscState->waveDecimationStepPhase = 0;
            }

            if (warpEnabled) {
                int size = max + 1;
                int halfSize = size >> 1;
                float slopeFirstHalf = 1.0f + warp;
                float slopeSecondHalf = 1.0f - warp;
                float secondHalfOffset = ((float)size) * warp;
                for (; k + 1 < BLOCK_SIZE; k += 2) {
                    fIndex += freq2;
                    iIndex = ((int)fIndex) & max;
                    iIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
                    float sample = quantizeWaveSample(oscState, wave[iIndex]);
                    oscState->waveDecimationHeldSample = sample;
                    oscValuesToFill[k] = sample;
                    oscValuesToFill[k + 1] = sample;
                }
                if (k < BLOCK_SIZE) {
                    fIndex += freq2;
                    iIndex = ((int)fIndex) & max;
                    iIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
                    float sample = quantizeWaveSample(oscState, wave[iIndex]);
                    oscState->waveDecimationHeldSample = sample;
                    oscValuesToFill[k] = sample;
                    oscState->waveDecimationStepPhase = 1;
                }
            } else {
                for (; k + 1 < BLOCK_SIZE; k += 2) {
                    fIndex += freq2;
                    iIndex = ((int)fIndex) & max;
                    float sample = quantizeWaveSample(oscState, wave[iIndex]);
                    oscState->waveDecimationHeldSample = sample;
                    oscValuesToFill[k] = sample;
                    oscValuesToFill[k + 1] = sample;
                }
                if (k < BLOCK_SIZE) {
                    fIndex += freq2;
                    iIndex = ((int)fIndex) & max;
                    float sample = quantizeWaveSample(oscState, wave[iIndex]);
                    oscState->waveDecimationHeldSample = sample;
                    oscValuesToFill[k] = sample;
                    oscState->waveDecimationStepPhase = 1;
                }
            }
            iIndex = (int)fIndex;
            oscState->index = (fIndex - iIndex) + (iIndex & max);
            return oscValuesToFill;
        }

        if (likely(!warpEnabled)) {
            for (int k = 0; k < 32; ) {
                fIndex += freq;
                oscValuesToFill[k++] = sampleLinearPositive(wave, fIndex, max);

                fIndex += freq;
                oscValuesToFill[k++] = sampleLinearPositive(wave, fIndex, max);

                fIndex += freq;
                oscValuesToFill[k++] = sampleLinearPositive(wave, fIndex, max);

                fIndex += freq;
                oscValuesToFill[k++] = sampleLinearPositive(wave, fIndex, max);
            }

            iIndex = (int)fIndex;
            oscState->index = (fIndex - iIndex) + (iIndex & max);
            return oscValuesToFill;
        }

        int size = max + 1;
        float fHalfSize = (float)(size >> 1);
        float slopeFirstHalf = 1.0f + warp;
        float slopeSecondHalf = 1.0f - warp;
        float secondHalfOffset = ((float)size) * warp;

        for (int k = 0; k < 32; ) {
            fIndex += freq;
            iIndex = (int)fIndex;
            float fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
            oscValuesToFill[k++] = sampleWarpedLinear(wave, fPhase, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);

            fIndex += freq;
            iIndex = (int)fIndex;
            fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
            oscValuesToFill[k++] = sampleWarpedLinear(wave, fPhase, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);

            fIndex += freq;
            iIndex = (int)fIndex;
            fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
            oscValuesToFill[k++] = sampleWarpedLinear(wave, fPhase, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);

            fIndex += freq;
            iIndex = (int)fIndex;
            fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
            oscValuesToFill[k++] = sampleWarpedLinear(wave, fPhase, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
        }

        iIndex = (int)fIndex;
        oscState->index = (fIndex - iIndex) + (iIndex & max);
        return oscValuesToFill;
    };

    // Feedback-capable block renderer used by feedback operator paths.
    // Decimation output quantization is preserved; interpolation is optional and can
    // be force-disabled with PFM3_FAST_FEEDBACK_INTERP for speed-focused builds.
    inline __attribute__((always_inline)) float* getNextBlockWithFeedbackAndEnveloppe(struct OscState *oscState, float feedback, float& env, float envInc, float freqMultiplier, float* lastValue) {
        int shape = (int) oscillator->shape;
        int max = waveTables[shape].max;
        float *wave = getWaveTableMipmap(shape, oscState->frequency);
        float freq = oscState->frequency * waveTables[shape].precomputedValue + waveTables[shape].floatToAdd;
        float fIndex = oscState->index;
        if (unlikely(freq < 0.0f)) {
            int shift = ((int)(-32.0f * freq) | max) + 1;
            fIndex += (float)shift;
        }
        int iIndex;
        float* oscValuesToFill = oscValues[4];

        float warp = oscState->effectiveWarp;
        bool warpEnabled = warp != 0.0f;
        bool decimationEnabled = oscState->waveDecimationEnabled != 0;
    #if PFM3_FAST_FEEDBACK_INTERP
        // Speed-first mode: skip feedback-path interpolation to reduce per-sample math.
        bool interpolationEnabled = false;
    #else
        bool interpolationEnabled = (oscState->waveInterpolationEnabled != 0) && !decimationEnabled;
    #endif

        lastValue[2] = .95f * lastValue[2] + feedback * .05f;
        float phaseModulationAmplitude = lastValue[2] * ((float) max) * .5f;

        float localLastValue0 = lastValue[0];
        float localLastValue1 = lastValue[1];

        // Optimisation to avoid multiple freqMultiplier in the loop
        float localEnvM = env * freqMultiplier;
        float envIncM   = envInc   * freqMultiplier;
        float halfFreq  = 0.5f * freq;
        if (likely(!warpEnabled)) {
            if (interpolationEnabled) {
                for (int k = 0; k < 32; k++) {
                    fIndex += halfFreq;
                    float s1 = sampleLinear(wave, fIndex + localLastValue0 * phaseModulationAmplitude, max);
                    float midLast0 = s1 - localLastValue1 + .99525f * localLastValue0;

                    fIndex += halfFreq;
                    float s2 = sampleLinear(wave, fIndex + midLast0 * phaseModulationAmplitude, max);

                    float newValue = 0.5f * s1 + 0.25f * (localLastValue1 + s2);
                    localLastValue0 = newValue - localLastValue1 + .99525f * localLastValue0;
                    localLastValue1 = newValue;

                    oscValuesToFill[k] = localLastValue0 * localEnvM;
                    localEnvM += envIncM;
                }
            } else if (unlikely(decimationEnabled)) {
                float decimationScale = oscState->waveDecimationScale;
                localEnvM *= oscState->waveDecimationInvScale;
                envIncM *= oscState->waveDecimationInvScale;
                for (int k = 0; k < 32; k++) {
                    fIndex += freq;
                    iIndex = ((int)fIndex) & max;

                    float phaseModulationOffset = localLastValue0 * phaseModulationAmplitude;
                    int index = (iIndex + (int)phaseModulationOffset) & max;

                    float newValue = wave[index];
                    localLastValue0 = newValue - localLastValue1 + .99525f * localLastValue0;
                    localLastValue1 = newValue;

                    int q = (int)(localLastValue0 * decimationScale);
                    oscValuesToFill[k] = (float)q * localEnvM;
                    localEnvM += envIncM;
                }
            } else {
                for (int k = 0; k < 32; k++) {
                    fIndex += freq;
                    float modIndex = fIndex + localLastValue0 * phaseModulationAmplitude;
                    float newValue = sampleLinear(wave, modIndex, max);
                    localLastValue0 = newValue - localLastValue1 + .99525f * localLastValue0;
                    localLastValue1 = newValue;

                    oscValuesToFill[k] = localLastValue0 * localEnvM;
                    localEnvM += envIncM;
                }
            }
        } else {
            int size = max + 1;
            int halfSize = size >> 1;
            float fHalfSize = (float)halfSize;
            float slopeFirstHalf = 1.0f + warp;
            float slopeSecondHalf = 1.0f - warp;
            float secondHalfOffset = ((float)size) * warp;
            if (interpolationEnabled) {
                for (int k = 0; k < 32; k++) {
                    fIndex += halfFreq;
                    iIndex = (int)fIndex;
                    float fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
                    float warpedPhase = fPhase < fHalfSize
                            ? fPhase * slopeFirstHalf
                            : fPhase * slopeSecondHalf + secondHalfOffset;
                    float s1 = sampleLinear(wave, warpedPhase + localLastValue0 * phaseModulationAmplitude, max);
                    float midLast0 = s1 - localLastValue1 + .99525f * localLastValue0;

                    fIndex += halfFreq;
                    iIndex = (int)fIndex;
                    fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
                    warpedPhase = fPhase < fHalfSize
                            ? fPhase * slopeFirstHalf
                            : fPhase * slopeSecondHalf + secondHalfOffset;
                    float s2 = sampleLinear(wave, warpedPhase + midLast0 * phaseModulationAmplitude, max);

                    float newValue = 0.5f * s1 + 0.25f * (localLastValue1 + s2);
                    localLastValue0 = newValue - localLastValue1 + .99525f * localLastValue0;
                    localLastValue1 = newValue;

                    oscValuesToFill[k] = localLastValue0 * localEnvM;
                    localEnvM += envIncM;
                }
            } else if (unlikely(decimationEnabled)) {
                float decimationScale = oscState->waveDecimationScale;
                localEnvM *= oscState->waveDecimationInvScale;
                envIncM *= oscState->waveDecimationInvScale;
                for (int k = 0; k < 32; k++) {
                    fIndex += freq;
                    iIndex = ((int)fIndex) & max;

                    int warpedIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);

                    float phaseModulationOffset = localLastValue0 * phaseModulationAmplitude;
                    int index = (warpedIndex + (int)phaseModulationOffset) & max;

                    float newValue = wave[index];
                    localLastValue0 = newValue - localLastValue1 + .99525f * localLastValue0;
                    localLastValue1 = newValue;

                    int q = (int)(localLastValue0 * decimationScale);
                    oscValuesToFill[k] = (float)q * localEnvM;
                    localEnvM += envIncM;
                }
            } else {
                for (int k = 0; k < 32; k++) {
                    fIndex += freq;
                    iIndex = (int)fIndex;
                    float fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
                    float warpedPhase = fPhase < fHalfSize
                            ? fPhase * slopeFirstHalf
                            : fPhase * slopeSecondHalf + secondHalfOffset;

                    float modIndex = warpedPhase + localLastValue0 * phaseModulationAmplitude;
                    float newValue = sampleLinear(wave, modIndex, max);
                    localLastValue0 = newValue - localLastValue1 + .99525f * localLastValue0;
                    localLastValue1 = newValue;

                    oscValuesToFill[k] = localLastValue0 * localEnvM;
                    localEnvM += envIncM;
                }
            }
        }
        lastValue[0] = localLastValue0;
        lastValue[1] = localLastValue1;

        // update env
        env += envInc * 32;

        iIndex = (int)fIndex;
        oscState->index = (fIndex - iIndex) + (iIndex & max);
        return oscValuesToFill;
    }

    // HQ block renderer. Called when waveInterpolationEnabled != 0 (decimation is off).
    inline __attribute__((always_inline)) float* getNextBlockHQ(struct OscState *oscState)  {
        int shape = (int) oscillator->shape;
        int max = waveTables[shape].max;
        float *wave = getWaveTableMipmap(shape, oscState->frequency);
        float freq = oscState->frequency * waveTables[shape].precomputedValue + waveTables[shape].floatToAdd;
        float warp = oscState->effectiveWarp;
        bool warpEnabled = warp != 0.0f;
        float fIndex = oscState->index;
        if (unlikely(freq < 0.0f)) {
            int shift = ((int)(-32.0f * freq) | max) + 1;
            fIndex += (float)shift;
        }
        int iIndex;
        float* oscValuesToFill = oscValues[oscValuesCpt];
        oscValuesCpt++;
        oscValuesCpt &= 0x3;

        float halfFreq = 0.5f * freq;
        if (likely(!warpEnabled)) {
            float previousValue = oscState->waveDecimationHeldSample;
            for (int k = 0; k < 32; ) {
                fIndex += halfFreq;
                float midValue = sampleLinearPositive(wave, fIndex, max);
                fIndex += halfFreq;
                float currentValue = sampleLinearPositive(wave, fIndex, max);
                oscValuesToFill[k++] = 0.5f * midValue + 0.25f * (previousValue + currentValue);
                previousValue = currentValue;

                fIndex += halfFreq;
                midValue = sampleLinearPositive(wave, fIndex, max);
                fIndex += halfFreq;
                currentValue = sampleLinearPositive(wave, fIndex, max);
                oscValuesToFill[k++] = 0.5f * midValue + 0.25f * (previousValue + currentValue);
                previousValue = currentValue;

                fIndex += halfFreq;
                midValue = sampleLinearPositive(wave, fIndex, max);
                fIndex += halfFreq;
                currentValue = sampleLinearPositive(wave, fIndex, max);
                oscValuesToFill[k++] = 0.5f * midValue + 0.25f * (previousValue + currentValue);
                previousValue = currentValue;

                fIndex += halfFreq;
                midValue = sampleLinearPositive(wave, fIndex, max);
                fIndex += halfFreq;
                currentValue = sampleLinearPositive(wave, fIndex, max);
                oscValuesToFill[k++] = 0.5f * midValue + 0.25f * (previousValue + currentValue);
                previousValue = currentValue;
            }

            oscState->waveDecimationHeldSample = previousValue;
            iIndex = (int)fIndex;
            oscState->index = (fIndex - iIndex) + (iIndex & max);
            return oscValuesToFill;
        }

        int size = max + 1;
        float fHalfSize = (float)(size >> 1);
        float slopeFirstHalf = 1.0f + warp;
        float slopeSecondHalf = 1.0f - warp;
        float secondHalfOffset = ((float)size) * warp;

        float previousValue = oscState->waveDecimationHeldSample;
        for (int k = 0; k < 32; ) {
            fIndex += halfFreq;
            iIndex = (int)fIndex;
            float fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
            float midValue = sampleWarpedLinear(wave, fPhase, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            fIndex += halfFreq;
            iIndex = (int)fIndex;
            fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
            float currentValue = sampleWarpedLinear(wave, fPhase, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            oscValuesToFill[k++] = 0.5f * midValue + 0.25f * (previousValue + currentValue);
            previousValue = currentValue;

            fIndex += halfFreq;
            iIndex = (int)fIndex;
            fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
            midValue = sampleWarpedLinear(wave, fPhase, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            fIndex += halfFreq;
            iIndex = (int)fIndex;
            fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
            currentValue = sampleWarpedLinear(wave, fPhase, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            oscValuesToFill[k++] = 0.5f * midValue + 0.25f * (previousValue + currentValue);
            previousValue = currentValue;

            fIndex += halfFreq;
            iIndex = (int)fIndex;
            fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
            midValue = sampleWarpedLinear(wave, fPhase, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            fIndex += halfFreq;
            iIndex = (int)fIndex;
            fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
            currentValue = sampleWarpedLinear(wave, fPhase, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            oscValuesToFill[k++] = 0.5f * midValue + 0.25f * (previousValue + currentValue);
            previousValue = currentValue;

            fIndex += halfFreq;
            iIndex = (int)fIndex;
            fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
            midValue = sampleWarpedLinear(wave, fPhase, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            fIndex += halfFreq;
            iIndex = (int)fIndex;
            fPhase = (float)(iIndex & max) + (fIndex - (float)iIndex);
            currentValue = sampleWarpedLinear(wave, fPhase, max, fHalfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            oscValuesToFill[k++] = 0.5f * midValue + 0.25f * (previousValue + currentValue);
            previousValue = currentValue;
        }

        oscState->waveDecimationHeldSample = previousValue;
        iIndex = (int)fIndex;
        oscState->index = (fIndex - iIndex) + (iIndex & max);
        return oscValuesToFill;
    };


private:
    DestinationEnum destFreq;
    DestinationEnum destWarp;
    static float* oscValues[5];
    static int oscValuesCpt;
    OscillatorParams* oscillator;
    float* phaseWarpParam;
    float phaseWarpFallback;
};
