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
    // Keeps D=2 sample-and-hold phase continuity across block boundaries.
    uint8_t waveDecimationStepPhase;
    // Precomputed scale factors for quantization/dequantization.
    float waveDecimationScale;
    float waveDecimationInvScale;
    // Held sample across D=2 decimation steps and block boundaries.
    float waveDecimationHeldSample;
    // Warp value after matrix modulation and clamping.
    float effectiveWarp;
    // Phase rate multiplier (1.0f at 1x, 0.5f when Voice runs 2x oversampling).
    float phaseRateScale;
    // Non-zero when Voice runs 2x oversampling (adapts feedback DC-blocker pole).
    uint8_t oversample2x;
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

    // C1-continuous periodic half-cycle phase warp in table-index space.
    // Preserves peak instantaneous slope (1 + |warp|) without derivative corners at 0 or halfSize.
    inline __attribute__((always_inline)) float warpPhaseSmooth(float fPhaseInTable, float halfSize, float invHalfSize, float warpScale) {
        float t;
        float delta;
        if (fPhaseInTable < halfSize) {
            t = fPhaseInTable * invHalfSize;
            delta = t * t * (1.5f - t);
        } else {
            t = (fPhaseInTable - halfSize) * invHalfSize;
            delta = 0.5f - t * t * (1.5f - t);
        }
        return fPhaseInTable + warpScale * delta;
    }

    // C1-continuous half-cycle phase warp with linear wavetable interpolation.
    inline __attribute__((always_inline)) float sampleWarpedLinear(const float* wave, float fPhaseInTable, int max, float halfSize, float invHalfSize, float warpScale) {
        float warped = warpPhaseSmooth(fPhaseInTable, halfSize, invHalfSize, warpScale);
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
        float warp = oscState->effectiveWarp;
        float mipmapFreq = (warp != 0.0f) ? oscState->frequency * (1.0f + fabsf(warp)) : oscState->frequency;
        float* wave = getWaveTableMipmap(shape, mipmapFreq);
        float phaseIncrement = (oscState->frequency * waveTable->precomputedValue + waveTable->floatToAdd) * oscState->phaseRateScale;
        bool waveDecimationEnabled = oscState->waveDecimationEnabled != 0;

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
            return s0 + frac * (s1 - s0);
        }

        int size = max + 1;
        int halfSize = size >> 1;
        float fHalfSize = (float)halfSize;
        float invHalfSize = 2.0f * waveTable->phaseMul;
        float warpScale = warp * (2.0f / 3.0f) * ((float)size);
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

        return sampleWarpedLinear(wave, fPhaseInTable, max, fHalfSize, invHalfSize, warpScale);
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
        return (oscState->frequency * waveTable->precomputedValue + waveTable->floatToAdd) * oscState->phaseRateScale;
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

    // Block renderer for Full / decimated modes.
    // Block renderer for Full / decimated modes (out-of-line in Osc.cpp to keep Voice::runAlgoBlock compact in I-Cache).
    float* getNextBlock(struct OscState *oscState);

    // Feedback-capable block renderer used by feedback operator paths (out-of-line in Osc.cpp).
    float* getNextBlockWithFeedbackAndEnveloppe(struct OscState *oscState, float feedback, float& env, float envInc, float freqMultiplier, float* lastValue);


private:
    DestinationEnum destFreq;
    DestinationEnum destWarp;
    static float* oscValues[5];
    static int oscValuesCpt;
    OscillatorParams* oscillator;
    float* phaseWarpParam;
    float phaseWarpFallback;
};
