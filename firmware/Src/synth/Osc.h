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
    void glideStep(struct OscState* oscState, float phase);

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

    // Fast piecewise-linear half-cycle warp remap in index space.
    inline __attribute__((always_inline)) int getWarpedIndexFast(int iIndex, int max, int halfSize, float slopeFirstHalf, float slopeSecondHalf, float secondHalfOffset) {
        float fi = (float)iIndex;
        float warped = iIndex < halfSize
                ? fi * slopeFirstHalf
                : fi * slopeSecondHalf + secondHalfOffset;
        int indexInteger = (int)warped;
        return indexInteger & max;
    }

    // Single-sample renderer used by non-block paths.
    // Mode behavior:
    // - Decimation enabled: D=2 sample-and-hold + bit-depth quantization.
    // - Decimation disabled, interpolation disabled: plain table lookup (Full).
    // - Decimation disabled, interpolation enabled: lightweight 2-point averaging (HQ).
    inline __attribute__((always_inline)) float getNextSample(struct OscState *oscState)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        int max = waveTable->max;
        float* wave = waveTable->table;
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

            oscState->index += phaseIncrement;

            int indexInteger = oscState->index;
            oscState->index -= indexInteger;
            indexInteger &= max;
            oscState->index += indexInteger;

            if (waveInterpolationEnabled) {
                float currentValue = wave[indexInteger];
                float sample = 0.5f * (currentValue + oscState->waveDecimationHeldSample);
                oscState->waveDecimationHeldSample = currentValue;
                return sample;
            }

            float sample = wave[indexInteger];
            return sample;
        }

        int size = max + 1;
        int halfSize = size >> 1;
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

        oscState->index += phaseIncrement;

        // convert to int;
        int indexInteger = oscState->index;
        // keep decimal part;
        oscState->index -= indexInteger;
        // Put it back inside the table
        indexInteger &= max;
        // Readjust the floating pont inside the table
        oscState->index += indexInteger;

        if (waveInterpolationEnabled) {
            indexInteger = getWarpedIndexFast(indexInteger, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            float currentValue = wave[indexInteger];
            float sample = 0.5f * (currentValue + oscState->waveDecimationHeldSample);
            oscState->waveDecimationHeldSample = currentValue;
            return sample;
        }

        indexInteger = getWarpedIndexFast(indexInteger, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
        float sample = wave[indexInteger];
        return sample;
    }

    inline __attribute__((always_inline)) float getPhase(struct OscState *oscState)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        return oscState->index * waveTable->phaseMul;
    }

    inline __attribute__((always_inline)) float geIndexFromtPhase(float phase)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        return phase * waveTable->max;
    }

    // Block renderer for Full mode. If HQ interpolation is requested, this function
    // delegates to getNextBlockHQ so the HQ path stays centralized.
	inline __attribute__((always_inline)) float* getNextBlock(struct OscState *oscState)  {

        if (unlikely(oscState->waveInterpolationEnabled != 0)) {
            return getNextBlockHQ(oscState);
        }

        int shape = (int) oscillator->shape;
   		int max = waveTables[shape].max;
        int size = max + 1;
        int halfSize = size >> 1;
   		float *wave = waveTables[shape].table;
        float freq = oscState->frequency * waveTables[shape].precomputedValue + waveTables[shape].floatToAdd;
		float freq2 = freq + freq;
        float warp = oscState->effectiveWarp;
        bool warpEnabled = warp != 0.0f;
        bool waveDecimationEnabled = oscState->waveDecimationEnabled != 0;
        float slopeFirstHalf = 1.0f + warp;
        float slopeSecondHalf = 1.0f - warp;
        float secondHalfOffset = ((float)size) * warp;
   		float fIndex = oscState->index;
   		int iIndex;
   		float* oscValuesToFill = oscValues[oscValuesCpt];
    	oscValuesCpt++;
    	oscValuesCpt &= 0x3;

        if (waveDecimationEnabled) {
            int k = 0;
            if (oscState->waveDecimationStepPhase != 0) {
                oscValuesToFill[k++] = oscState->waveDecimationHeldSample;
                oscState->waveDecimationStepPhase = 0;
            }

            if (warpEnabled) {
                for (; k + 1 < BLOCK_SIZE; k += 2) {
                    fIndex += freq2;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &= max;
                    fIndex += iIndex;
                    iIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
                    float sample = quantizeWaveSample(oscState, wave[iIndex]);
                    oscState->waveDecimationHeldSample = sample;
                    oscValuesToFill[k] = sample;
                    oscValuesToFill[k + 1] = sample;
                }
                if (k < BLOCK_SIZE) {
                    fIndex += freq2;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &= max;
                    fIndex += iIndex;
                    iIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
                    float sample = quantizeWaveSample(oscState, wave[iIndex]);
                    oscState->waveDecimationHeldSample = sample;
                    oscValuesToFill[k] = sample;
                    oscState->waveDecimationStepPhase = 1;
                }
            } else {
                for (; k + 1 < BLOCK_SIZE; k += 2) {
                    fIndex += freq2;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &= max;
                    fIndex += iIndex;
                    float sample = quantizeWaveSample(oscState, wave[iIndex]);
                    oscState->waveDecimationHeldSample = sample;
                    oscValuesToFill[k] = sample;
                    oscValuesToFill[k + 1] = sample;
                }
                if (k < BLOCK_SIZE) {
                    fIndex += freq2;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &= max;
                    fIndex += iIndex;
                    float sample = quantizeWaveSample(oscState, wave[iIndex]);
                    oscState->waveDecimationHeldSample = sample;
                    oscValuesToFill[k] = sample;
                    oscState->waveDecimationStepPhase = 1;
                }
            }
			oscState->index = fIndex;
			return oscValuesToFill;
		}

		if (!warpEnabled) {
			for (int k=0; k<32; ) {
                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                float sample = wave[iIndex];
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                sample = wave[iIndex];
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                sample = wave[iIndex];
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                sample = wave[iIndex];
                oscValuesToFill[k++] = sample;

			}

	    	oscState->index = fIndex;
	    	return oscValuesToFill;
		}

		for (int k=0; k<32; ) {
            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            int tableIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            float sample = wave[tableIndex];
            oscValuesToFill[k++] = sample;

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            tableIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            sample = wave[tableIndex];
            oscValuesToFill[k++] = sample;

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            tableIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            sample = wave[tableIndex];
            oscValuesToFill[k++] = sample;

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            tableIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
            sample = wave[tableIndex];
            oscValuesToFill[k++] = sample;

   		}

    	oscState->index = fIndex;
    	return oscValuesToFill;
    };

    // Feedback-capable block renderer used by feedback operator paths.
    // Decimation output quantization is preserved; interpolation is optional and can
    // be force-disabled with PFM3_FAST_FEEDBACK_INTERP for speed-focused builds.
    inline __attribute__((always_inline)) float* getNextBlockWithFeedbackAndEnveloppe(struct OscState *oscState, float feedback, float& env, float envInc, float freqMultiplier, float* lastValue) {
        int shape = (int) oscillator->shape;
        int max = waveTables[shape].max;
        int size = max + 1;
        int halfSize = size >> 1;
        float *wave = waveTables[shape].table;
        float freq = oscState->frequency * waveTables[shape].precomputedValue + waveTables[shape].floatToAdd;
        float fIndex = oscState->index;
        int iIndex;
        float* oscValuesToFill = oscValues[4];

        float warp = oscState->effectiveWarp;
        bool warpEnabled = warp != 0.0f;
    #if PFM3_FAST_FEEDBACK_INTERP
        // Speed-first mode: skip feedback-path interpolation to reduce per-sample math.
        bool interpolationEnabled = false;
    #else
        bool interpolationEnabled = (oscState->waveInterpolationEnabled != 0) && (oscState->waveDecimationEnabled == 0);
    #endif
        float slopeFirstHalf = 1.0f + warp;
        float slopeSecondHalf = 1.0f - warp;
        float secondHalfOffset = ((float)size) * warp;

        lastValue[2] = .95f * lastValue[2] + feedback * .05f;
        float phaseModulationAmplitude = lastValue[2] * ((float) max) * .5f;

        float localLastValue0 = lastValue[0];
        float localLastValue1 = lastValue[1];

        // Optimisation to avoid multiple freqMultiplier in the loop
        float localEnvM = env * freqMultiplier;
        float envIncM   = envInc   * freqMultiplier;
        if (!warpEnabled) {
            if (interpolationEnabled) {
                for (int k = 0; k < 32; k++) {
                    fIndex += freq;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &= max;
                    fIndex += iIndex;

                    float phaseModulationOffset = localLastValue0 * phaseModulationAmplitude;
                    int index = iIndex + (int)phaseModulationOffset;
                    index &= max;

                    float currentValue = wave[index];
                    float newValue = 0.5f * (currentValue + localLastValue1);
                    localLastValue0 = newValue - localLastValue1 + .99525f * localLastValue0;
                    localLastValue1 = newValue;

                    float outputSample = quantizeOscOutputBeforeEnvelope(oscState, localLastValue0);
                    oscValuesToFill[k] = outputSample * localEnvM;
                    localEnvM += envIncM;
                }
            } else {
                for (int k = 0; k < 32; k++) {
                    fIndex += freq;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &= max;
                    fIndex += iIndex;

                    float phaseModulationOffset = localLastValue0 * phaseModulationAmplitude;
                    int index = iIndex + (int)phaseModulationOffset;
                    index &= max;

                    float newValue = wave[index];
                    localLastValue0 = newValue - localLastValue1 + .99525f * localLastValue0;
                    localLastValue1 = newValue;

                    float outputSample = quantizeOscOutputBeforeEnvelope(oscState, localLastValue0);
                    oscValuesToFill[k] = outputSample * localEnvM;
                    localEnvM += envIncM;
                }
            }
        } else {
            if (interpolationEnabled) {
                for (int k = 0; k < 32; k++) {
                    fIndex += freq;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &= max;
                    fIndex += iIndex;

                    int warpedIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);

                    float phaseModulationOffset = localLastValue0 * phaseModulationAmplitude;
                    int index = warpedIndex + (int)phaseModulationOffset;
                    index &= max;

                    float currentValue = wave[index];
                    float newValue = 0.5f * (currentValue + localLastValue1);
                    localLastValue0 = newValue - localLastValue1 + .99525f * localLastValue0;
                    localLastValue1 = newValue;

                    float outputSample = quantizeOscOutputBeforeEnvelope(oscState, localLastValue0);
                    oscValuesToFill[k] = outputSample * localEnvM;
                    localEnvM += envIncM;
                }
            } else {
                for (int k = 0; k < 32; k++) {
                    fIndex += freq;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &= max;
                    fIndex += iIndex;

                    int warpedIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);

                    float phaseModulationOffset = localLastValue0 * phaseModulationAmplitude;
                    int index = warpedIndex + (int)phaseModulationOffset;
                    index &= max;

                    float newValue = wave[index];
                    localLastValue0 = newValue - localLastValue1 + .99525f * localLastValue0;
                    localLastValue1 = newValue;

                    float outputSample = quantizeOscOutputBeforeEnvelope(oscState, localLastValue0);
                    oscValuesToFill[k] = outputSample * localEnvM;
                    localEnvM += envIncM;
                }
            }
        }
        lastValue[0] = localLastValue0;
        lastValue[1] = localLastValue1;

        // update env
        env += envInc * 32;

        oscState->index = fIndex;
        return oscValuesToFill;
    }

    // HQ block renderer. Uses interpolation when decimation is off; when decimation
    // is on it intentionally falls back to decimated hold+quantized output behavior.
	inline __attribute__((always_inline)) float* getNextBlockHQ(struct OscState *oscState)  {
        int shape = (int) oscillator->shape;
   		int max = waveTables[shape].max;
        int size = max + 1;
        int halfSize = size >> 1;
   		float *wave = waveTables[shape].table;
        float freq = oscState->frequency * waveTables[shape].precomputedValue + waveTables[shape].floatToAdd;
        float freq2 = freq + freq;
        float warp = oscState->effectiveWarp;
        bool warpEnabled = warp != 0.0f;
        bool waveDecimationEnabled = oscState->waveDecimationEnabled != 0;
        float slopeFirstHalf = 1.0f + warp;
        float slopeSecondHalf = 1.0f - warp;
        float secondHalfOffset = ((float)size) * warp;
   		float fIndex = oscState->index;
   		int iIndex;
   		float fp;
   		float* oscValuesToFill = oscValues[oscValuesCpt];
    	oscValuesCpt++;
    	oscValuesCpt &= 0x3;

        if (waveDecimationEnabled) {
            int k = 0;
            if (oscState->waveDecimationStepPhase != 0) {
                oscValuesToFill[k++] = oscState->waveDecimationHeldSample;
                oscState->waveDecimationStepPhase = 0;
            }

            if (warpEnabled) {
                for (; k + 1 < BLOCK_SIZE; k += 2) {
                    fIndex += freq2;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &= max;
                    fIndex += iIndex;
                    iIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
                    float sample = quantizeWaveSample(oscState, wave[iIndex]);
                    oscState->waveDecimationHeldSample = sample;
                    oscValuesToFill[k] = sample;
                    oscValuesToFill[k + 1] = sample;
                }
                if (k < BLOCK_SIZE) {
                    fIndex += freq2;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &= max;
                    fIndex += iIndex;
                    iIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
                    float sample = quantizeWaveSample(oscState, wave[iIndex]);
                    oscState->waveDecimationHeldSample = sample;
                    oscValuesToFill[k] = sample;
                    oscState->waveDecimationStepPhase = 1;
                }
            } else {
                for (; k + 1 < BLOCK_SIZE; k += 2) {
                    fIndex += freq2;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &= max;
                    fIndex += iIndex;
                    float sample = quantizeWaveSample(oscState, wave[iIndex]);
                    oscState->waveDecimationHeldSample = sample;
                    oscValuesToFill[k] = sample;
                    oscValuesToFill[k + 1] = sample;
                }
                if (k < BLOCK_SIZE) {
                    fIndex += freq2;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &= max;
                    fIndex += iIndex;
                    float sample = quantizeWaveSample(oscState, wave[iIndex]);
                    oscState->waveDecimationHeldSample = sample;
                    oscValuesToFill[k] = sample;
                    oscState->waveDecimationStepPhase = 1;
                }
            }
            oscState->index = fIndex;
            return oscValuesToFill;
        }

		if (!warpEnabled) {
			float previousValue = oscState->waveDecimationHeldSample;
			for (int k=0; k<32; ) {
                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
        float currentValue = wave[iIndex];
        float sample = 0.5f * (currentValue + previousValue);
        previousValue = currentValue;
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                currentValue = wave[iIndex];
                sample = 0.5f * (currentValue + previousValue);
                previousValue = currentValue;
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                currentValue = wave[iIndex];
                sample = 0.5f * (currentValue + previousValue);
                previousValue = currentValue;
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                currentValue = wave[iIndex];
                sample = 0.5f * (currentValue + previousValue);
                previousValue = currentValue;
                oscValuesToFill[k++] = sample;
			}

			oscState->waveDecimationHeldSample = previousValue;

	    	oscState->index = fIndex;
	    	return oscValuesToFill;
		}

			float previousValue = oscState->waveDecimationHeldSample;
			for (int k=0; k<32; ) {
                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                int tableIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
                float currentValue = wave[tableIndex];
                float sample = 0.5f * (currentValue + previousValue);
                previousValue = currentValue;
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                tableIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
                currentValue = wave[tableIndex];
                sample = 0.5f * (currentValue + previousValue);
                previousValue = currentValue;
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                tableIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
                currentValue = wave[tableIndex];
                sample = 0.5f * (currentValue + previousValue);
                previousValue = currentValue;
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                tableIndex = getWarpedIndexFast(iIndex, max, halfSize, slopeFirstHalf, slopeSecondHalf, secondHalfOffset);
                currentValue = wave[tableIndex];
                sample = 0.5f * (currentValue + previousValue);
                previousValue = currentValue;
                oscValuesToFill[k++] = sample;
			}

			oscState->waveDecimationHeldSample = previousValue;

    	oscState->index = fIndex;
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
