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
    float index;
    float frequency;
    float mainFrequencyPlusMatrix;
    float mainFrequency;
    float fromFrequency;
    float nextFrequency;
    uint8_t waveDecimationBits;
    uint8_t waveDecimationEnabled;
    uint8_t waveDecimationStepPhase;
    float waveDecimationScale;
    float waveDecimationInvScale;
    float waveDecimationHeldSample;
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

    inline __attribute__((always_inline)) int getWarpedIndex(int iIndex, int max, int halfSize, int size, float slopeFirstHalf, float slopeSecondHalf) {
        float warped = iIndex < halfSize
                ? ((float)iIndex) * slopeFirstHalf
                : (float)size - ((float)(size - iIndex)) * slopeSecondHalf;
        int indexInteger = (int)warped;
        indexInteger &= max;
        return indexInteger;
    }

    inline __attribute__((always_inline)) float getNextSample(struct OscState *oscState)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        float phaseIncrement = oscState->frequency * waveTable->precomputedValue + waveTable->floatToAdd;
        float warp = oscState->effectiveWarp;
        bool waveDecimationEnabled = oscState->waveDecimationEnabled != 0;

        if (likely(warp == 0.0f)) {
            if (waveDecimationEnabled) {
                if (oscState->waveDecimationStepPhase != 0) {
                    oscState->waveDecimationStepPhase = 0;
                    return oscState->waveDecimationHeldSample;
                }

                oscState->waveDecimationStepPhase = 1;
                float phaseIncrement2 = phaseIncrement + phaseIncrement;
                oscState->index += phaseIncrement2;

                int indexInteger = oscState->index;
                oscState->index -= indexInteger;
                indexInteger &= waveTable->max;
                oscState->index += indexInteger;

                float sample = waveTable->table[indexInteger];
                sample = quantizeOscOutputBeforeEnvelope(oscState, sample);
                oscState->waveDecimationHeldSample = sample;
                return sample;
            }

            oscState->index += phaseIncrement;

            int indexInteger = oscState->index;
            oscState->index -= indexInteger;
            indexInteger &= waveTable->max;
            oscState->index += indexInteger;
            float sample = waveTable->table[indexInteger];
            return sample;
        }

        int size = waveTable->max + 1;
        int halfSize = size >> 1;
        float slopeFirstHalf = 1.0f + warp;
        float slopeSecondHalf = 1.0f - warp;

        if (waveDecimationEnabled) {
            if (oscState->waveDecimationStepPhase != 0) {
                oscState->waveDecimationStepPhase = 0;
                return oscState->waveDecimationHeldSample;
            }

            oscState->waveDecimationStepPhase = 1;
            float phaseIncrement2 = phaseIncrement + phaseIncrement;
            oscState->index += phaseIncrement2;

            int indexInteger = oscState->index;
            oscState->index -= indexInteger;
            indexInteger &= waveTable->max;
            oscState->index += indexInteger;

            indexInteger = getWarpedIndex(indexInteger, waveTable->max, halfSize, size, slopeFirstHalf, slopeSecondHalf);

            float sample = waveTable->table[indexInteger];
            sample = quantizeOscOutputBeforeEnvelope(oscState, sample);
            oscState->waveDecimationHeldSample = sample;
            return sample;
        }

        oscState->index += phaseIncrement;

        // convert to int;
        int indexInteger = oscState->index;
        // keep decimal part;
        oscState->index -= indexInteger;
        // Put it back inside the table
        indexInteger &= waveTable->max;
        // Readjust the floating pont inside the table
        oscState->index += indexInteger;
        indexInteger = getWarpedIndex(indexInteger, waveTable->max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
        float sample = waveTable->table[indexInteger];
        return sample;
    }

    inline __attribute__((always_inline)) float getPhase(struct OscState *oscState)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        return oscState->index * waveTable->phaseMul;
    }

    inline __attribute__((always_inline)) float getNextDecimatedSampleNoWarp(struct OscState *oscState, float& fIndex, float freq2, int max, float *wave) {
        fIndex += freq2;
        int iIndex = fIndex;
        fIndex -= iIndex;
        iIndex &= max;
        fIndex += iIndex;

        float sample = wave[iIndex];
        sample = quantizeOscOutputBeforeEnvelope(oscState, sample);
        oscState->waveDecimationHeldSample = sample;
        return sample;
    }

    inline __attribute__((always_inline)) float getNextDecimatedSampleWarp(struct OscState *oscState, float& fIndex, float freq2, int max, float *wave,
            int halfSize, int size, float slopeFirstHalf, float slopeSecondHalf) {
        fIndex += freq2;
        int iIndex = fIndex;
        fIndex -= iIndex;
        iIndex &= max;
        fIndex += iIndex;

        iIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);

        float sample = wave[iIndex];
        sample = quantizeOscOutputBeforeEnvelope(oscState, sample);
        oscState->waveDecimationHeldSample = sample;
        return sample;
    }

    inline __attribute__((always_inline)) float* fillDecimatedBlockNoWarp(struct OscState *oscState, float *oscValuesToFill, float& fIndex, float freq2, int max, float *wave) {
        int k = 0;

        if (oscState->waveDecimationStepPhase != 0) {
            oscValuesToFill[k++] = oscState->waveDecimationHeldSample;
            oscState->waveDecimationStepPhase = 0;
        }

        for (; k + 1 < BLOCK_SIZE; k += 2) {
            float sample = getNextDecimatedSampleNoWarp(oscState, fIndex, freq2, max, wave);
            oscValuesToFill[k] = sample;
            oscValuesToFill[k + 1] = sample;
        }

        if (k < BLOCK_SIZE) {
            float sample = getNextDecimatedSampleNoWarp(oscState, fIndex, freq2, max, wave);
            oscValuesToFill[k] = sample;
            oscState->waveDecimationStepPhase = 1;
        }

        return oscValuesToFill;
    }

    inline __attribute__((always_inline)) float* fillDecimatedBlockWarp(struct OscState *oscState, float *oscValuesToFill, float& fIndex, float freq2, int max, float *wave,
            int halfSize, int size, float slopeFirstHalf, float slopeSecondHalf) {
        int k = 0;

        if (oscState->waveDecimationStepPhase != 0) {
            oscValuesToFill[k++] = oscState->waveDecimationHeldSample;
            oscState->waveDecimationStepPhase = 0;
        }

        for (; k + 1 < BLOCK_SIZE; k += 2) {
            float sample = getNextDecimatedSampleWarp(oscState, fIndex, freq2, max, wave, halfSize, size, slopeFirstHalf, slopeSecondHalf);
            oscValuesToFill[k] = sample;
            oscValuesToFill[k + 1] = sample;
        }

        if (k < BLOCK_SIZE) {
            float sample = getNextDecimatedSampleWarp(oscState, fIndex, freq2, max, wave, halfSize, size, slopeFirstHalf, slopeSecondHalf);
            oscValuesToFill[k] = sample;
            oscState->waveDecimationStepPhase = 1;
        }

        return oscValuesToFill;
    }

    inline __attribute__((always_inline)) float geIndexFromtPhase(float phase)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        return phase * waveTable->max;
    }


   	inline __attribute__((always_inline)) float* getNextBlock(struct OscState *oscState)  {
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
                    iIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
                    float sample = quantizeOscOutputBeforeEnvelope(oscState, wave[iIndex]);
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
                    iIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
                    float sample = quantizeOscOutputBeforeEnvelope(oscState, wave[iIndex]);
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
                    float sample = quantizeOscOutputBeforeEnvelope(oscState, wave[iIndex]);
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
                    float sample = quantizeOscOutputBeforeEnvelope(oscState, wave[iIndex]);
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
            int tableIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
            float sample = wave[tableIndex];
            oscValuesToFill[k++] = sample;

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            tableIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
            sample = wave[tableIndex];
            oscValuesToFill[k++] = sample;

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            tableIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
            sample = wave[tableIndex];
            oscValuesToFill[k++] = sample;

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            tableIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
            sample = wave[tableIndex];
            oscValuesToFill[k++] = sample;

   		}

    	oscState->index = fIndex;
    	return oscValuesToFill;
    };


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
        float slopeFirstHalf = 1.0f + warp;
        float slopeSecondHalf = 1.0f - warp;

        lastValue[2] = .95f * lastValue[2] + feedback * .05f;
        float phaseModulationAmplitude = lastValue[2] * ((float) max) * .5f;

        float localLastValue0 = lastValue[0];
        float localLastValue1 = lastValue[1];

        // Optimisation to avoid multiple freqMultiplier in the loop
        float localEnvM = env * freqMultiplier;
        float envIncM   = envInc   * freqMultiplier;
        if (!warpEnabled) {
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
        } else {
            for (int k = 0; k < 32; k++) {
                fIndex += freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &= max;
                fIndex += iIndex;

                iIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);

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
        lastValue[0] = localLastValue0;
        lastValue[1] = localLastValue1;

        // update env
        env += envInc * 32;

        oscState->index = fIndex;
        return oscValuesToFill;
    }


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
                    iIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
                    float sample = quantizeOscOutputBeforeEnvelope(oscState, wave[iIndex]);
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
                    iIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
                    float sample = quantizeOscOutputBeforeEnvelope(oscState, wave[iIndex]);
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
                    float sample = quantizeOscOutputBeforeEnvelope(oscState, wave[iIndex]);
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
                    float sample = quantizeOscOutputBeforeEnvelope(oscState, wave[iIndex]);
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
                fp = fIndex - (float)iIndex;
                float sample;
                int iIndexNext = (iIndex + 1) & max;
                sample = wave[iIndex] * (1-fp) + wave[iIndexNext] * fp;
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                fp = fIndex - (float)iIndex;
                iIndexNext = (iIndex + 1) & max;
                sample = wave[iIndex] * (1-fp) + wave[iIndexNext] * fp;
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                fp = fIndex - (float)iIndex;
                iIndexNext = (iIndex + 1) & max;
                sample = wave[iIndex] * (1-fp) + wave[iIndexNext] * fp;
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                fp = fIndex - (float)iIndex;
                iIndexNext = (iIndex + 1) & max;
                sample = wave[iIndex] * (1-fp) + wave[iIndexNext] * fp;
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
                fp = fIndex - (float)iIndex;
                float sample;
                int tableIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
                int iIndexNext = (iIndex + 1) & max;
                int tableIndexNext = getWarpedIndex(iIndexNext, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
                sample = wave[tableIndex] * (1-fp) + wave[tableIndexNext] * fp;
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                fp = fIndex - (float)iIndex;
                tableIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
                iIndexNext = (iIndex + 1) & max;
                tableIndexNext = getWarpedIndex(iIndexNext, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
                sample = wave[tableIndex] * (1-fp) + wave[tableIndexNext] * fp;
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                fp = fIndex - (float)iIndex;
                tableIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
                iIndexNext = (iIndex + 1) & max;
                tableIndexNext = getWarpedIndex(iIndexNext, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
                sample = wave[tableIndex] * (1-fp) + wave[tableIndexNext] * fp;
                oscValuesToFill[k++] = sample;

                fIndex +=  freq;
                iIndex = fIndex;
                fIndex -= iIndex;
                iIndex &=  max;
                fIndex += iIndex;
                fp = fIndex - (float)iIndex;
                tableIndex = getWarpedIndex(iIndex, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
                iIndexNext = (iIndex + 1) & max;
                tableIndexNext = getWarpedIndex(iIndexNext, max, halfSize, size, slopeFirstHalf, slopeSecondHalf);
                sample = wave[tableIndex] * (1-fp) + wave[tableIndexNext] * fp;
                oscValuesToFill[k++] = sample;
			}

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
