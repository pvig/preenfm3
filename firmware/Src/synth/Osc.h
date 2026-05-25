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
};


extern struct WaveTable waveTables[];
extern float exp2_harm[];


class Osc : public SynthStateAware
{
public:
    Osc() {};
    virtual ~Osc() {};

    void init(SynthState* synthState, struct OscillatorParams *oscParams, DestinationEnum df);

    void newNote(struct OscState* oscState, float newNoteFrequency, float phase);
    float getNoteRealFrequencyEstimation(struct OscState* oscState, float newNoteFrequency);
    void glideToNote(struct OscState* oscState, float newNoteFrequency);
    void glideStep(struct OscState* oscState, float phase);

    inline void calculateFrequencyWithMatrix(struct OscState *oscState, Matrix* matrix, float expHarm) {
        oscState->mainFrequencyPlusMatrix = oscState->mainFrequency;
        oscState->mainFrequencyPlusMatrix *= expHarm;
        oscState->mainFrequencyPlusMatrix +=  (oscState->mainFrequency  * (matrix->getDestination(destFreq) + matrix->getDestination(ALL_OSC_FREQ)) * .1f);
    }

    inline float quantizeWaveSample(float value, uint8_t bits) {
        float scale = (float)(1u << bits);
        float scaled = value * scale;
        // Keep rounding in float domain to avoid integer overflow for large FM sums.
        float q = (scaled >= 0.0f) ? floorf(scaled + 0.5f) : ceilf(scaled - 0.5f);
        return q * (1.0f / scale);
    }

    inline float quantizeFrequencyWithFmSum(struct OscState *oscState) {
        if (!oscState->waveDecimationEnabled) {
            return oscState->frequency;
        }

        float fmSum = oscState->frequency - oscState->mainFrequencyPlusMatrix;
        float quantizedFmSum = quantizeWaveSample(fmSum, oscState->waveDecimationBits);
        return oscState->mainFrequencyPlusMatrix + quantizedFmSum;
    }

    inline float quantizePhaseIncrement(struct OscState *oscState, float phaseIncrement) {
        if (!oscState->waveDecimationEnabled) {
            return phaseIncrement;
        }
        return quantizeWaveSample(phaseIncrement, oscState->waveDecimationBits);
    }

    inline float quantizePhaseAccumulator(struct OscState *oscState, float index, int max) {
        if (!oscState->waveDecimationEnabled) {
            return index;
        }

        int indexInteger = (int)floorf(index);
        float frac = index - (float)indexInteger;
        frac = quantizeWaveSample(frac, oscState->waveDecimationBits);
        if (unlikely(frac >= 1.0f)) {
            // Keep integer index stable for current-sample lookup; clamp frac to max representable step.
            frac = 1.0f - (1.0f / (float)(1u << oscState->waveDecimationBits));
        }
        indexInteger &= max;
        return (float)indexInteger + frac;
    }

    inline float quantizeWaveInputSample(struct OscState *oscState, float sample) {
        if (!oscState->waveDecimationEnabled) {
            return sample;
        }
        return quantizeWaveSample(sample, oscState->waveDecimationBits);
    }

    inline float quantizePhaseModulationOffset(struct OscState *oscState, float phaseModulationOffset) {
        if (!oscState->waveDecimationEnabled) {
            return phaseModulationOffset;
        }
        return quantizeWaveSample(phaseModulationOffset, oscState->waveDecimationBits);
    }

    inline float quantizeOscOutputBeforeEnvelope(struct OscState *oscState, float outputSample) {
        if (!oscState->waveDecimationEnabled) {
            return outputSample;
        }
        return quantizeWaveSample(outputSample, oscState->waveDecimationBits);
    }

    inline float getNextSample(struct OscState *oscState)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        float quantizedFrequency = quantizeFrequencyWithFmSum(oscState);
        float phaseIncrement = quantizedFrequency * waveTable->precomputedValue + waveTable->floatToAdd;
        phaseIncrement = quantizePhaseIncrement(oscState, phaseIncrement);

        oscState->index += phaseIncrement;

        // convert to int;
        int indexInteger = oscState->index;
        // keep decimal part;
        oscState->index -= indexInteger;
        // Put it back inside the table
        indexInteger &= waveTable->max;
        // Readjust the floating pont inside the table
        oscState->index += indexInteger;
        oscState->index = quantizePhaseAccumulator(oscState, oscState->index, waveTable->max);

        float sample = waveTable->table[indexInteger];
        sample = quantizeWaveInputSample(oscState, sample);
        return quantizeOscOutputBeforeEnvelope(oscState, sample);
    }

    inline float getPhase(struct OscState *oscState)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        return oscState->index * waveTable->phaseMul;
    }

    inline float geIndexFromtPhase(float phase)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        return phase * waveTable->max;
    }


   	inline float* getNextBlock(struct OscState *oscState)  {
        int shape = (int) oscillator->shape;
   		int max = waveTables[shape].max;
   		float *wave = waveTables[shape].table;
        float quantizedFrequency = quantizeFrequencyWithFmSum(oscState);
        float freq = quantizedFrequency * waveTables[shape].precomputedValue + waveTables[shape].floatToAdd;
        freq = quantizePhaseIncrement(oscState, freq);
   		float fIndex = oscState->index;
   		int iIndex;
   		float* oscValuesToFill = oscValues[oscValuesCpt];
    	oscValuesCpt++;
    	oscValuesCpt &= 0x3;

   		for (int k=0; k<32; ) {
            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            fIndex = quantizePhaseAccumulator(oscState, fIndex, max);
            float sample = quantizeWaveInputSample(oscState, wave[iIndex]);
            oscValuesToFill[k++] = quantizeOscOutputBeforeEnvelope(oscState, sample);

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            fIndex = quantizePhaseAccumulator(oscState, fIndex, max);
            sample = quantizeWaveInputSample(oscState, wave[iIndex]);
            oscValuesToFill[k++] = quantizeOscOutputBeforeEnvelope(oscState, sample);

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            fIndex = quantizePhaseAccumulator(oscState, fIndex, max);
            sample = quantizeWaveInputSample(oscState, wave[iIndex]);
            oscValuesToFill[k++] = quantizeOscOutputBeforeEnvelope(oscState, sample);

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            fIndex = quantizePhaseAccumulator(oscState, fIndex, max);
            sample = quantizeWaveInputSample(oscState, wave[iIndex]);
            oscValuesToFill[k++] = quantizeOscOutputBeforeEnvelope(oscState, sample);

   		}
    	oscState->index = fIndex;
    	return oscValuesToFill;
    };


    inline float* getNextBlockWithFeedbackAndEnveloppe(struct OscState *oscState, float feedback, float& env, float envInc, float freqMultiplier, float* lastValue) {
        int shape = (int) oscillator->shape;
        int max = waveTables[shape].max;
        float *wave = waveTables[shape].table;
        float quantizedFrequency = quantizeFrequencyWithFmSum(oscState);
        float freq = quantizedFrequency * waveTables[shape].precomputedValue + waveTables[shape].floatToAdd;
        freq = quantizePhaseIncrement(oscState, freq);
        float fIndex = oscState->index;
        int iIndex;
        float* oscValuesToFill = oscValues[4];

        lastValue[2] = .95f * lastValue[2] + feedback * .05f;
        float phaseModulationAmplitude = lastValue[2] * ((float) max) * .5f;

        float localLastValue0 = lastValue[0];
        float localLastValue1 = lastValue[1];

        // Optimisation to avoid multiple freqMultiplier in the loop
        float localEnvM = env * freqMultiplier;
        float envIncM   = envInc   * freqMultiplier;
        for (int k = 0; k < 32; k++) {
            fIndex += freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &= max;
            fIndex += iIndex;
            fIndex = quantizePhaseAccumulator(oscState, fIndex, max);

            float phaseModulationOffset = quantizePhaseModulationOffset(oscState, localLastValue0 * phaseModulationAmplitude);
            int index = iIndex + (int)phaseModulationOffset;
            index &= max;

            // Get rid of DC offset
            float newValue = quantizeWaveInputSample(oscState, wave[index]);
            localLastValue0 = newValue - localLastValue1 + .99525f * localLastValue0;
            localLastValue1 = newValue;

            float outputSample = quantizeOscOutputBeforeEnvelope(oscState, localLastValue0);
            oscValuesToFill[k] = outputSample * localEnvM;
            localEnvM += envIncM;
        }
        lastValue[0] = localLastValue0;
        lastValue[1] = localLastValue1;

        // update env
        env += envInc * 32;

        oscState->index = fIndex;
        return oscValuesToFill;
    }


   	float* getNextBlockHQ(struct OscState *oscState)  {
        int shape = (int) oscillator->shape;
   		int max = waveTables[shape].max;
   		float *wave = waveTables[shape].table;
        float quantizedFrequency = quantizeFrequencyWithFmSum(oscState);
        float freq = quantizedFrequency * waveTables[shape].precomputedValue + waveTables[shape].floatToAdd;
        freq = quantizePhaseIncrement(oscState, freq);
   		float fIndex = oscState->index;
   		int iIndex;
   		float fp;
   		float* oscValuesToFill = oscValues[oscValuesCpt];
    	oscValuesCpt++;
    	oscValuesCpt &= 0x3;
   		for (int k=0; k<32; ) {
            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            fIndex = quantizePhaseAccumulator(oscState, fIndex, max);
            fp = fIndex - floorf(fIndex);
            float sample;
            if (oscState->waveDecimationEnabled) {
                sample = quantizeWaveInputSample(oscState, wave[iIndex]);
            } else {
                int iIndexNext = (iIndex + 1) & max;
                sample = wave[iIndex] * (1-fp) + wave[iIndexNext] * fp;
            }
            oscValuesToFill[k++] = quantizeOscOutputBeforeEnvelope(oscState, sample);

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            fIndex = quantizePhaseAccumulator(oscState, fIndex, max);
            fp = fIndex - floorf(fIndex);
            if (oscState->waveDecimationEnabled) {
                sample = quantizeWaveInputSample(oscState, wave[iIndex]);
            } else {
                int iIndexNext = (iIndex + 1) & max;
                sample = wave[iIndex] * (1-fp) + wave[iIndexNext] * fp;
            }
            oscValuesToFill[k++] = quantizeOscOutputBeforeEnvelope(oscState, sample);

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            fIndex = quantizePhaseAccumulator(oscState, fIndex, max);
            fp = fIndex - floorf(fIndex);
            if (oscState->waveDecimationEnabled) {
                sample = quantizeWaveInputSample(oscState, wave[iIndex]);
            } else {
                int iIndexNext = (iIndex + 1) & max;
                sample = wave[iIndex] * (1-fp) + wave[iIndexNext] * fp;
            }
            oscValuesToFill[k++] = quantizeOscOutputBeforeEnvelope(oscState, sample);

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            fIndex = quantizePhaseAccumulator(oscState, fIndex, max);
            fp = fIndex - floorf(fIndex);
            if (oscState->waveDecimationEnabled) {
                sample = quantizeWaveInputSample(oscState, wave[iIndex]);
            } else {
                int iIndexNext = (iIndex + 1) & max;
                sample = wave[iIndex] * (1-fp) + wave[iIndexNext] * fp;
            }
            oscValuesToFill[k++] = quantizeOscOutputBeforeEnvelope(oscState, sample);
   		}
    	oscState->index = fIndex;
    	return oscValuesToFill;
    };


private:
    DestinationEnum destFreq;
    static float* oscValues[5];
    static int oscValuesCpt;
    OscillatorParams* oscillator;
};
