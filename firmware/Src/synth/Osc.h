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

    inline float quantizeWaveSample(struct OscState *oscState, float value) {
        // Quantize with fixed-point style fractional precision to preserve FM behavior.
        float scaled = value * oscState->waveDecimationScale;
        if (scaled > 2147483520.0f) {
            scaled = 2147483520.0f;
        } else if (scaled < -2147483520.0f) {
            scaled = -2147483520.0f;
        }
        int q = (int)(scaled + (scaled >= 0.0f ? 0.5f : -0.5f));
        return (float)q * oscState->waveDecimationInvScale;
    }

    inline float quantizeOscOutputBeforeEnvelope(struct OscState *oscState, float outputSample) {
        if (!oscState->waveDecimationEnabled) {
            return outputSample;
        }
        return quantizeWaveSample(oscState, outputSample);
    }

    inline float getNextSample(struct OscState *oscState)  {
        struct WaveTable* waveTable = &waveTables[(int) oscillator->shape];
        float phaseIncrement = oscState->frequency * waveTable->precomputedValue + waveTable->floatToAdd;

        if (oscState->waveDecimationEnabled) {
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

        // convert to int;
        int indexInteger = oscState->index;
        // keep decimal part;
        oscState->index -= indexInteger;
        // Put it back inside the table
        indexInteger &= waveTable->max;
        // Readjust the floating pont inside the table
        oscState->index += indexInteger;
        float sample = waveTable->table[indexInteger];
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
        float freq = oscState->frequency * waveTables[shape].precomputedValue + waveTables[shape].floatToAdd;
		float freq2 = freq + freq;
   		float fIndex = oscState->index;
   		int iIndex;
   		float* oscValuesToFill = oscValues[oscValuesCpt];
    	oscValuesCpt++;
    	oscValuesCpt &= 0x3;

		if (oscState->waveDecimationEnabled) {
			for (int k=0; k<32; k++) {
                if (oscState->waveDecimationStepPhase == 0) {
                    fIndex += freq2;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &=  max;
                    fIndex += iIndex;
                    float sample = wave[iIndex];
                    sample = quantizeOscOutputBeforeEnvelope(oscState, sample);
                    oscState->waveDecimationHeldSample = sample;
                    oscState->waveDecimationStepPhase = 1;
                    oscValuesToFill[k] = sample;
                } else {
                    oscState->waveDecimationStepPhase = 0;
                    oscValuesToFill[k] = oscState->waveDecimationHeldSample;
                }
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
            float sample = wave[iIndex];
            oscValuesToFill[k++] = quantizeOscOutputBeforeEnvelope(oscState, sample);

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            sample = wave[iIndex];
            oscValuesToFill[k++] = quantizeOscOutputBeforeEnvelope(oscState, sample);

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            sample = wave[iIndex];
            oscValuesToFill[k++] = quantizeOscOutputBeforeEnvelope(oscState, sample);

            fIndex +=  freq;
            iIndex = fIndex;
            fIndex -= iIndex;
            iIndex &=  max;
            fIndex += iIndex;
            sample = wave[iIndex];
            oscValuesToFill[k++] = quantizeOscOutputBeforeEnvelope(oscState, sample);

   		}

    	oscState->index = fIndex;
    	return oscValuesToFill;
    };


    inline float* getNextBlockWithFeedbackAndEnveloppe(struct OscState *oscState, float feedback, float& env, float envInc, float freqMultiplier, float* lastValue) {
        int shape = (int) oscillator->shape;
        int max = waveTables[shape].max;
        float *wave = waveTables[shape].table;
        float freq = oscState->frequency * waveTables[shape].precomputedValue + waveTables[shape].floatToAdd;
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
        float freq = oscState->frequency * waveTables[shape].precomputedValue + waveTables[shape].floatToAdd;
        float freq2 = freq + freq;
   		float fIndex = oscState->index;
   		int iIndex;
   		float fp;
   		float* oscValuesToFill = oscValues[oscValuesCpt];
    	oscValuesCpt++;
    	oscValuesCpt &= 0x3;

        if (oscState->waveDecimationEnabled) {
            for (int k=0; k<32; k++) {
                if (oscState->waveDecimationStepPhase == 0) {
                    fIndex +=  freq2;
                    iIndex = fIndex;
                    fIndex -= iIndex;
                    iIndex &=  max;
                    fIndex += iIndex;
                    float sample = wave[iIndex];
                    sample = quantizeOscOutputBeforeEnvelope(oscState, sample);
                    oscState->waveDecimationHeldSample = sample;
                    oscState->waveDecimationStepPhase = 1;
                    oscValuesToFill[k] = sample;
                } else {
                    oscState->waveDecimationStepPhase = 0;
                    oscValuesToFill[k] = oscState->waveDecimationHeldSample;
                }
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
                if (oscState->waveDecimationEnabled) {
                    sample = wave[iIndex];
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
                fp = fIndex - (float)iIndex;
                if (oscState->waveDecimationEnabled) {
                    sample = wave[iIndex];
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
                fp = fIndex - (float)iIndex;
                if (oscState->waveDecimationEnabled) {
                    sample = wave[iIndex];
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
                fp = fIndex - (float)iIndex;
                if (oscState->waveDecimationEnabled) {
                    sample = wave[iIndex];
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
