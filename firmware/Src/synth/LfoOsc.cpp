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

#include "LfoOsc.h"
#include "LfoCurveTables.h"
#include <math.h>



extern float noise[32];

namespace {
float tableLookupLinear(const float* table, float x) {
    if (x <= 0.0f) {
        return table[0];
    }
    if (x >= 1.0f) {
        return table[LFO_CURVE_TABLE_SIZE];
    }

    float tablePos = x * (float)LFO_CURVE_TABLE_SIZE;
    int index = (int)tablePos;
    float frac = tablePos - (float)index;
    return table[index] + (table[index + 1] - table[index]) * frac;
}

float waveLookupLinear(const float* table, int max, float phase01) {
    phase01 -= (int)phase01;
    if (phase01 < 0.0f) {
        phase01 += 1.0f;
    }

    float tablePos = phase01 * (float)max;
    int index0 = (int)tablePos;
    float frac = tablePos - (float)index0;
    int index1 = (index0 + 1) & max;

    float sample0 = table[index0];
    return sample0 + (table[index1] - sample0) * frac;
}

float sinLookup(float phase01) {
    return waveLookupLinear(sinTable, waveTables[0].max, phase01);
}
}

void LfoOsc::executeTriangle(float& lfoValue, float phase, bool phaseWrapped) {
    (void)phaseWrapped;
    if (phase < .5f) {
        lfoValue = phase * 4.0f - 1.0f;
    } else {
        lfoValue = 1.0f - (phase - .5f) * 4.0f;
    }
}

void LfoOsc::executeSawUp(float& lfoValue, float phase, bool phaseWrapped) {
    (void)phaseWrapped;
    lfoValue = -1.0f + phase * 2.0f;
}

void LfoOsc::executeSawDown(float& lfoValue, float phase, bool phaseWrapped) {
    (void)phaseWrapped;
    lfoValue = 1.0f - phase * 2.0f;
}

void LfoOsc::executeDecayExp(float& lfoValue, float phase, bool phaseWrapped) {
    (void)phaseWrapped;
    float norm = tableLookupLinear(lfoDecayExpK6Table, phase);
    lfoValue = norm * 2.0f - 1.0f;
}

void LfoOsc::executeDecayLog(float& lfoValue, float phase, bool phaseWrapped) {
    (void)phaseWrapped;
    float norm = tableLookupLinear(lfoDecayLogA31Table, phase);
    lfoValue = norm * 2.0f - 1.0f;
}

void LfoOsc::executeRiseExp(float& lfoValue, float phase, bool phaseWrapped) {
    (void)phaseWrapped;
    float norm = 1.0f - tableLookupLinear(lfoDecayExpK6Table, phase);
    lfoValue = norm * 2.0f - 1.0f;
}

void LfoOsc::executeRiseLog(float& lfoValue, float phase, bool phaseWrapped) {
    (void)phaseWrapped;
    float norm = 1.0f - tableLookupLinear(lfoDecayLogA31Table, phase);
    lfoValue = norm * 2.0f - 1.0f;
}

void LfoOsc::executeAttackDecay(float& lfoValue, float phase, bool phaseWrapped) {
    (void)phaseWrapped;
    if (phase < 0.5f) {
        float t = phase * adAttackPhaseScale;
        float s = t * t * (3.0f - 2.0f * t);
        lfoValue = -1.0f + 2.0f * s;
    } else {
        float t = (phase - adDecayPhaseOffset) * adDecayPhaseScale;
        float s = t * t * (3.0f - 2.0f * t);
        lfoValue = 1.0f - 2.0f * s;
    }
}

void LfoOsc::executeAttackHoldDecay(float& lfoValue, float phase, bool phaseWrapped) {
    (void)phaseWrapped;
    if (phase < 0.25f) {
        float t = phase * ahdAttackPhaseScale;
        float s = t * t * (3.0f - 2.0f * t);
        lfoValue = -1.0f + 2.0f * s;
    } else if (phase < 0.5f) {
        lfoValue = 1.0f;
    } else {
        float t = (phase - ahdDecayPhaseOffset) * ahdDecayPhaseScale;
        float s = t * t * (3.0f - 2.0f * t);
        lfoValue = 1.0f - 2.0f * s;
    }
}

void LfoOsc::executeDecayS(float& lfoValue, float phase, bool phaseWrapped) {
    (void)phaseWrapped;
    float s = phase * phase * (3.0f - 2.0f * phase);
    lfoValue = 1.0f - 2.0f * s;
}

void LfoOsc::executeBuchla(float& lfoValue, float phase, bool phaseWrapped) {
    (void)phaseWrapped;
    if (phase < shapeAttack) {
        float t = phase * shapeInvAttack;
        float s = t * t * (3.0f - 2.0f * t);
        lfoValue = -1.0f + 2.0f * s;
    } else if (phase < (shapeAttack + shapeHold)) {
        lfoValue = 1.0f;
    } else {
        float d = (phase - shapeAttack - shapeHold) * shapeInvDecay;
        float norm = tableLookupLinear(lfoDecayExpK7Table, d);
        lfoValue = norm * 2.0f - 1.0f;
        float bump = 0.10f * tableLookupLinear(lfoDampingExp18Table, d) * sinLookup(d * 3.0f);
        lfoValue += bump;
        if (lfoValue > 1.0f) {
            lfoValue = 1.0f;
        } else if (lfoValue < -1.0f) {
            lfoValue = -1.0f;
        }
    }
}

void LfoOsc::executeWaveTable(float& lfoValue, float phase, bool phaseWrapped) {
    (void)phaseWrapped;
    lfoValue = waveLookupLinear(shapeWaveTable, shapeWaveMax, phase);
}

void LfoOsc::executeSquare(float& lfoValue, float phase, bool phaseWrapped) {
    (void)phaseWrapped;
    lfoValue = phase < .5f ? -1.0f : 1.0f;
}

void LfoOsc::executeRandomFamily(float& lfoValue, float phase, bool phaseWrapped) {
    switch (randomRuntimeMode) {
    case 1: // random
        if (phaseWrapped) {
            currentRandomValue = noise[0];
        }
        lfoValue = currentRandomValue;
        break;
    case 2: // brownian
        if (phaseWrapped) {
            noiseLp = noise[0] * 0.4f + noiseLp * 0.6f;
            currentRandomValue = noiseLp;
        }
        lfoValue = currentRandomValue;
        break;
    case 3: // wandering
        if (phaseWrapped) {
            currentRandomValue = nextRandomValue;
            nextRandomValue = noise[0];
        }
        lfoValue = phase * (nextRandomValue - currentRandomValue) + currentRandomValue;
        break;
    case 4: // flow
        if (phaseWrapped) {
            noiseLp = noise[0] * 0.4f + noiseLp * 0.6f;
            currentRandomValue = nextRandomValue;
            nextRandomValue = noiseLp;
        }
        lfoValue = phase * (nextRandomValue - currentRandomValue) + currentRandomValue;
        break;
    default:
        lfoValue = -1.0f;
        break;
    }
}

void LfoOsc::init(struct LfoParams *lfoParams, float* lfoSyncMode, float* phase, Matrix *matrix, SourceEnum source, DestinationEnum dest) {
    Lfo::init(matrix, source, dest);
    this->type = LFO_TRIANGLE;
    this->ramp = 0;
    this->initPhase = phase;
    this->syncMode = lfoSyncMode;
    this->rampInv = 10000000 ;
    this->currentRamp = 0;
    this->lfo = lfoParams;
    valueChanged(ENCODER_LFO_SHAPE);
    valueChanged(3);
    this->destination = dest;
    this->currentRandomValue = 0.0f;
    this->startupDelaySeconds = 0.0f;

    ticks = 1536;
    midiClock(0, true);
}


void LfoOsc::midiClock(int songPosition, bool computeStep) {

    ticks &= 0x7ff;
    float phaseOffset = *this->initPhase;
    if (phaseOffset < 0.0f) {
        // Negative phase is used as note-on startup delay, not as running phase offset.
        phaseOffset = 0.0f;
    }

    switch ((int)(lfo->freq * 10.0f + .05f)) {
    case LFO_MIDICLOCK_MC_DIV_16:
        if ((songPosition & 0x1)==0) {
            if (computeStep) {
                currentFreq = PREENFM_FREQUENCY / BLOCK_SIZE / 32.0f * invTab[ticks];
                ticks = 0;
            }
            phase = (songPosition & 0x3E) * 0.015625f + phaseOffset;
        }
        break;
    case LFO_MIDICLOCK_MC_DIV_8:
        if ((songPosition & 0x1)==0) {
            if (computeStep) {
                currentFreq = PREENFM_FREQUENCY / BLOCK_SIZE / 16.0f * invTab[ticks];
                ticks = 0;
            }
            phase = (songPosition & 0x1E) * 0.03125f + phaseOffset;
        }
        break;
    case LFO_MIDICLOCK_MC_DIV_4:
        if ((songPosition & 0x1)==0) {
            if (computeStep) {
                currentFreq = PREENFM_FREQUENCY / BLOCK_SIZE / 8.0f * invTab[ticks];
                ticks = 0;
            }
            phase = (songPosition & 0xE) * 0.0625f + phaseOffset;
        }
        break;
    case LFO_MIDICLOCK_MC_DIV_2:
        if ((songPosition & 0x1)==0) {
            if (computeStep) {
                currentFreq = PREENFM_FREQUENCY / BLOCK_SIZE / 4.0f * invTab[ticks];
                ticks = 0;
            }
            // 0,2,4,6
            phase = (songPosition & 0x6) * .125f + phaseOffset;
        }
        break;
    case LFO_MIDICLOCK_MC:
        // Midi Clock
        if ((songPosition & 0x1) == 0) {
            if (computeStep) {
                currentFreq = PREENFM_FREQUENCY / BLOCK_SIZE / 2.0f * invTab[ticks];
                ticks = 0;
            }
            // 0 or 2 -> 0 ou .5
            phase = (songPosition & 0x2) * .25f + phaseOffset;
        }
        break;
    case LFO_MIDICLOCK_MC_TIME_2:
        if ((songPosition & 0x1)==0) {
            if (computeStep) {
                currentFreq = PREENFM_FREQUENCY / BLOCK_SIZE * invTab[ticks];
                ticks = 0;
            }
            phase =  phaseOffset;
        }
        break;
    case LFO_MIDICLOCK_MC_TIME_3:
        if ((songPosition & 0x3)==0) {
            if (computeStep) {
                currentFreq = PREENFM_FREQUENCY / BLOCK_SIZE * invTab[ticks] * 3.0;
                ticks = 0;
            }
            phase = phaseOffset;
        }
        break;
    case LFO_MIDICLOCK_MC_TIME_4:
        if ((songPosition & 0x1)==0) {
            if (computeStep) {
                currentFreq = PREENFM_FREQUENCY / BLOCK_SIZE * invTab[ticks] * 2.0f;
                ticks = 0;
            }
            phase = phaseOffset;
        }
        break;
    case LFO_MIDICLOCK_MC_TIME_8:
        if ((songPosition & 0x1)==0) {
            if (computeStep) {
                currentFreq = PREENFM_FREQUENCY / BLOCK_SIZE * invTab[ticks] * 4.0f;
                ticks = 0;
            }
            phase = phaseOffset;
        }
        break;
    }
}


void LfoOsc::nextValueInMatrix() {
    if (unlikely(oneShotLimit > 0 && !oneShotActive)) {
        matrix->setSource((enum SourceEnum)source, oneShotHoldValue);
        return;
    }

    float lfoValue = 0.0f;

    if (startupDelaySeconds > 0.0f) {
        startupDelaySeconds -= PREENFM_FREQUENCY_INVERSED_LFO;
        if (startupDelaySeconds > 0.0f) {
            (this->*shapeExecutor)(lfoValue, 0.0f, false);
            lfoValue += lfo->bias;
            matrix->setSource((enum SourceEnum)source, lfoValue);
            return;
        }
        startupDelaySeconds = 0.0f;
    }

    ticks ++;

    if (this->isNotMidiSynchronized) {
        currentFreq = lfo->freq + this->matrix->getDestination(destination);
    }
    phase += currentFreq * PREENFM_FREQUENCY_INVERSED_LFO;

    if (unlikely(oneShotActive && phase >= 1.0f)) {
        while (phase >= 1.0f) {
            phase -= 1.0f;

            oneShotRemaining--;
            if (oneShotRemaining <= 0) {
                oneShotActive = false;
                switch (oneShotTerminalMode) {
                case 1:
                    lfoValue = currentRandomValue;
                    break;
                case 2:
                    lfoValue = nextRandomValue;
                    break;
                default:
                    lfoValue = oneShotTerminalShapeValue;
                    break;
                }
                lfoValue += lfo->bias;
                oneShotHoldValue = lfoValue;
                matrix->setSource((enum SourceEnum)source, lfoValue);
                return;
            }

            if (lfo->shape == LFO_RANDOM) {
                currentRandomValue = noise[0];
            } else if (lfo->shape == LFO_BROWNIAN) {
                noiseLp = noise[0] * 0.4f + noiseLp * 0.6f;
                currentRandomValue = noiseLp;
            } else if (lfo->shape == LFO_WANDERING) {
                currentRandomValue = nextRandomValue;
                nextRandomValue = noise[0];
            } else if (lfo->shape == LFO_FLOW) {
                noiseLp = noise[0] * 0.4f + noiseLp * 0.6f;
                currentRandomValue = nextRandomValue;
                nextRandomValue = noiseLp;
            }
        }
    }

    bool phaseWrapped = false;
    if (unlikely(phase >= 1.0f)) {
        // One wrap normalization for all non-one-shot shape paths.
        while (phase >= 1.0f) {
            phase -= 1.0f;
            phaseWrapped = true;
        }
    }

    (this->*shapeExecutor)(lfoValue, phase, phaseWrapped);


    if (unlikely(currentRamp < ramp)) {
        lfoValue = lfoValue * currentRamp  * rampInv ;
        currentRamp += PREENFM_FREQUENCY_INVERSED_LFO;
    }

    lfoValue += lfo->bias;

    matrix->setSource((enum SourceEnum)source, lfoValue);
}


void LfoOsc::noteOn() {
    auto setupNegativeDelay = [this]() {
        if (*this->initPhase < 0.0f) {
            this->phase = 0.0f;
            this->startupDelaySeconds = -*this->initPhase;
            return true;
        }
        this->startupDelaySeconds = 0.0f;
        return false;
    };

    bool hasNegativeDelay = setupNegativeDelay();

    if (oneShotLimit > 0) {
        if (!hasNegativeDelay && (lfo->freq * 10.0f) < LFO_MIDICLOCK_MC_DIV_16) {
            phase = *this->initPhase;
        }

        oneShotRemaining = oneShotLimit;
        oneShotActive = true;
        currentRamp = ramp >= 0.0f ? 0.0f : 1.0f;

        if (unlikely(lfo->shape == LFO_RANDOM)) {
            currentRandomValue = noise[0];
        } else if (unlikely(lfo->shape == LFO_BROWNIAN || lfo->shape == LFO_WANDERING || lfo->shape == LFO_FLOW)) {
            noiseLp = noise[0] * 0.4f + noiseLp * 0.6f;
            currentRandomValue = noiseLp;
        }
        return;
    }

    if (ramp >= 0.0f) {
        currentRamp = 0.0f;
        if (!hasNegativeDelay && (lfo->freq * 10.0f) < LFO_MIDICLOCK_MC_DIV_16) {
            phase = *this->initPhase;
        }
        // Retriger value if random...
        if (unlikely(lfo->shape == LFO_RANDOM)) {
            currentRandomValue = noise[0];
        } else if (unlikely(lfo->shape == LFO_BROWNIAN || lfo->shape == LFO_WANDERING || lfo->shape == LFO_FLOW)) {
            noiseLp = noise[0] * 0.4f + noiseLp * 0.6f;
            currentRandomValue = noiseLp;
        }
    } else {
        // For KSyn Off
        currentRamp = 1; // greater than 0 :
    }
}

