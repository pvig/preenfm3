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

#ifndef LFOOSC_H_
#define LFOOSC_H_

#include "Lfo.h"
#include "Osc.h"



class LfoOsc: public Lfo {
public:
    virtual ~LfoOsc() {};

    using Lfo::init;

    void init(struct LfoParams *lfoParams, float* lfoSyncMode, float* lfoPhase, Matrix* matrix, SourceEnum source, DestinationEnum dest);

    inline void updateShapeConstants() {
        // Precompute coefficients used in shape runtime code.
        adAttackPhaseScale = 0.0f;
        adDecayPhaseOffset = 0.0f;
        adDecayPhaseScale = 0.0f;
        ahdAttackPhaseScale = 0.0f;
        ahdDecayPhaseOffset = 0.0f;
        ahdDecayPhaseScale = 0.0f;
        shapeAttack = 0.0f;
        shapeHold = 0.0f;
        shapeInvAttack = 0.0f;
        shapeInvDecay = 0.0f;
        shapeWaveTable = 0;
        shapeWaveMax = 0;
        randomRuntimeMode = 0;

        switch ((int)lfo->shape) {
        case LFO_TRIANGLE:
            shapeExecutor = &LfoOsc::executeTriangle;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = -1.0f;
            break;
        case LFO_SAW:
            shapeExecutor = &LfoOsc::executeSawUp;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = 1.0f;
            break;
        case LFO_SAW_DOWN:
            shapeExecutor = &LfoOsc::executeSawDown;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = -1.0f;
            break;
        case LFO_DECAY_EXP:
            shapeExecutor = &LfoOsc::executeDecayExp;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = -1.0f;
            break;
        case LFO_DECAY_LOG:
            shapeExecutor = &LfoOsc::executeDecayLog;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = -1.0f;
            break;
        case LFO_DECAY_S:
            shapeExecutor = &LfoOsc::executeDecayS;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = -1.0f;
            break;
        case LFO_RISE_EXP:
            shapeExecutor = &LfoOsc::executeRiseExp;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = 1.0f;
            break;
        case LFO_RISE_LOG:
            shapeExecutor = &LfoOsc::executeRiseLog;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = 1.0f;
            break;
        case LFO_ATTACK_DECAY:
            shapeExecutor = &LfoOsc::executeAttackDecay;
            adAttackPhaseScale = 2.0f;
            adDecayPhaseOffset = 0.5f;
            adDecayPhaseScale = 2.0f;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = -1.0f;
            break;
        case LFO_ATTACK_HOLD_DECAY:
            shapeExecutor = &LfoOsc::executeAttackHoldDecay;
            ahdAttackPhaseScale = 4.0f;
            ahdDecayPhaseOffset = 0.5f;
            ahdDecayPhaseScale = 2.0f;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = -1.0f;
            break;
        case LFO_BUCHLA_PLONG:
            shapeExecutor = &LfoOsc::executeBuchla;
            shapeAttack = 0.05f;
            shapeHold = 0.06f;
            shapeInvAttack = 1.0f / shapeAttack;
            shapeInvDecay = 1.0f / (1.0f - shapeAttack - shapeHold);
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = -1.0f;
            break;
        case LFO_BUCHLA_PLONG2:
            shapeExecutor = &LfoOsc::executeBuchla;
            shapeAttack = 0.05f;
            shapeHold = 0.12f;
            shapeInvAttack = 1.0f / shapeAttack;
            shapeInvDecay = 1.0f / (1.0f - shapeAttack - shapeHold);
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = -1.0f;
            break;
        case LFO_SIN:
            shapeExecutor = &LfoOsc::executeWaveTable;
            shapeWaveTable = waveTables[OSC_SHAPE_SIN].table;
            shapeWaveMax = waveTables[OSC_SHAPE_SIN].max;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = 0.0f;
            break;
        case LFO_SIN_SQUARE:
            shapeExecutor = &LfoOsc::executeWaveTable;
            shapeWaveTable = waveTables[OSC_SHAPE_SIN_SQUARE].table;
            shapeWaveMax = waveTables[OSC_SHAPE_SIN_SQUARE].max;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = waveTables[OSC_SHAPE_SIN_SQUARE].table[0];
            break;
        case LFO_SIN_ZERO:
            shapeExecutor = &LfoOsc::executeWaveTable;
            shapeWaveTable = waveTables[OSC_SHAPE_SIN_ZERO].table;
            shapeWaveMax = waveTables[OSC_SHAPE_SIN_ZERO].max;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = waveTables[OSC_SHAPE_SIN_ZERO].table[0];
            break;
        case LFO_SIN_POS:
            shapeExecutor = &LfoOsc::executeWaveTable;
            shapeWaveTable = waveTables[OSC_SHAPE_SIN_POS].table;
            shapeWaveMax = waveTables[OSC_SHAPE_SIN_POS].max;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = waveTables[OSC_SHAPE_SIN_POS].table[0];
            break;
        case LFO_USER1:
            shapeExecutor = &LfoOsc::executeWaveTable;
            shapeWaveTable = waveTables[OSC_SHAPE_USER1].table;
            shapeWaveMax = waveTables[OSC_SHAPE_USER1].max;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = waveTables[OSC_SHAPE_USER1].table[0];
            break;
        case LFO_USER2:
            shapeExecutor = &LfoOsc::executeWaveTable;
            shapeWaveTable = waveTables[OSC_SHAPE_USER2].table;
            shapeWaveMax = waveTables[OSC_SHAPE_USER2].max;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = waveTables[OSC_SHAPE_USER2].table[0];
            break;
        case LFO_USER3:
            shapeExecutor = &LfoOsc::executeWaveTable;
            shapeWaveTable = waveTables[OSC_SHAPE_USER3].table;
            shapeWaveMax = waveTables[OSC_SHAPE_USER3].max;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = waveTables[OSC_SHAPE_USER3].table[0];
            break;
        case LFO_USER4:
            shapeExecutor = &LfoOsc::executeWaveTable;
            shapeWaveTable = waveTables[OSC_SHAPE_USER4].table;
            shapeWaveMax = waveTables[OSC_SHAPE_USER4].max;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = waveTables[OSC_SHAPE_USER4].table[0];
            break;
        case LFO_USER5:
            shapeExecutor = &LfoOsc::executeWaveTable;
            shapeWaveTable = waveTables[OSC_SHAPE_USER5].table;
            shapeWaveMax = waveTables[OSC_SHAPE_USER5].max;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = waveTables[OSC_SHAPE_USER5].table[0];
            break;
        case LFO_USER6:
            shapeExecutor = &LfoOsc::executeWaveTable;
            shapeWaveTable = waveTables[OSC_SHAPE_USER6].table;
            shapeWaveMax = waveTables[OSC_SHAPE_USER6].max;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = waveTables[OSC_SHAPE_USER6].table[0];
            break;
        case LFO_SQUARE:
            shapeExecutor = &LfoOsc::executeSquare;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = 1.0f;
            break;
        case LFO_RANDOM:
            shapeExecutor = &LfoOsc::executeRandomFamily;
            randomRuntimeMode = 1;
            oneShotTerminalMode = 1;
            oneShotTerminalShapeValue = 0.0f;
            break;
        case LFO_BROWNIAN:
            shapeExecutor = &LfoOsc::executeRandomFamily;
            randomRuntimeMode = 2;
            oneShotTerminalMode = 1;
            oneShotTerminalShapeValue = 0.0f;
            break;
        case LFO_WANDERING:
            shapeExecutor = &LfoOsc::executeRandomFamily;
            randomRuntimeMode = 3;
            oneShotTerminalMode = 2;
            oneShotTerminalShapeValue = 0.0f;
            break;
        case LFO_FLOW:
            shapeExecutor = &LfoOsc::executeRandomFamily;
            randomRuntimeMode = 4;
            oneShotTerminalMode = 2;
            oneShotTerminalShapeValue = 0.0f;
            break;
        default:
            shapeExecutor = &LfoOsc::executeTriangle;
            oneShotTerminalMode = 0;
            oneShotTerminalShapeValue = -1.0f;
            break;
        }
    }

	void valueChanged(int encoder) {
	    switch (encoder) {
        case ENCODER_LFO_SHAPE:
            updateShapeConstants();
            break;
        case ENCODER_LFO_KSYNC: {
            this->ramp = lfo->keybRamp;

            int mode = (int)(*syncMode + 0.5f);
            if (mode >= LFO_SYNC_ONESHOT_INTERNAL_1 && mode <= LFO_SYNC_ONESHOT_INTERNAL_8) {
                oneShotLimit = mode - LFO_SYNC_ONESHOT_INTERNAL_1 + 1;
            } else if (mode >= LFO_SYNC_ONESHOT_EXTERNAL_1 && mode <= LFO_SYNC_ONESHOT_EXTERNAL_8) {
                oneShotLimit = mode - LFO_SYNC_ONESHOT_EXTERNAL_1 + 1;
            } else {
                oneShotLimit = 0;
            }

            float keybRampAbs = this->ramp < 0.0f ? -this->ramp : this->ramp;
            this->rampInv = 50 * invTab[(int)(keybRampAbs * 50.0f)];

            if (this->ramp < 0 && oneShotLimit == 0) {
                // resync all LFO
                phase = 0;
            }
            break;
        }
	    case ENCODER_LFO_FREQ:
	        isNotMidiSynchronized = ((lfo->freq * 10.0f) < LFO_MIDICLOCK_MC_DIV_16);
	        break;
	    }
	}


	void midiClock(int songPosition, bool computeStep);

	void nextValueInMatrix();

	void noteOn();

	void noteOff() {
		// Nothing to do
	}



private:
    typedef void (LfoOsc::*ShapeExecutorFn)(float& lfoValue, float phase, bool phaseWrapped);

    void executeTriangle(float& lfoValue, float phase, bool phaseWrapped);
    void executeSawUp(float& lfoValue, float phase, bool phaseWrapped);
    void executeSawDown(float& lfoValue, float phase, bool phaseWrapped);
    void executeDecayExp(float& lfoValue, float phase, bool phaseWrapped);
    void executeDecayLog(float& lfoValue, float phase, bool phaseWrapped);
    void executeRiseExp(float& lfoValue, float phase, bool phaseWrapped);
    void executeRiseLog(float& lfoValue, float phase, bool phaseWrapped);
    void executeAttackDecay(float& lfoValue, float phase, bool phaseWrapped);
    void executeAttackHoldDecay(float& lfoValue, float phase, bool phaseWrapped);
    void executeDecayS(float& lfoValue, float phase, bool phaseWrapped);
    void executeBuchla(float& lfoValue, float phase, bool phaseWrapped);
    void executeWaveTable(float& lfoValue, float phase, bool phaseWrapped);
    void executeSquare(float& lfoValue, float phase, bool phaseWrapped);
    void executeRandomFamily(float& lfoValue, float phase, bool phaseWrapped);

	LfoType type;
	LfoParams* lfo ;
    float* syncMode;
    float currentRamp, ramp, rampInv;
    float phase;
    float* initPhase;
    DestinationEnum destination;
    float currentRandomValue;
    float nextRandomValue = 0;
    float noiseLp = 0;
    float currentFreq ;
    int oneShotLimit = 0;
    int oneShotRemaining = 0;
    bool oneShotActive = false;
    float oneShotHoldValue = 0.0f;
    float shapeAttack = 0.0f;
    float shapeHold = 0.0f;
    float shapeInvAttack = 0.0f;
    float shapeInvDecay = 0.0f;
    int oneShotTerminalMode = 0;
    float oneShotTerminalShapeValue = -1.0f;
    float adAttackPhaseScale = 0.0f;
    float adDecayPhaseOffset = 0.0f;
    float adDecayPhaseScale = 0.0f;
    float ahdAttackPhaseScale = 0.0f;
    float ahdDecayPhaseOffset = 0.0f;
    float ahdDecayPhaseScale = 0.0f;
    float startupDelaySeconds = 0.0f;
    const float* shapeWaveTable = 0;
    int shapeWaveMax = 0;
    int randomRuntimeMode = 0;
    ShapeExecutorFn shapeExecutor = &LfoOsc::executeTriangle;
    //
    bool isNotMidiSynchronized;

};

#endif /* LFOOSC_H_ */
