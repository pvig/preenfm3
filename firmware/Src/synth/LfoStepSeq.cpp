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

#include <math.h>
#include "LfoStepSeq.h"
#include "Osc.h"

// Linear step values (12 = 1.0f -> 1 octave / 12 degrees in Note+ with mul = 1.0)
static const float linValues[16] = {
    0.0f / 12.0f, 1.0f / 12.0f, 2.0f / 12.0f, 3.0f / 12.0f,
    4.0f / 12.0f, 5.0f / 12.0f, 6.0f / 12.0f, 7.0f / 12.0f,
    8.0f / 12.0f, 9.0f / 12.0f, 10.0f / 12.0f, 11.0f / 12.0f,
    12.0f / 12.0f, 13.0f / 12.0f, 14.0f / 12.0f, 15.0f / 12.0f
};

void LfoStepSeq::init(struct StepSequencerParams* stepSeqParam, struct StepSequencerSteps *stepSeqSteps, Matrix *matrix, SourceEnum source, DestinationEnum dest, DestinationEnum lenDest) {
	Lfo::init(matrix, source, (DestinationEnum)0);
    this->seqParams = stepSeqParam;
    this->seqSteps = stepSeqSteps;
    this->matrixGateDestination = dest;
    this->matrixLengthDestination = lenDest;
	gated = false;
	stepOnNoteStarted = false;
	loopCount = 0;
	targetOutput = 0.0f;
	currentOutput = 0.0f;
	valueChanged(0);
	ticks = 1536;
	midiClock(0, true);
}

int LfoStepSeq::getActiveSteps(uint8_t activeSteps[16]) const {
	int maxSteps = (int) seqParams->loopLength;
	if (maxSteps < 1 || maxSteps > 16) {
		maxSteps = 16;
	}
	if (this->matrix != NULL && this->matrixLengthDestination != DESTINATION_NONE) {
		maxSteps += (int) roundf(this->matrix->getDestination(this->matrixLengthDestination) * 12.0f);
		if (maxSteps < 1) {
			maxSteps = 1;
		} else if (maxSteps > 16) {
			maxSteps = 16;
		}
	}
	int count = 0;
	for (int k = 0; k < maxSteps; k++) {
		int mode = ((uint8_t) seqSteps->steps[k] >> 4) & 0x0F;
		if (mode != STEPSEQ_MODE_PASS) {
			activeSteps[count++] = (uint8_t) k;
		}
	}
	return count;
}

void LfoStepSeq::midiClock(int songPosition, bool computeStep) {

	ticks &= 0x7ff;
	uint8_t activeSteps[16];
	int loopLen = getActiveSteps(activeSteps);
	if (loopLen == 0) {
		loopLen = 1;
	}

    switch ((int)this->seqParams->bpm)  {
	case LFO_SEQ_MIDICLOCK_DIV_4:
		// Midi Clock  / 4
		if ((songPosition & 0x1)==0) {
			if (computeStep) {
                phaseStep = 0.5f * invTab[ticks];
                ticks = 0;
			}
            loopCount = (songPosition >> 2) / loopLen;
            phase = (float)((songPosition >> 2) % loopLen) + (float)(songPosition & 0x3) * .25f;
		}
		break;
	case LFO_SEQ_MIDICLOCK_DIV_2:
		if ((songPosition & 0x1)==0) {
			if (computeStep) {
                phaseStep = 1.0f * invTab[ticks];
				ticks = 0;
			}
			loopCount = (songPosition >> 1) / loopLen;
			phase = (float)((songPosition >> 1) % loopLen) + (float)(songPosition & 0x1) * .5f;
		}
		break;
	case LFO_SEQ_MIDICLOCK:
		if ((songPosition & 0x1)==0) {
			if (computeStep) {
			    // We're called evey 2 ticks which is half a quarter
			    // We want to move forward 4 beats per quarter so : 2 beats per half a quarter
                phaseStep = 2.0f * invTab[ticks];
				ticks = 0;
			}
			loopCount = songPosition / loopLen;
			phase = (float)(songPosition % loopLen);
		}
		break;
	case LFO_SEQ_MIDICLOCK_TIME_2:
		if ((songPosition & 0x1)==0) {
			if (computeStep) {
                phaseStep = 4.0f * invTab[ticks];
				ticks = 0;
			}
			loopCount = (songPosition << 1) / loopLen;
			phase = (float)((songPosition << 1) % loopLen);
		}
		break;
	case LFO_SEQ_MIDICLOCK_TIME_4:
		if ((songPosition & 0x1)==0) {
			if (computeStep) {
                phaseStep = 8.0f * invTab[ticks];
				ticks = 0;
			}
            loopCount = (songPosition << 1) / loopLen;
            phase = (float)((songPosition << 1) % loopLen);
		}
		break;
	}
}


void LfoStepSeq::valueChanged(int encoder) {
	if (encoder < 2) {
		if (seqParams->bpm < LFO_SEQ_MIDICLOCK_DIV_4) {
		    // bpm / 60 = Herz
		    // We have 4 tick per beat so * 4
		    // 4/60 = .066666
			phaseStep = seqParams->bpm * 0.066666666666666f * PREENFM_FREQUENCY_INVERSED_LFO;
		} else if ((int) seqParams->bpm == LFO_SEQ_STEP_ON_NOTE) {
			phaseStep = 0.0f;
			phase = (float) ((int) phase);
			stepOnNoteStarted = false;
		}
	}
}

void LfoStepSeq::nextValueInMatrix() {

	ticks++;

	phase += phaseStep;
	int phaseInteger = phase;
	phase -= phaseInteger;
	float phaseDecimal = phase;

	uint8_t activeSteps[16];
	int count = getActiveSteps(activeSteps);
	if (count > 0 && phaseInteger >= count) {
		if (seqParams->bpm < LFO_SEQ_MIDICLOCK_DIV_4) {
			loopCount += phaseInteger / count;
		}
		phaseInteger %= count;
	} else if (count == 0) {
		phaseInteger = 0;
	}
	phase += phaseInteger;

	// Add gate and matrix value
	float gatePlusMatrix = seqParams->gate + (this->matrix->getDestination(matrixGateDestination));

	int stepMode = STEPSEQ_MODE_NORM;
	bool isGlide = false;
	if (count == 0) {
		targetOutput = 0.0f;
		gated = true;
	} else {
		uint8_t stepByte = (uint8_t) seqSteps->steps[activeSteps[phaseInteger]];
		int stepVal = stepByte & 0x0F;
		stepMode = (stepByte >> 4) & 0x0F;
		isGlide = (stepMode >= STEPSEQ_MODE_GLIDE1 && stepMode <= STEPSEQ_MODE_GLIDE3);

		bool altMuted = (stepMode == STEPSEQ_MODE_ALT2 && (loopCount & 1) != 0)
			|| (stepMode == STEPSEQ_MODE_ALT3 && (loopCount % 3) != 0);

		if (stepMode == STEPSEQ_MODE_MUTE || altMuted) {
			targetOutput = 0.0f;
			gated = true;
		} else if (isGlide) {
			// Option A (Slide/Tie): keep gate open across glided steps.
			targetOutput = linValues[stepVal];
			gated = false;
		} else if (stepMode >= STEPSEQ_MODE_LFO1 && stepMode <= STEPSEQ_MODE_LFO3) {
			// LFO-controlled trig/ratchet: gate opens when the selected LFO > 0 (and step gate is open).
			bool lfoTrigOpen = this->matrix->getSource((SourceEnum) (MATRIX_SOURCE_LFO1 + (stepMode - STEPSEQ_MODE_LFO1))) > 0.0f;
			bool stepGateOpen = (gatePlusMatrix >= 1.0f) || (gatePlusMatrix > 0.0f && phaseDecimal <= gatePlusMatrix);
			if (lfoTrigOpen && stepGateOpen) {
				targetOutput = linValues[stepVal];
				gated = false;
			} else {
				targetOutput = 0.0f;
				gated = true;
			}
		} else if (stepMode >= STEPSEQ_MODE_RATCHET2 && stepMode <= STEPSEQ_MODE_RATCHET4) {
			// Classic subdivisions per 1/16 step: 2 (1/32), 3 (1/48 triplet), 4 (1/64).
			int subDiv = 2 + (stepMode - STEPSEQ_MODE_RATCHET2);
			float subPhase = phaseDecimal * (float) subDiv;
			subPhase -= (int) subPhase;
			float ratchetGate = (gatePlusMatrix >= 1.0f) ? 0.5f : gatePlusMatrix;
			if (ratchetGate <= 0.0f || subPhase > ratchetGate) {
				targetOutput = 0.0f;
				gated = true;
			} else {
				targetOutput = linValues[stepVal];
				gated = false;
			}
		} else if (gatePlusMatrix <= 0.0f) {
			targetOutput = 0.0f;
		} else if (gatePlusMatrix < 1.0f) {
			if (phaseDecimal > gatePlusMatrix) {
				targetOutput = 0.0f;
				gated = true;
			} else {
				targetOutput = linValues[stepVal];
				gated = false;
			}
		} else {
			targetOutput = linValues[stepVal];
			gated = false;
		}
	}

	if (isGlide) {
		// Exponential portamento scaled by phaseStep (3 speeds: Gld1 fast, Gld2 medium, Gld3 slow).
		static const float glideRates[3] = { 8.0f, 2.6667f, 1.1429f };
		float effectivePhaseStep = (phaseStep > 0.0f) ? phaseStep : (120.0f * 0.066666666666666f * PREENFM_FREQUENCY_INVERSED_LFO);
		float coeff = effectivePhaseStep * glideRates[stepMode - STEPSEQ_MODE_GLIDE1];
		if (coeff > 1.0f) {
			coeff = 1.0f;
		}
		currentOutput += coeff * (targetOutput - currentOutput);
	} else {
		// Fast anti-click ramp on non-glided steps.
		float diff = targetOutput - currentOutput;
		if (diff > 0.1f) {
			currentOutput += 0.1f;
		} else if (diff < -0.1f) {
			currentOutput -= 0.1f;
		} else {
			currentOutput = targetOutput;
		}
	}
	matrix->setSource((enum SourceEnum)source, currentOutput);
}

void LfoStepSeq::noteOn() {
	if (seqParams->bpm < LFO_SEQ_MIDICLOCK_DIV_4) {
		phase = 0;
		loopCount = 0;
		uint8_t activeSteps[16];
		int count = getActiveSteps(activeSteps);
		if (count > 0) {
			uint8_t stepByte = (uint8_t) seqSteps->steps[activeSteps[0]];
			int stepMode = (stepByte >> 4) & 0x0F;
			targetOutput = (stepMode == STEPSEQ_MODE_MUTE) ? 0.0f : linValues[stepByte & 0x0F];
		} else {
			targetOutput = 0.0f;
		}
	} else if ((int) seqParams->bpm == LFO_SEQ_STEP_ON_NOTE) {
		uint8_t activeSteps[16];
		int count = getActiveSteps(activeSteps);
		if (count == 0) {
			phase = 0;
			loopCount = 0;
			targetOutput = 0.0f;
			return;
		}
		if (!stepOnNoteStarted) {
			phase = 0;
			loopCount = 0;
			stepOnNoteStarted = true;
		} else {
			int nextStep = (int) phase + 1;
			if (nextStep >= count) {
				nextStep = 0;
				loopCount++;
			}
			phase = (float) nextStep;
		}
		uint8_t stepByte = (uint8_t) seqSteps->steps[activeSteps[(int) phase]];
		int stepMode = (stepByte >> 4) & 0x0F;
		bool altMuted = (stepMode == STEPSEQ_MODE_ALT2 && (loopCount & 1) != 0)
			|| (stepMode == STEPSEQ_MODE_ALT3 && (loopCount % 3) != 0);
		if (stepMode == STEPSEQ_MODE_MUTE || altMuted) {
			targetOutput = 0.0f;
		} else if (stepMode >= STEPSEQ_MODE_LFO1 && stepMode <= STEPSEQ_MODE_LFO3) {
			bool lfoTrigOpen = this->matrix->getSource((SourceEnum) (MATRIX_SOURCE_LFO1 + (stepMode - STEPSEQ_MODE_LFO1))) > 0.0f;
			targetOutput = lfoTrigOpen ? linValues[stepByte & 0x0F] : 0.0f;
		} else {
			targetOutput = linValues[stepByte & 0x0F];
		}
	}
}

void LfoStepSeq::noteOff() {
}

