/*
 * Copyright 2013 Xavier Hosxe
 *
 * Author: Xavier Hosxe (xavier <dot> hosxe (at) g m a i l <dot> com)
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

#include "stm32h7xx_hal.h"
#include "FMDisplayMixer.h"
#include "FMDisplayMenu.h"
#include "FMDisplayEditor.h"
#include "FMDisplaySequencer.h"
#include "SynthState.h"
#include "Hexter.h"
#include "Timbre.h"
#include "Synth.h"
#include "preenfm3lib.h"
#include "preenfm3.h"

extern Synth synth;
extern RNG_HandleTypeDef hrng;
extern float diatonicScaleFrequency[];


SynthState::SynthState() {
    operatorNumber = 0;

    // First default preset
    fullState.synthMode = SYNTH_MODE_MIXER;
    fullState.preenFMBankNumber = 0;
    fullState.preenFMPresetNumber = 0;
    fullState.preenFMMixerNumber = 0;
    fullState.preenFMMixerPresetNumber = 0;
    fullState.dx7BankNumber = 0;
    fullState.dx7PresetNumber = 0;
    // Default pfm3 edit page
    fullState.mainPage = -1;
    fullState.editPage = 0;

    // https://en.wikipedia.org/wiki/C_(musical_note)
    // Frequency of note 60 C4 :
    fullState.midiConfigValue[MIDICONFIG_USB] = 2;
    fullState.midiConfigValue[MIDICONFIG_RECEIVES] = 3;
    fullState.midiConfigValue[MIDICONFIG_SENDS] = 1;
    fullState.midiConfigValue[MIDICONFIG_PROGRAM_CHANGE] = 1;
    fullState.midiConfigValue[MIDICONFIG_TEST_NOTE] = 60;
    fullState.midiConfigValue[MIDICONFIG_TEST_VELOCITY] = 120;
    fullState.midiConfigValue[MIDICONFIG_ENCODER] = 1;
    fullState.midiConfigValue[MIDICONFIG_ARPEGGIATOR_IN_PRESET] = 1;
    fullState.midiConfigValue[MIDICONFIG_CPU_USAGE] = 0;
    fullState.midiConfigValue[MIDICONFIG_TFT_BACKLIGHT] = 100;
    fullState.midiConfigValue[MIDICONFIG_TFT_AUTO_REINIT] = 0;
    fullState.midiConfigValue[MIDICONFIG_ENCODER_PUSH] = 0;
    fullState.midiConfigValue[MIDICONFIG_REVERB_PARAMS] = 0;
    // Init randomizer values to 1
    fullState.randomizer.Oper = 1;
    fullState.randomizer.EnvT = 1;
    fullState.randomizer.IM = 1;
    fullState.randomizer.Modl = 1;

    // Mixer
    fullState.mixerCurrentEdit = 0;
    fullState.menuCurrentEdit = 0;

    for (int k = 0; k < 12; k++) {
        fullState.name[k] = 0;
    }

    // edit with timbre 0
    currentTimbre = 0;
    stepSelect[0] = 0;
    stepSelect[1] = 0;
    patternSelect = 0;

    isPlayingNote = false;

    for (int row = 0; row < NUMBER_OF_ROWS; row++) {
        for (int param = 0; param < NUMBER_OF_ENCODERS_PFM2; param++) {
            struct ParameterDisplay* paramDisplay = &(allParameterRows.row[row]->params[param]);
            if (paramDisplay->numberOfValues > 1.0) {
                paramDisplay->incValue = ((paramDisplay->maxValue - paramDisplay->minValue) / (paramDisplay->numberOfValues - 1.0f));
            } else {
                paramDisplay->incValue = 0.0f;
            }
        }
    }

    // Init mixer state with default values
    char mixerStateChars[sizeof(mixerState)];
    uint32_t size;
    mixerState.getFullDefaultState(mixerStateChars, &size, 0);
    mixerState.restoreFullState(mixerStateChars);

    for (int b = 0; b < NUMBER_OF_BUTTONIDS; b++) {
        fullState.buttonState[b] = 0;
    }
    fullState.operatorNumber = 0;
}

void SynthState::init(FMDisplayMixer* displayMixer, FMDisplayEditor* displayEditor, FMDisplayMenu* displayMenu, FMDisplaySequencer* displaySequencer) {
    this->displayMixer = displayMixer;
    this->displayEditor = displayEditor;
    this->displayMenu = displayMenu;
    this->displaySequencer = displaySequencer;

    mixerState.fxBus_.init();
}


void SynthState::encoderTurnedForStepSequencer(int row, int encoder4, int encoder6, int ticks) {
    int whichStepSeq = row - ROW_LFOSEQ1;
    StepSequencerSteps * seqSteps = &((StepSequencerSteps *) (&params->lfoSteps1))[whichStepSeq];

    if (encoder6 == 3) {
        int oldPos = stepSelect[whichStepSeq];
        int newValue = stepSelect[whichStepSeq] + (ticks > 0 ? 1 : -1);

        if (newValue > 15) {
            newValue = 0;
        } else if (newValue < 0) {
            newValue = 15;
        }

        stepSelect[whichStepSeq] = newValue;
        propagateNewParamValue(currentTimbre, row, encoder4, (ParameterDisplay*) NULL, oldPos, stepSelect[whichStepSeq]);

    } else if (encoder6 == 4) {
        char *step = &seqSteps->steps[(int) stepSelect[whichStepSeq]];
        int oldValue = (int) (*step);

        int newValue = (*step) + ticks;

        if (newValue > 15) {
            newValue = 15;
        } else if (newValue < 0) {
            newValue = 0;
        }
        if ((*step) != newValue) {
            (*step) = newValue;
            propagateNewParamValue(currentTimbre, row, encoder4, (ParameterDisplay*) NULL, oldValue, (int) (*step));
        }
    }
}

void SynthState::encoderTurnedForArpPattern(int row, int encoder, int ticks) {
    if (encoder == 0) {
        // Encoder 0: move cursor
        int oldPos = patternSelect;
        int newPos = oldPos;
        newPos += (ticks > 0 ? 1 : -1);

        if (newPos > 15) {
            newPos = 0;
        } else if (newPos < 0) {
            newPos = 15;
        }
        if (newPos != oldPos) {
            patternSelect = newPos;
            propagateNewParamValue(currentTimbre, row, encoder, (ParameterDisplay*) NULL, oldPos, patternSelect);
        }

    } else if (encoder == 1) {
        // Change value(s)
        if (params->engineArp2.pattern < ARPEGGIATOR_PRESET_PATTERN_COUNT) {
            return;
        }
        arp_pattern_t pattern = params->engineArpUserPatterns.patterns[(int) params->engineArp2.pattern - ARPEGGIATOR_PRESET_PATTERN_COUNT];
        const uint16_t oldMask = ARP_PATTERN_GETMASK(pattern);
        uint16_t newMask = oldMask;

        uint16_t bitsToModify = 0;
        switch (encoder) {
        case 1:
            bitsToModify = 0x1 << patternSelect;
            break;       // modify single note
        case 2:
            bitsToModify = 0x1111 << (patternSelect & 3);
            break; // modify all
        case 3:
            bitsToModify = 0xf << ((patternSelect >> 2) << 2);
            break; // modify entire bar
        }
        if (ticks > 0) {
            newMask |= bitsToModify;
        } else {
            newMask &= ~bitsToModify;
        }

        if (oldMask != newMask) {
            ARP_PATTERN_SETMASK(pattern, newMask);
            params->engineArpUserPatterns.patterns[(int) params->engineArp2.pattern - ARPEGGIATOR_PRESET_PATTERN_COUNT] = pattern;
            propagateNewParamValue(currentTimbre, row, encoder, (ParameterDisplay*) NULL, oldMask, newMask);
        }
    }
}

void SynthState::twoButtonsPressed(int button1, int button2) {
    switch (button1) {
    case BUTTON_PFM3_MENU:
        switch (button2) {
            case BUTTON_PFM3_1:
                propagateNoteOn(-2);
                break;
            case BUTTON_PFM3_2:
                propagateNoteOn(0);
                break;
            case BUTTON_PFM3_3:
                propagateNoteOn(2);
                break;
            case BUTTON_PFM3_4:
                propagateNoteOn(-14);
                break;
            case BUTTON_PFM3_5:
                propagateNoteOn(-12);
                break;
            case BUTTON_PFM3_6:
                propagateNoteOn(-10);
                break;
            case BUTTON_PREVIOUS_INSTRUMENT: {
                SynthEditMode previousMode = fullState.synthMode;
                fullState.synthMode = SYNTH_MODE_REINIT_TFT;
                propagateNewPfm3Page();
                fullState.synthMode = previousMode;
                propagateNewPfm3Page();
                break;
            }
            case BUTTON_PFM3_SEQUENCER:
                preenfm3SwitchToMidiController();
                break;
            case BUTTON_PFM3_MIXER:
                break;
        }
        break;
    case BUTTON_NEXT_INSTRUMENT:
        if (button2 == BUTTON_PREVIOUS_INSTRUMENT) {
            propagateNoteOff();
            for (int t = 0; t < NUMBER_OF_TIMBRES; t++) {
                propagateBeforeNewParamsLoad(t);
            }
            propagateAfterNewMixerLoad();
        }
        break;
    case BUTTON_PREVIOUS_INSTRUMENT:
        if (button2 == BUTTON_NEXT_INSTRUMENT) {
            propagateNoteOff();
            for (int t = 0; t < NUMBER_OF_TIMBRES; t++) {
                propagateBeforeNewParamsLoad(t);
            }
            propagateAfterNewMixerLoad();
        }
        break;
    }

#ifdef PPMIMAGE_ENABLE
    // Screenshot !!
    if (button1 == BUTTON_PFM3_MENU && button2 == BUTTON_NEXT_INSTRUMENT) {
        storage->getPPMImage()->saveImage();
        propagateNewPfm3Page();
    }
#endif

#ifdef DEBUG__KO
    if (button1 == BUTTON_LFO) {
        if (button2 == BUTTON_MATRIX) {
            synth.debugVoice();
        }
        if (button2 == BUTTON_BACK) {
            synth.showCycles();
        }
        if (button2 == BUTTON_MENUSELECT) {
            storage->getPatchBank()->testMemoryPreset();
        }
    }
#endif

}

void SynthState::encoderTurnedWhileButtonPressed(int encoder, int ticks, int button) {
    if (button == BUTTON_NEXT_INSTRUMENT) {
        encoderTurned(encoder, ticks * 10);
        return;
    }

    if (fullState.synthMode == SYNTH_MODE_EDIT_PFM3) {
        displayEditor->encoderTurnedWhileButtonPressed(encoder, ticks, button);
    } else if (fullState.synthMode == SYNTH_MODE_MIXER) {
        displayMixer->encoderTurnedWhileButtonPressed(encoder, ticks, button);
    } else if (fullState.synthMode == SYNTH_MODE_SEQUENCER) {
        displaySequencer->encoderTurnedWhileButtonPressed(encoder, ticks, button);
    }

}

bool SynthState::newRandomizerValue(int encoder, int ticks) {
    int8_t oldValue = 0;
    int8_t newValue = 0;
    switch (encoder) {
    case 0:
        oldValue = fullState.randomizer.Oper;
        fullState.randomizer.Oper += (ticks > 0 ? 1 : -1);
        fullState.randomizer.Oper = fullState.randomizer.Oper > 3 ? 3 : fullState.randomizer.Oper;
        fullState.randomizer.Oper = fullState.randomizer.Oper < 0 ? 0 : fullState.randomizer.Oper;
        newValue = fullState.randomizer.Oper;
        break;
    case 1:
    case 2:
        oldValue = fullState.randomizer.EnvT;
        fullState.randomizer.EnvT += (ticks > 0 ? 1 : -1);
        fullState.randomizer.EnvT = fullState.randomizer.EnvT > 3 ? 3 : fullState.randomizer.EnvT;
        fullState.randomizer.EnvT = fullState.randomizer.EnvT < 0 ? 0 : fullState.randomizer.EnvT;
        newValue = fullState.randomizer.EnvT;
        break;
    case 3:
        oldValue = fullState.randomizer.IM;
        fullState.randomizer.IM += (ticks > 0 ? 1 : -1);
        fullState.randomizer.IM = fullState.randomizer.IM > 3 ? 3 : fullState.randomizer.IM;
        fullState.randomizer.IM = fullState.randomizer.IM < 0 ? 0 : fullState.randomizer.IM;
        newValue = fullState.randomizer.IM;
        break;
    case 4:
    case 5:
        oldValue = fullState.randomizer.Modl;
        fullState.randomizer.Modl += (ticks > 0 ? 1 : -1);
        fullState.randomizer.Modl = fullState.randomizer.Modl > 3 ? 3 : fullState.randomizer.Modl;
        fullState.randomizer.Modl = fullState.randomizer.Modl < 0 ? 0 : fullState.randomizer.Modl;
        newValue = fullState.randomizer.Modl;
        break;
    }
    return newValue != oldValue;
}

void SynthState::encoderTurned(int encoder, int ticks) {
    switch (fullState.synthMode) {
    case SYNTH_MODE_EDIT_PFM3:
        displayEditor->encoderTurnedPfm3(encoder, ticks);
        break;
    case SYNTH_MODE_MENU:
        displayMenu->encoderTurned(currentTimbre, encoder, ticks);
        // Did we modify reverb visibility ?
        if (unlikely(fullState.currentMenuItem->menuState == MENU_CONFIG_SETTINGS && encoder == 3 &&  fullState.menuSelect == MIDICONFIG_REVERB_PARAMS)) {
            displayMixer->setReverbParamVisible(fullState.midiConfigValue[MIDICONFIG_REVERB_PARAMS] > 0);
        }
        break;
    case SYNTH_MODE_MIXER:
        displayMixer->encoderTurned(encoder, ticks);
        break;
    case SYNTH_MODE_SEQUENCER: {
        if (encoder == 2) {
            setCurrentInstrument(ticks > 0 ? 0 : -1);
        }
        displaySequencer->encoderTurned(currentTimbre, encoder, ticks);
        break;
    }
    default:
        break;
    }
}

void SynthState::loadNewPreset(int timbre) {
    storeTestNote();
    propagateNoteOff();
    propagateBeforeNewParamsLoad(timbre);
    storage->getPatchBank()->copyNewPreset(params);
    propagateAfterNewParamsLoad(timbre);
    restoreTestNote();
}

void SynthState::loadPreset(int timbre, PFM3File const *bank, int patchNumber, struct OneSynthParams* params) {
    storeTestNote();
    propagateNoteOff();
    propagateBeforeNewParamsLoad(timbre);
    storage->getPatchBank()->loadPatch(bank, patchNumber, params);
    propagateAfterNewParamsLoad(timbre);
    restoreTestNote();
}

void SynthState::loadDx7Patch(int timbre, PFM3File const *bank, int patchNumber, struct OneSynthParams* params) {
    storeTestNote();
    propagateNoteOff();
    propagateBeforeNewParamsLoad(timbre);
    uint8_t* packedPatch = storage->getDX7SysexFile()->dx7LoadPatch(bank, patchNumber);
    if (packedPatch != 0) {
        hexter->loadHexterPatch(packedPatch, params);
    }
    propagateAfterNewParamsLoad(timbre);
    restoreTestNote();
}

void SynthState::loadMixer(PFM3File const *bank, int patchNumber) {
    propagateBeforeNewParamsLoad(currentTimbre);
    storage->getMixerBank()->loadMixer(bank, patchNumber);
    // Update and clean all timbres
    this->currentTimbre = 0;
    propagateNewTimbre(currentTimbre);
    propagateAfterNewMixerLoad();
}

void SynthState::loadPresetFromMidi(int timbre, int bank, int bankLSB, int patchNumber, struct OneSynthParams *params) {
    switch (bank) {
        case 0: {
            PFM3File const *bank = storage->getPatchBank()->getFile(bankLSB);
            if (bank->fileType != FILE_EMPTY) {
                fullState.preenFMBankNumber = bankLSB;
                fullState.preenFMPresetNumber = patchNumber;
                loadPreset(timbre, bank, patchNumber, params);
            }
            break;
        }
        case 1: {
            PFM3File const *bank = storage->getMixerBank()->getFile(bankLSB);
            if (bank->fileType != FILE_EMPTY) {
                fullState.preenFMMixerNumber = bankLSB;
                fullState.preenFMMixerPresetNumber = patchNumber;
                loadMixer(bank, patchNumber);
            }
            break;
        }
        case 2:
        case 3:
        case 4: {
            int dx7bank = bank - 2;
            PFM3File const *bank = storage->getDX7SysexFile()->getFile(bankLSB + dx7bank * 128);
            if (bank->fileType != FILE_EMPTY) {
                fullState.dx7BankNumber = bankLSB;
                fullState.dx7PresetNumber = patchNumber;
                loadDx7Patch(timbre, bank, patchNumber, params);
            }
            break;
        }
    }
}



void SynthState::setCurrentInstrument(int value) {
    if (value == -1) {
        currentTimbre = currentTimbre - 1;
        if (unlikely(currentTimbre < 0)) {
            currentTimbre = NUMBER_OF_TIMBRES - 1;
        }
    } else if (value == 0) {
        currentTimbre = (currentTimbre + 1) % NUMBER_OF_TIMBRES;
    } else if (value <= (NUMBER_OF_TIMBRES + 1)) {
        currentTimbre = value - 1;
    } else {
        return;
    }
    propagateNewTimbre(currentTimbre);
}

void SynthState::buttonLongPressed(int button) {
    switch (fullState.synthMode) {
        case SYNTH_MODE_SEQUENCER:
            displaySequencer->buttonLongPressed(currentTimbre, button);
            break;
        case SYNTH_MODE_EDIT_PFM3:
            displayEditor->buttonLongPressed(currentTimbre, button);
            break;
        default:
            break;
    }
}

void SynthState::buttonPressed(int button) {
    SynthEditMode oldSynthMode = fullState.synthMode;

    bool swithToMenuIfBackButtonPressed = true;

    if (fullState.synthMode == SYNTH_MODE_EDIT_PFM3) {
        if (fullState.mainPage != -1) {
            // Back must not go to MENU mode if not in the main page
            swithToMenuIfBackButtonPressed = false;
        }
        displayEditor->buttonPressed(button);
    } else if (fullState.synthMode == SYNTH_MODE_MENU) {
        displayMenu->buttonPressed(currentTimbre, button);
        swithToMenuIfBackButtonPressed = false;
    } else if (fullState.synthMode == SYNTH_MODE_MIXER) {
        displayMixer->buttonPressed(button);
    } else if (fullState.synthMode == SYNTH_MODE_SEQUENCER) {
        if (displaySequencer->getSequencerMode() != SEQ_MODE_NORMAL
            && displaySequencer->getSequencerMode() != SEQ_MODE_STEP) {
            swithToMenuIfBackButtonPressed = false;
        }
        displaySequencer->buttonPressed(currentTimbre, button);
    }

    // Anywhere these are main functions !!!!
    switch (button) {
    case BUTTON_PFM3_EDIT:
        fullState.synthMode = SYNTH_MODE_EDIT_PFM3;
        break;
    case BUTTON_PFM3_SEQUENCER:
        fullState.synthMode = SYNTH_MODE_SEQUENCER;
        break;
    case BUTTON_PFM3_MIXER:
        fullState.synthMode = SYNTH_MODE_MIXER;
        break;
    case BUTTON_PFM3_MENU:
        if (swithToMenuIfBackButtonPressed) {
            if (fullState.synthModeBeforeMenu != SYNTH_MODE_MENU) {
                fullState.synthModeBeforeMenu = fullState.synthMode;
            }
            // Make sure to propagate synth mode
            oldSynthMode = SYNTH_MODE_EDIT_PFM3;
            fullState.synthMode = SYNTH_MODE_MENU;
            fullState.previousChoice = fullState.previousMenuChoice.main;
            fullState.currentMenuItem = MenuItemUtil::getMenuItem(MAIN_MENU);
            propagateNewSynthMode();
            return;
        }
        break;
    case BUTTON_PREVIOUS_INSTRUMENT:
        // select next instrument as current one
        if (fullState.synthMode != SYNTH_MODE_MENU) {
            setCurrentInstrument(-1);
        }
        break;
    case BUTTON_NEXT_INSTRUMENT:
        // select next instrument as current one
        if (fullState.synthMode != SYNTH_MODE_MENU) {
            setCurrentInstrument(0);
        }
        break;
    case BUTTON_ENCODER_1:
    case BUTTON_ENCODER_2:
    case BUTTON_ENCODER_3:
    case BUTTON_ENCODER_4:
    case BUTTON_ENCODER_5:
    case BUTTON_ENCODER_6:
        setCurrentInstrument(button - BUTTON_ENCODER_1 + 1);
        break;
    }

    if (oldSynthMode != fullState.synthMode) {
        propagateNewSynthMode();
        return;
    }
}


void SynthState::propagateAfterNewParamsLoad(int timbre) {
    for (SynthParamListener* listener = firstParamListener; listener != 0; listener = listener->nextListener) {
        listener->afterNewParamsLoad(timbre);
    }
}

void SynthState::propagateAfterNewMixerLoad() {
    for (SynthParamListener* listener = firstParamListener; listener != 0; listener = listener->nextListener) {
        listener->afterNewMixerLoad();
    }
}

void SynthState::propagateNewTimbre(int timbre) {
    propagateNoteOff();
    for (SynthParamListener* listener = firstParamListener; listener != 0; listener = listener->nextListener) {
        listener->newTimbre(timbre);
    }
}

void SynthState::tempoClick() {
    if (fullState.synthMode == SYNTH_MODE_MENU) {
        if (fullState.currentMenuItem->menuType == MENUTYPE_TEMPORARY) {
            if (doneClick > 4) {
                fullState.synthMode = fullState.synthModeBeforeMenu;
                propagateNewSynthMode();
            }
            doneClick++;
        }
    } else {
        doneClick = 0;
    }
}

void SynthState::setParamsAndTimbre(struct OneSynthParams *newParams, int newCurrentTimbre) {
    this->params = newParams;
    this->currentTimbre = newCurrentTimbre;
}




/*
 * Randomizer
 */

int getRandomInt(int max) {
    if (max <= 1) {
        return 0;
    }
    uint32_t rnd;
    if (HAL_RNG_GenerateRandomNumber(&hrng, &rnd) != HAL_OK) {
        // Keep randomizer responsive even if hardware RNG is temporarily unavailable.
        rnd = (HAL_GetTick() * 214013u) + 2531011u;
    }
    return rnd % (uint32_t) max;
}

float getRandomFloat(float min, float max) {
    uint32_t rnd;
    if (HAL_RNG_GenerateRandomNumber(&hrng, &rnd) != HAL_OK) {
        rnd = (HAL_GetTick() * 1103515245u) + 12345u;
    }
    float f = ((float) (rnd % 100000)) / 100000.0f;
    return f * (max - min) + min;
}

/*
 * Returns a random oscillator waveform shape index based on Oper level.
 *   soft(1): 0..6, forced to 0 if >2     → sine / triangle / sawtooth only
 *   medi(2): 0..7, skips index 6 (noise) → varied palette, no pure random wave
 *   high(3): 0..6, unrestricted           → full range
 * Pad EnvT additionally re-rolls shapes >2 toward 0..2 with 75% probability.
 */
float getRandomShape(int operatorRandom) {
    int shape = getRandomInt(7);
    switch (operatorRandom) {
    case 1:
        if (shape > 2) {
            shape = 0;
        }
        break;
    case 2:
        shape = getRandomInt(8);
        if (shape == 6) { // Rand
            shape = 0;
        }
        break;
    case 3:
        break;
    }
    return shape;
}

/*
 * Returns frequency type: 0 = keyboard-tracked, 1 = fixed pitch.
 * Probability of keyboard tracking by Oper level:
 *   soft(1): 7/8 = 87.5%   medi(2): 6/8 = 75%   high(3): 4/8 = 50%
 * Higher Oper allows more fixed-frequency modulators for richer FM texture.
 */
float getRandomFrequencyType(int operatorRandom) {
    int keyboardRatio = 7;
    switch (operatorRandom) {
    case 1:
        keyboardRatio = 7; // 7/8 keyboard
        break;
    case 2:
        keyboardRatio = 6; // 6/8 keyboard
        break;
    case 3:
        keyboardRatio = 4; // 4/8 keyboard
        break;
    }
    int freqType = (getRandomInt(8) < keyboardRatio) ? 0 : 1;
    return freqType;
}

/*
 * Returns a frequency multiplier for the oscillator based on Oper level.
 *   soft(1): {0.5, 1, 2, 4}              octave steps only (4 choices)
 *   medi(2): {0.25, 0.5, 1, 1.5, 2, 3, 4}  adds 5th (1.5) and 3rd (3) (7 choices)
 *   high(3): 0.25 to 6.25 in 0.25 steps  full inharmonic palette (24 choices)
 */
float getRandomFrequency(int operatorRandom) {
    float random1Frequency[] = { .5f, 1.0f, 2.0f, 4.0f };
    float random2Frequency[] = { .25, .5f, 1.0f, 1.5, 2.0f, 3.0f, 4.0f };
    float freq = 0;
    switch (operatorRandom) {
    case 1:
        freq = random1Frequency[getRandomInt(4)];
        break;
    case 2:
        freq = random2Frequency[getRandomInt(7)];
        break;
    case 3:
        freq = getRandomInt(24) * .25 + .25;
        break;
    }
    return freq;
}

/*
 * Returns a detune offset per Oper level (soft always returns 0).
 * Values outside the inner band are zeroed so most operators land near unison:
 *   medi(2): sample −0.05..+0.04 in 0.01 steps; keep only if within ±0.03
 *   high(3): sample −0.05..+0.14 in 0.01 steps; keep only if within ±0.07
 */
float getFineTune(int operatorRandom) {
    float fineTune = 0;
    switch (operatorRandom) {
    case 2:
        fineTune = getRandomInt(10) * .01 - .05;
        if (fineTune < -0.03f || fineTune > 0.03f) {
            fineTune = 0;
        }
        break;
    case 3:
        fineTune = getRandomInt(20) * .01 - .05;
        if (fineTune < -0.07f || fineTune > 0.07f) {
            fineTune = 0;
        }
        break;
    }
    return fineTune;
}

/*
 * Pick a matrix modulation destination appropriate for the given Modl level.
 *
 * safeDestinations (24):     IM indices, pan, mix, filter params, env times,
 *                             matrix multipliers, LFO frequencies.
 * advancedDestinations (16):  per-oscillator pitch, phase, warp, feedback —
 *                             unlocked at Modl=3 with 1/3 probability.
 *
 * Destinations targeting inactive oscillators (mix + weighted IM < 0.03) are
 * skipped via a wrapping scan; INDEX_ALL_MODULATION is the final fallback.
 */
DestinationEnum getRandomModDestination(int modulationRandom, const OneSynthParams* params) {
    // Pick oscillator destinations only if that oscillator is currently audible.
    auto absf = [](float v) { return v < 0.0f ? -v : v; };

    auto getOscContribution = [absf](const OneSynthParams* params, int osc) {
        float im1 = params->engineIm1.modulationIndex1 + 0.25f * params->engineIm1.modulationIndexVelo1;
        float im2 = params->engineIm1.modulationIndex2 + 0.25f * params->engineIm1.modulationIndexVelo2;
        float im3 = params->engineIm2.modulationIndex3 + 0.25f * params->engineIm2.modulationIndexVelo3;
        float im4 = params->engineIm2.modulationIndex4 + 0.25f * params->engineIm2.modulationIndexVelo4;
        float im5 = params->engineIm3.modulationIndex5 + 0.25f * params->engineIm3.modulationIndexVelo5;
        float feedback = params->engineIm3.modulationIndex6 + 0.25f * params->engineIm3.modulationIndexVelo6;

        float modContribution = 0.0f;
        switch (osc) {
        case 1:
            modContribution = 0.0f;
            return absf(params->engineMix1.mixOsc1) + modContribution;
        case 2:
            modContribution = absf(im1);
            return absf(params->engineMix1.mixOsc2) + modContribution;
        case 3:
            modContribution = absf(im2);
            return absf(params->engineMix2.mixOsc3) + modContribution;
        case 4:
            modContribution = absf(im3);
            return absf(params->engineMix2.mixOsc4) + modContribution;
        case 5:
            modContribution = absf(im4);
            return absf(params->engineMix3.mixOsc5) + modContribution;
        case 6:
            modContribution = absf(im5) + (0.5f * absf(feedback));
            return absf(params->engineMix3.mixOsc6) + modContribution;
        default: return 0.0f;
        }
    };

    auto isDestinationRelevant = [&](DestinationEnum dest, const OneSynthParams* params) {
        bool osc1Active = getOscContribution(params, 1) > 0.03f;
        bool osc2Active = getOscContribution(params, 2) > 0.03f;
        bool osc3Active = getOscContribution(params, 3) > 0.03f;
        bool osc4Active = getOscContribution(params, 4) > 0.03f;
        bool osc5Active = getOscContribution(params, 5) > 0.03f;
        bool osc6Active = getOscContribution(params, 6) > 0.03f;
        bool anyOscActive = osc1Active || osc2Active || osc3Active || osc4Active || osc5Active || osc6Active;

        switch (dest) {
        case OSC1_FREQ:
        case PAN_OSC1:
        case MIX_OSC1:
        case OSC1_PHASE:
        case OSC1_WARP:
            return osc1Active;
        case OSC2_FREQ:
        case PAN_OSC2:
        case MIX_OSC2:
        case OSC2_PHASE:
        case OSC2_WARP:
            return osc2Active;
        case OSC3_FREQ:
        case PAN_OSC3:
        case MIX_OSC3:
        case OSC3_PHASE:
        case OSC3_WARP:
            return osc3Active;
        case OSC4_FREQ:
        case PAN_OSC4:
        case MIX_OSC4:
        case OSC4_PHASE:
        case OSC4_WARP:
            return osc4Active;
        case OSC5_FREQ:
        case OSC5_PHASE:
        case OSC5_WARP:
            return osc5Active;
        case OSC6_FREQ:
        case OSC6_PHASE:
        case OSC6_WARP:
            return osc6Active;
        case ALL_OSC_FREQ:
        case ALL_OSC_FREQ_HARM:
        case ALL_PAN:
        case ALL_MIX:
            return anyOscActive;
        case FILTER1_PARAM1:
        case FILTER1_PARAM2:
        case FILTER1_AMP:
            // Only route to effect1 params when a filter is actually active
            return params->effect1.type != (float) FILTER_OFF;
        case FILTER2_PARAM1:
        case FILTER2_PARAM2:
        case FILTER2_AMP:
            // Only route to effect2 params when a modulation effect is active
            return params->effect2.type != (float) FILTER2_OFF;
        default:
            return true;
        }
    };

    auto pickRelevantDestination = [&](const DestinationEnum* pool, int size, const OneSynthParams* params) {
        int start = getRandomInt(size);
        for (int i = 0; i < size; i++) {
            DestinationEnum candidate = pool[(start + i) % size];
            if (isDestinationRelevant(candidate, params)) {
                return candidate;
            }
        }
        return INDEX_ALL_MODULATION;
    };

    static const DestinationEnum safeDestinations[] = {
            INDEX_ALL_MODULATION, INDEX_MODULATION1, INDEX_MODULATION2, INDEX_MODULATION3,
            INDEX_MODULATION4, PAN_OSC1, PAN_OSC2, ALL_PAN,
            FILTER1_PARAM1, FILTER1_PARAM2, FILTER2_PARAM1, FILTER2_PARAM2,
            FILTER1_AMP, FILTER2_AMP, ALL_ENV_ATTACK, ALL_ENV_DECAY,
            ALL_ENV_RELEASE, MTX1_MUL, MTX2_MUL, MTX3_MUL,
            MTX4_MUL, LFO1_FREQ, LFO2_FREQ, LFO3_FREQ
    };

    static const DestinationEnum advancedDestinations[] = {
            OSC1_FREQ, OSC2_FREQ, OSC3_FREQ, OSC4_FREQ,
            OSC5_FREQ, OSC6_FREQ, ALL_OSC_FREQ_HARM,
            OSC1_PHASE, OSC2_PHASE, OSC3_PHASE, OSC4_PHASE,
            OSC1_WARP, OSC2_WARP, OSC3_WARP, OSC4_WARP,
            MTX_DEST_FEEDBACK
    };

    int safeSize = sizeof(safeDestinations) / sizeof(safeDestinations[0]);
    int advSize = sizeof(advancedDestinations) / sizeof(advancedDestinations[0]);

    if (modulationRandom >= 3 && getRandomInt(3) == 0) {
        return pickRelevantDestination(advancedDestinations, advSize, params);
    }
    return pickRelevantDestination(safeDestinations, safeSize, params);
}

/*
 * Randomize the current preset based on the four encoder choices on the
 * MENU_PRESET_RANDOMIZER screen.  Each control is 0-3; 0 (--) skips that
 * section.  Controls also cross-influence each other for coherent results.
 *
 * Oper  0=--  1=soft  2=medi  3=high   oscillator topology
 *   mix:      IM=0: 0.60-1.00  IM=1: 0.53-0.93  IM=2: 0.46-0.86  IM=3: 0.39-0.79
 *   algo:     perc biases toward lower indices; pad toward upper; else full range
 *   shape:    soft=0..2  medi=0..7(no6)  high=0..6
 *   freqMul:  soft={0.5,1,2,4}  medi=7 steps 0.25..4  high=0.25..6.25 (24 steps)
 *   freqType: soft=7/8 kbd  medi=6/8 kbd  high=4/8 kbd  (remainder = fixed pitch)
 *   detune:   soft=0  medi=0..±0.03  high=0..±0.07  (zeroed if outside inner band)
 *   pan:      osc1=0  osc2=±0.3  osc3=±0.3  osc4=±0.5  osc5=±0.5  osc6=±0.7
 *
 * EnvT  0=--  1=perc  2=pad   3=rand   envelope character
 *   perc: attack 0..0.3s   decay 0.05..0.5s   sustain 0..1   release 0.2..5s
 *   pad:  attack 0.5..3s   decay 1..5s        sustain 0..1   release 1..8s
 *   rand: all stages fully random  (times 0..1s; release 0..4s)
 *   Cross-effects: biases algo range, osc shapes, LFO speed/shape, FX pool
 *
 * IM    0=--  1=soft  2=medi  3=high   FM modulation depth
 *   im1-5:    soft=0.25..2   medi=0.5..3   high=1..5
 *   velo1-5:  soft=0.2..1    medi=1..2     high=1..4
 *   feedback: soft=0..0.8    medi=0..1.5   high=0..2.5  (velo up to 80% of max)
 *   Cross-effect: carrier mix and FX gain scale down by 0.07 per IM step
 *
 * Modl  0=--  1=soft  2=medi  3=high   matrix routing
 *   Modl=1:  2 rows  level-1 depth   unique destinations  safe pool only
 *   Modl=2: +3 rows  level-2 depth   unique destinations  safe pool
 *   Modl=3: +6 rows  level-3 depth   reuse allowed        advanced pool (1/3)
 *   Mul ranges  pitch/phase/warp: 0.03..0.45 (0.75 at high)
 *               pan/mix:          0.08..1.2  (1.8 at high)
 *               other:            0.2..1.6   (2.4 medi, 3.5 high); neg allowed at high
 *   Rows 10-12 fixed: modwheel→allIM  pitchbend→allFreq  aftertouch→IM1
 */
void SynthState::randomizePreset() {
    int operatorRandom = fullState.randomizer.Oper;
    int envelopeTypeRandom = fullState.randomizer.EnvT;
    int imRandom = fullState.randomizer.IM;
    int modulationRandom = fullState.randomizer.Modl;

    // general

    params->engine1.velocity = 8;

    // --- Oper: oscillator mix, pan, algorithm, shapes, frequencies, FX ---
    if (operatorRandom > 0) {
        // High IM generates louder FM; scale carrier mix down to avoid clipping
        float mixMin = 0.6f - imRandom * 0.07f;
        float mixMax = 1.0f - imRandom * 0.07f;
        params->engineMix1.mixOsc1 = getRandomFloat(mixMin, mixMax);
        params->engineMix1.mixOsc2 = getRandomFloat(mixMin, mixMax);
        params->engineMix2.mixOsc3 = getRandomFloat(mixMin, mixMax);
        params->engineMix2.mixOsc4 = getRandomFloat(mixMin, mixMax);
        params->engineMix3.mixOsc5 = getRandomFloat(mixMin, mixMax);
        params->engineMix3.mixOsc6 = getRandomFloat(mixMin, mixMax);

        params->engineMix1.panOsc1 = 0.0;
        params->engineMix1.panOsc2 = getRandomFloat(-0.3f, 0.3f);
        params->engineMix2.panOsc3 = getRandomFloat(-0.3f, 0.3f);
        params->engineMix2.panOsc4 = getRandomFloat(-0.5f, 0.5f);
        params->engineMix3.panOsc5 = getRandomFloat(-0.5f, 0.5f);
        params->engineMix3.panOsc6 = getRandomFloat(-0.7f, 0.7f);

        // Bias algorithm toward envelope character:
        // Perc → deep FM chains (lower indices, fewer carriers)
        // Pad  → stacked carriers for lush output (higher indices)
        int algoLow = 0, algoHigh = ALGO_END;
        if (envelopeTypeRandom == 1) {
            algoHigh = ALGO_END * 2 / 3;
        } else if (envelopeTypeRandom == 2) {
            algoLow = ALGO_END / 3;
        }
        params->engine1.algo = algoLow + getRandomInt(algoHigh - algoLow);

        for (int o = 0; o < 6; o++) {
            struct OscillatorParams* currentOsc = &((struct OscillatorParams*) &params->osc1)[o];
            struct OperatorPhaseRowParams* currentPhase = &((struct OperatorPhaseRowParams*) &params->phaseOp1)[o];
            currentOsc->shape = getRandomShape(operatorRandom);
            // Pad mode: pull shapes back toward smooth (sine/tri/saw) 75% of the time
            if (envelopeTypeRandom == 2 && currentOsc->shape > 2 && getRandomInt(4) > 0) {
                currentOsc->shape = getRandomInt(3);
            }
            currentOsc->frequencyMul = getRandomFrequency(operatorRandom);
            currentOsc->frequencyType = getRandomFrequencyType(operatorRandom);
            currentOsc->detune = getFineTune(operatorRandom);
            currentPhase->unused1 = 0.0f;
        }

        // FX - filter and modulation effects tuned to envelope character and IM level
        {
            // Normalize output gain by carrier count.
            // effect1.param3 drives mixerGain (applied to every voice sample), so
            // dividing by numCarriers keeps the summed output at roughly fxBase*avgMix
            // regardless of whether the algorithm has 1 or 6 carrier operators.
            int numCarriers = algoInformation[(int)params->engine1.algo].mix;
            float fxBase = 0.8f - imRandom * 0.07f;
            float fxGain = fxBase / (float)(numCarriers > 0 ? numCarriers : 1);

            // Curated filter pools per envelope type
            // percFilters: HP/BP for transient clarity, distortion for punch,
            //              and open LP variants to warm drums without muffling
            const int percFilters[]    = { FILTER_HP, FILTER_HP2, FILTER_HP3,
                                           FILTER_BASS, FILTER_BP,
                                           FILTER_CRUSHER, FILTER_SAT, FILTER_FOLD,
                                           FILTER_LP, FILTER_LP2, FILTER_LP3 };
            const int padFilters[]     = { FILTER_LP, FILTER_LP2, FILTER_LP3,
                                           FILTER_LPHP, FILTER_LOWSHELF,
                                           FILTER_TILT, FILTER_STEREO, FILTER_BP };

            int chosenType = FILTER_OFF;
            float p1Min = 0.2f, p1Max = 0.8f;
            float p2Min = 0.1f, p2Max = 0.6f;

            if (envelopeTypeRandom == 1) {                   // perc: HP/distortion for punch + open LP
                if (getRandomInt(10) < 6) {
                    chosenType = percFilters[getRandomInt(11)];
                    // LP variants get a high cutoff (0.55+) so drums stay open
                    p1Min = (chosenType == FILTER_LP || chosenType == FILTER_LP2 || chosenType == FILTER_LP3)
                            ? 0.55f : 0.3f;
                    p1Max = 0.9f;
                    p2Min = 0.1f; p2Max = 0.6f;
                }
            } else if (envelopeTypeRandom == 2) {            // pad: LP/shelf to smooth FM partials
                if (getRandomInt(10) < 8) {
                    chosenType = padFilters[getRandomInt(8)];
                    p1Min = 0.25f; p1Max = 0.65f;
                    p2Min = 0.0f;  p2Max = 0.35f;
                }
            } else {                                          // rand: any available filter
                if (getRandomInt(10) < 5) {
                    // Full range: skip FILTER_OFF(0) and FILTER_MIXER(1)
                    chosenType = FILTER_MIXER + 1 + getRandomInt(FILTER_LAST - FILTER_MIXER - 1);
                }
            }

            params->effect1.type   = (float) chosenType;
            params->effect1.param1 = getRandomFloat(p1Min, p1Max);
            params->effect1.param2 = getRandomFloat(p2Min, p2Max);
            params->effect1.param3 = (chosenType == FILTER_OFF) ? 0.9f : fxGain;

            // effect2: modulation effects (chorus/flange/stereo) scaled to EnvT and Modl
            int chosenType2 = FILTER2_OFF;
            if (envelopeTypeRandom == 2) {                   // pad: chorus/ensemble/widener
                const int padFx2[] = { FILTER2_CHORUS, FILTER2_DIMENSION,
                                       FILTER2_WIDE, FILTER2_DIFFUSER, FILTER2_DOUBLER };
                if (getRandomInt(10) < (modulationRandom > 0 ? 9 : 6)) {
                    chosenType2 = padFx2[getRandomInt(5)];
                }
            } else if (envelopeTypeRandom == 1) {            // perc: occasional flange/grain
                if (getRandomInt(10) < 3) {
                    const int percFx2[] = { FILTER2_FLANGE, FILTER2_GRAIN1 };
                    chosenType2 = percFx2[getRandomInt(2)];
                }
            } else {                                          // rand: any modulation effect
                if (getRandomInt(10) < (2 + modulationRandom * 2)) {
                    // Full range: skip FILTER2_OFF(0)
                    chosenType2 = FILTER2_FLANGE + getRandomInt(FILTER2_LAST - FILTER2_FLANGE);
                }
            }
            params->effect2.type   = (float) chosenType2;
            params->effect2.param1 = getRandomFloat(0.15f, 0.55f);
            params->effect2.param2 = getRandomFloat(0.3f,  0.7f);
            params->effect2.param3 = fxGain;
        }
    }

    // --- EnvT: envelope shape applied uniformly to all 6 operators ----------
    for (int e = 0; e < 6; e++) {
        struct EnvelopeParamsA* enva = &((struct EnvelopeParamsA*) &params->env1Time)[e * 2];
        struct EnvelopeParamsB* envb = &((struct EnvelopeParamsB*) &params->env1Level)[e * 2];

        switch (envelopeTypeRandom) {
        case 1:
            enva->attackLevel = 1.0;
            enva->attackTime = getRandomFloat(0, 0.3f);
            enva->decayLevel = getRandomFloat(0.5f, 1.0f);
            enva->decayTime = getRandomFloat(0.05, 0.5f);

            envb->sustainLevel = getRandomFloat(0.0f, 1.0f);
            envb->sustainTime = getRandomFloat(0.02, 1.0f);
            envb->releaseLevel = 0.0f;
            envb->releaseTime = getRandomFloat(0.2, 5.0f);
            ;

            break;
        case 2:
            enva->attackLevel = getRandomFloat(0.25f, 1.0f);
            enva->attackTime = getRandomFloat(0.5f, 3.0f);
            enva->decayLevel = getRandomFloat(0.5f, 1.0f);
            enva->decayTime = getRandomFloat(0.5f, 3.0f);

            envb->sustainLevel = getRandomFloat(0.0f, 1.0f);
            envb->sustainTime = getRandomFloat(0.5f, 2.0f);
            envb->releaseLevel = 0.0f;
            envb->releaseTime = getRandomFloat(1.0f, 5.0f);
            break;
        case 3:
            enva->attackLevel = getRandomFloat(0, 1.0f);
            enva->attackTime = getRandomFloat(0, 1.0f);
            enva->decayLevel = getRandomFloat(0, 1.0f);
            enva->decayTime = getRandomFloat(0, 1.0f);

            envb->sustainLevel = getRandomFloat(0, 1.0f);
            envb->sustainTime = getRandomFloat(0, 1.0f);
            // Always decay fully to zero: a non-zero releaseLevel makes the voice die
            // at non-zero amplitude (click) and can trigger the modulator loop trick
            // (releaseLevel==1 && releaseTime==0 → envelope loops forever).
            envb->releaseLevel = 0.0f;
            envb->releaseTime = getRandomFloat(0, 4.0f);
            break;
        }
    }

    // --- IM: FM modulation indices (im1-im5) and feedback (im6) ------------
    if (imRandom > 0) {
        struct EngineIm1* im1 = (struct EngineIm1*) &params->engineIm1;
        struct EngineIm2* im2 = (struct EngineIm2*) &params->engineIm2;
        struct EngineIm3* im3 = (struct EngineIm3*) &params->engineIm3;

        float min = 0;
        float max = 0;
        float minVelo = 0;
        float maxVelo = 0;

        switch (imRandom) {
        case 1:
            min = 0.25f;
            max = 2.0f;
            minVelo = 0.2f;
            maxVelo = 1.0f;
            break;
        case 2:
            min = .5f;
            max = 3.0f;
            minVelo = 1.0f;
            maxVelo = 2.0f;
            break;
        case 3:
            min = 1.0f;
            max = 5.0f;
            minVelo = 1.0f;
            maxVelo = 4.0f;
            break;
        }
        im1->modulationIndex1 = getRandomFloat(min, max);
        im1->modulationIndexVelo1 = getRandomFloat(minVelo, maxVelo);
        im1->modulationIndex2 = getRandomFloat(min, max);
        im1->modulationIndexVelo2 = getRandomFloat(minVelo, maxVelo);
        im2->modulationIndex3 = getRandomFloat(min, max);
        im2->modulationIndexVelo3 = getRandomFloat(minVelo, maxVelo);
        im2->modulationIndex4 = getRandomFloat(min, max);
        im2->modulationIndexVelo4 = getRandomFloat(minVelo, maxVelo);
        im3->modulationIndex5 = getRandomFloat(min, max);
        im3->modulationIndexVelo5 = getRandomFloat(minVelo, maxVelo);

        // feedback / self-mod path used by many algorithms
        float maxFeedback = 0.6f;
        switch (imRandom) {
        case 1:
            maxFeedback = 0.8f;
            break;
        case 2:
            maxFeedback = 1.5f;
            break;
        case 3:
            maxFeedback = 2.5f;
            break;
        }
        im3->modulationIndex6 = getRandomFloat(0.0f, maxFeedback);
        im3->modulationIndexVelo6 = getRandomFloat(0.0f, maxFeedback * 0.8f);
    }

    // --- Modl: matrix routing, LFOs, LFO envelopes, step sequencers --------
    if (modulationRandom > 0) {
        bool percussiveEnvelope = envelopeTypeRandom == 1;

        params->matrixRowState1.source = MATRIX_SOURCE_LFO1;
        params->matrixRowState2.source = MATRIX_SOURCE_LFO1;
        params->matrixRowState3.source = MATRIX_SOURCE_LFO2;
        params->matrixRowState4.source = MATRIX_SOURCE_LFO2;
        params->matrixRowState5.source = MATRIX_SOURCE_LFO3;
        params->matrixRowState6.source = MATRIX_SOURCE_LFOENV1;
        params->matrixRowState7.source = MATRIX_SOURCE_LFOENV2;
        params->matrixRowState8.source = MATRIX_SOURCE_LFOSEQ1;
        params->matrixRowState9.source = MATRIX_SOURCE_LFOSEQ2;

        params->matrixRowState10.source = MATRIX_SOURCE_MODWHEEL;
        params->matrixRowState10.mul = 2.0f;
        params->matrixRowState10.dest1 = INDEX_ALL_MODULATION;

        params->matrixRowState11.source = MATRIX_SOURCE_PITCHBEND;
        params->matrixRowState11.mul = 1.0f;
        params->matrixRowState11.dest1 = ALL_OSC_FREQ;

        params->matrixRowState12.source = MATRIX_SOURCE_AFTERTOUCH;
        params->matrixRowState12.mul = 1.0f;
        params->matrixRowState12.dest1 = INDEX_MODULATION1;

        for (int m = 1; m <= 9; m++) {
            struct MatrixRowParams* matrixRow = &((struct MatrixRowParams*) &params->matrixRowState1)[m - 1];
            matrixRow->mul = 0;
            matrixRow->dest1 = 0;
        }

        bool rowUsed[9] = { false, false, false, false, false, false, false, false, false };
        bool destinationUsed[DESTINATION_MAX] = { false };

        auto destinationIndex = [](DestinationEnum dest) {
            int idx = (int) dest;
            if (idx < 0 || idx >= DESTINATION_MAX) {
                return 0;
            }
            return idx;
        };

        auto chooseDestination = [&](int level, bool preferUnique) {
            DestinationEnum chosen = getRandomModDestination(level, params);
            if (!preferUnique) {
                destinationUsed[destinationIndex(chosen)] = true;
                return chosen;
            }
            for (int attempt = 0; attempt < 6; attempt++) {
                DestinationEnum candidate = getRandomModDestination(level, params);
                if (!destinationUsed[destinationIndex(candidate)]) {
                    chosen = candidate;
                    break;
                }
            }
            destinationUsed[destinationIndex(chosen)] = true;
            return chosen;
        };

        auto getMatrixMulForDestination = [&](DestinationEnum dest, int level) {
            bool pitchOrPhase = dest == OSC1_FREQ || dest == OSC2_FREQ || dest == OSC3_FREQ || dest == OSC4_FREQ
                    || dest == OSC5_FREQ || dest == OSC6_FREQ || dest == ALL_OSC_FREQ || dest == ALL_OSC_FREQ_HARM
                    || dest == OSC1_PHASE || dest == OSC2_PHASE || dest == OSC3_PHASE || dest == OSC4_PHASE
                    || dest == OSC5_PHASE || dest == OSC6_PHASE || dest == OSC1_WARP || dest == OSC2_WARP
                    || dest == OSC3_WARP || dest == OSC4_WARP || dest == OSC5_WARP || dest == OSC6_WARP;

            bool mixOrPan = dest == PAN_OSC1 || dest == PAN_OSC2 || dest == PAN_OSC3 || dest == PAN_OSC4
                    || dest == ALL_PAN || dest == MIX_OSC1 || dest == MIX_OSC2 || dest == MIX_OSC3
                    || dest == MIX_OSC4 || dest == ALL_MIX;

            float minAbs = 0.2f;
            float maxAbs = 1.6f;
            if (level >= 2) {
                maxAbs = 2.4f;
            }
            if (level >= 3) {
                maxAbs = 3.5f;
            }

            if (pitchOrPhase) {
                minAbs = 0.03f;
                maxAbs = (level >= 3) ? 0.75f : 0.45f;
            } else if (mixOrPan) {
                minAbs = 0.08f;
                maxAbs = (level >= 3) ? 1.8f : 1.2f;
            }

            float mul = getRandomFloat(minAbs, maxAbs);
            if (level >= 3 && !percussiveEnvelope && getRandomInt(4) == 0) {
                mul = -mul;
            }
            return mul;
        };

        auto getRandomMatrixRow = [&](bool avoidReuse) {
            if (!avoidReuse) {
                return getRandomInt(9);
            }

            int available = 0;
            for (int i = 0; i < 9; i++) {
                if (!rowUsed[i]) {
                    available++;
                }
            }
            if (available == 0) {
                return getRandomInt(9);
            }

            int nth = getRandomInt(available);
            for (int i = 0; i < 9; i++) {
                if (!rowUsed[i]) {
                    if (nth == 0) {
                        rowUsed[i] = true;
                        return i;
                    }
                    nth--;
                }
            }
            return getRandomInt(9);
        };

        // Modl=1: 2 base rows; level-1 depth; unique destination preference
        for (int i = 0; i < 2; i++) {
            struct MatrixRowParams* matrixRow = &((struct MatrixRowParams*) &params->matrixRowState1)[getRandomMatrixRow(true)];
            matrixRow->dest1 = chooseDestination(1, true);
            matrixRow->mul = getMatrixMulForDestination((DestinationEnum) matrixRow->dest1, 1);
        }

        // Modl>=2: +3 rows; Modl-level depth; unique destinations; safe pool
        if (modulationRandom >= 2) {
            for (int i = 0; i < 3; i++) {
                struct MatrixRowParams* matrixRow = &((struct MatrixRowParams*) &params->matrixRowState1)[getRandomMatrixRow(true)];
                matrixRow->dest1 = chooseDestination(modulationRandom, true);
                matrixRow->mul = getMatrixMulForDestination((DestinationEnum) matrixRow->dest1, modulationRandom);
            }
        }

        // Modl=3: +6 rows; reuse allowed; advanced destinations unlocked (1/3 chance)
        if (modulationRandom >= 3) {
            for (int i = 0; i < 6; i++) {
                struct MatrixRowParams* matrixRow = &((struct MatrixRowParams*) &params->matrixRowState1)[getRandomMatrixRow(false)];
                matrixRow->dest1 = chooseDestination(modulationRandom, false);
                matrixRow->mul = getMatrixMulForDestination((DestinationEnum) matrixRow->dest1, modulationRandom);
            }
        }

        if (percussiveEnvelope) {
            // Force a couple of transient-oriented LFO routes in percussive mode.
            params->matrixRowState1.mul = getRandomFloat(0.7f, 2.0f);
            params->matrixRowState1.dest1 = INDEX_ALL_MODULATION;
            params->matrixRowState3.mul = getRandomFloat(0.3f, 1.2f);
            params->matrixRowState3.dest1 = FILTER1_PARAM1;
        }

        // LFOs 1-3: shape, frequency, sync mode and phase vary by EnvT
        for (int o = 0; o < 3; o++) {
            struct LfoParams* osc = &((struct LfoParams*) &params->lfoOsc1)[o];
            float* syncModes = &params->lfoSyncModes.lfo1;
            float* lfoPhases = &params->lfoPhases.phaseLfo1;
            if (percussiveEnvelope) {
                // One-shot shapes for transient automation; 1-8 trigger cycles
                // Speed: 1.5..(7+Modl×1.5) Hz; phase 0..0.12 staggers triggers
                // All shapes here have oneShotTerminalShapeValue = -1 (freeze at minimum)
                // so they taper cleanly to zero when the shot completes (with bias=0).
                // RISE_EXP and RISE_LOG are excluded: they freeze at +1, leaving a
                // permanent positive DC offset on the destination after the shot.
                const int percussiveShapes[] = {
                    LFO_DECAY_EXP,         // exponential decay
                    LFO_DECAY_LOG,         // logarithmic decay
                    LFO_DECAY_S,           // S-curve (sigmoid) decay
                    LFO_ATTACK_DECAY,      // attack then decay
                    LFO_ATTACK_HOLD_DECAY, // attack, hold, then decay
                    LFO_BUCHLA_PLONG,      // Buchla-style long pluck
                    LFO_BUCHLA_PLONG2      // Buchla-style long pluck variant
                };
                osc->shape = percussiveShapes[getRandomInt(7)];
                osc->freq = getRandomFloat(1.5f, 7.0f + modulationRandom * 1.5f);
                // Bias must be 0: one-shot shapes taper to 0 then freeze;
                // a non-zero bias leaves a permanent DC offset on the destination.
                osc->bias = 0.0f;
                osc->keybRamp = getRandomFloat(0.0f, 0.8f);
                syncModes[o] = (float) (LFO_SYNC_ONESHOT_INTERNAL_1 + getRandomInt(8));
                lfoPhases[o] = getRandomFloat(0.0f, 0.12f);
            } else {
                if (envelopeTypeRandom == 2) {  // pad: slow smooth LFOs for subtle movement
                    // Smooth, zero-mean shapes only. Excluded: SIN_POS/SIN_ZERO/SIN_SQUARE
                    // (always-positive output = DC offset on destination).
                    const int padShapes[] = {
                        LFO_SIN,       // sine
                        LFO_TRIANGLE,  // triangle
                        LFO_SAW,       // ramp up
                        LFO_SAW_DOWN,  // ramp down
                        LFO_BROWNIAN,  // organic random walk
                        LFO_WANDERING, // smooth wandering random
                        LFO_FLOW       // smooth flowing random
                    };
                    osc->shape = padShapes[getRandomInt(7)];
                    osc->freq = getRandomFloat(0.03f, 1.0f + modulationRandom * 0.5f);
                    osc->bias = 0;
                    osc->keybRamp = getRandomFloat(0.0f, 4.0f);
                } else {
                    // --/rand: wide variety including stochastic shapes
                    const int randShapes[] = {
                        LFO_SIN,       // sine
                        LFO_SAW,       // ramp up
                        LFO_TRIANGLE,  // triangle
                        LFO_SQUARE,    // square
                        LFO_RANDOM,    // sample-and-hold random
                        LFO_BROWNIAN,  // random walk
                        LFO_WANDERING, // smooth wandering
                        LFO_SAW_DOWN,  // ramp down
                        LFO_FLOW       // smooth flowing random
                    };
                    osc->shape = randShapes[getRandomInt(9)];
                    osc->freq = getRandomFloat(0.2f, 3.0f + modulationRandom * 2.0f);
                    // Bias 0: a non-zero bias creates a DC offset on the destination.
                    // When routed to ALL_ENV_RELEASE or an amplitude destination,
                    // a positive bias can make notes sound much longer than intended.
                    osc->bias = 0.0f;
                    osc->keybRamp = getRandomFloat(0.0f, 1.0f);
                }
                syncModes[o] = LFO_SYNC_INTERNAL;
                lfoPhases[o] = 0.0f;
            }
        }
        params->lfoSyncModes.unused1 = 0.0f;
        // LFO envelopes 1-2: all four stages 0.05..1s; lfoEnv2 loop count 1 or 2
        for (int e = 0; e < 2; e++) {
            struct EnvelopeLfoParams* env = &((struct EnvelopeLfoParams*) &params->lfoEnv1)[e];
            env->attack = getRandomFloat(0.05f, 1.0f);
            env->decay = getRandomFloat(0.05f, 1.0f);
            env->sustain = getRandomFloat(0.05f, 1.0f);
            if (e == 0) {
                env->release = getRandomFloat(0.05f, 1.0f);
            } else {
                params->lfoEnv2.loop = getRandomInt(2) + 1;
            }
        }

        // Step sequencers 1-2: shared BPM 60-180; gate 0.25-1; accents (11-15) on every 4th step
        int bpm = getRandomInt(120) + 60;
        for (int s = 0; s < 2; s++) {
            struct StepSequencerParams* stepSeq = &((struct StepSequencerParams*) &params->lfoSeq1)[s];
            stepSeq->bpm = bpm;
            stepSeq->gate = getRandomFloat(0.25, 1);

            struct StepSequencerSteps* steps = &((struct StepSequencerSteps*) &params->lfoSteps1)[s];
            for (int k = 0; k < 16; k++) {
                if ((k % 4) == 0) {
                    steps->steps[k] = getRandomInt(5) + 11;
                } else {
                    steps->steps[k] = getRandomInt(10);
                }
            }
        }
    }
}

char* SynthState::getTimbreName(int t) {
    return timbres[t].getPresetName();
}

uint8_t SynthState::getTimbrePlayMode(int t) {
    return timbres[t].getParamRaw()->engine1.playMode;
}


/*
 * When a scala scale file does not exist anymore ?
 * If both number and name are wrong, scala is disable
 */
bool SynthState::scalaSettingsChanged(int timbre) {
    if (!mixerState.instrumentState_[timbre].scalaEnable) {
        mixerState.instrumentState_[timbre].scaleFrequencies = diatonicScaleFrequency;
    } else {
        const struct PFM3File*  scalaFile = storage->getScalaFile()->getFile(mixerState.instrumentState_[timbre].scaleScaleNumber);
        if (scalaFile->fileType == FILE_EMPTY) {
            // Lest file
            return false;
        }
        // Let's fill scala Name (used in loadScalaScale)
        for (int c = 0; c < 12; c++) {
            mixerState.instrumentState_[timbre].scalaScaleFileName[c] = scalaFile->name[c];
        }
        float* newScaleFrequencies = storage->getScalaFile()->loadScalaScale(&mixerState, timbre);
        if (newScaleFrequencies != 0) {
            mixerState.instrumentState_[timbre].scaleFrequencies = newScaleFrequencies;
        } else {
            return false;
        }
    }
    return true;
}

const char* SynthState::getSequenceName() {
    return displaySequencer->getSequenceName();
}


