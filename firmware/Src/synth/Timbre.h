/*
 * Copyright 2013 Xavier Hosxe
 *
 * Author: Xavier Hosxe (xavier <.> hosxe < a t > gmail.com)
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

#ifndef TIMBRE_H_
#define TIMBRE_H_

#include "Common.h"

#include "Osc.h"
#include "Env.h"
#include "Lfo.h"
#include "LfoOsc.h"
#include "LfoEnv.h"
#include "LfoEnv2.h"
#include "LfoStepSeq.h"
#include "Matrix.h"
#include "note_stack.h"
#include "event_scheduler.h"

extern float panTable[];
class Voice;

enum {
    CLOCK_OFF,
    CLOCK_INTERNAL,
    CLOCK_EXTERNAL
};

class Timbre {
    friend class Synth;
    friend class Voice;
public:
    Timbre();
    virtual ~Timbre();
    void init(SynthState *synthState, int timbreNumber);
    void setVoiceNumber(int v, int n);
    void initVoicePointer(int n, Voice *voice);
    void updateArpegiatorInternalClock();
    void cleanNextBlock();
    void prepareMatrixForNewBlock();
    uint8_t voicesNextBlock();
    void glide();
    void voicesToTimbre(float volumeGain);
    void gateFx();
    void afterNewParamsLoad();
    void setNewValue(int index, struct ParameterDisplay *param, float newValue);
    void setNewEffecParam(int encoder);
    int getSeqStepValue(int whichStepSeq, int step);
    void setSeqStepValue(int whichStepSeq, int step, int value);
    // Arpegiator
    void arpeggiatorNoteOn(char note, char velocity);
    void arpeggiatorNoteOff(char note);
    void StartArpeggio();
    void StepArpeggio();
    void Start();
    void arpeggiatorSetHoldPedal(uint8_t value);
    void setLatchMode(uint8_t value);
    void setDirection(uint8_t value);
    void setNewBPMValue(float bpm);
    void setArpeggiatorClock(float bpm);
    void resetArpeggiator();
    uint16_t getArpeggiatorPattern() const;

    void noteOn(char note, char velocity);
    void noteOff(char note);

    void noteOnMPE(uint8_t channel, uint8_t note, uint8_t velocity);
    void noteOffMPE(uint8_t channel, uint8_t note, uint8_t velocityOff);

    void preenNoteOn(char note, char velocity);
    inline void preenNoteOnUpdateMatrix(int voiceToUse, int note, int velocity);
    void preenNoteOff(char note);

    void numberOfVoicesChanged(uint8_t newNumberOfVoices) {
        if (likely(newNumberOfVoices > 0)) {
            numberOfVoiceInverse_ = 1.0f / (float) newNumberOfVoices;
        } else {
            numberOfVoiceInverse_ = 1.0f;
        }
        numberOfVoices_ = newNumberOfVoices;
    }

    void lfoValueChange(int currentRow, int encoder, float newValue);

    void setHoldPedal(int value);

    void resetMatrixDestination(float oldValue);

    void setMatrixSource(enum SourceEnum source, float newValue);
    void setMatrixSourceMPE(uint8_t channel, enum SourceEnum source, float newValue);

    void setMatrixPolyAfterTouch(uint8_t note, float newValue);
    void verifyLfoUsed(int encoder, float oldValue, float newValue);

    void midiClockStop() {
        OnMidiStop();
    }

    void midiClockContinue(int songPosition);
    void midiClockStart();
    void midiClockSongPositionStep(int songPosition);

    struct OneSynthParams* getParamRaw() {
        return &params_;
    }

    float* getSampleBlock() {
        return sampleBlock_;
    }

    const float* getSampleBlock() const {
        return sampleBlock_;
    }

    int8_t voiceNumber_[MAX_NUMBER_OF_VOICES];

    // Midi note response
    // Midi Note Scale
    void updateMidiNoteScale(int scale);

    // Do matrix use LFO
    bool isLfoUsed(int lfo) {
        return lfoUSed_[lfo] > 0;
    }

    uint8_t getLowerNote() {
        return lowerNote_;
    }
    float getLowerNoteFrequency() {
        return lowerNoteFrequency;
    }

    char* getPresetName() {
        return params_.presetName;
    }

    float getNumberOfVoiceInverse() {
        return numberOfVoiceInverse_;
    }

    bool isUnisonMode() {
    	return params_.engine1.playMode == PLAY_MODE_UNISON;
    }

    void stopPlayingNow();


    void setMPESetting(uint8_t mpeSetting) {
        mpeSetting_ = mpeSetting;
    }

    uint8_t getMPESetting() {
        return mpeSetting_;
    }

private:

    // MiniPal Arpegiator
    void SendLater(uint8_t note, uint8_t velocity, uint8_t when, uint8_t tag);
    void SendScheduledNotes();
    void FlushQueue();
    void Tick();
    void OnMidiContinue();
    void OnMidiStart();
    void OnMidiStop();
    void OnMidiClock();
    void SendNote(uint8_t note, uint8_t velocity);

    /** --------------FX conf--------------  */
    void fxAfterBlock();
    float delayInterpolation(float readPos, float buffer[], int bufferLenM1);
    float delayInterpolation2(float readPos, float buffer[], int bufferLenM1, int offset);
    float iirFilter(float x, float a0, float *yn1, float *yn2, float *xn1, float *xn2) ;

    #define delayBufferSize 2048
    const float delayBufferSizeF       = delayBufferSize;
    const float delayBufferSize120     = delayBufferSize * 0.3333f;
    const float delayBufferSize240     = delayBufferSize * 0.6666f;
    const float delayBufferSize90      = delayBufferSize * 0.25f;
    const float delayBufferSize180     = delayBufferSize * 0.5f;
    const float delayBufferSize270     = delayBufferSize * 0.75f;
    const int delayBufferSizeM1   = delayBufferSize - 1;
    const int delayBufferSizeM4   = delayBufferSize - 4;
    const float delayBufferSizeInv = 1.0f / delayBufferSize;

    static float delayBuffer[NUMBER_OF_TIMBRES][delayBufferSize];
    float *delayBuffer_;

    float param1S = 0;
    float matrixFilterFrequencyS = 0;
    float param2S = 0;
    float delaySize1 = 0, delaySize2 = 0, delaySize3 = 0, delaySize4 = 0;
    float delaySizeInc1 = 0, delaySizeInc2 = 0, delaySizeInc3 = 0, delaySizeInc4 = 0;
    float delayOut1 = 0, delayOut2 = 0, delayOut3 = 0, delayOut4 = 0, delayOut5 = 0, delayOut6 = 0;
    float feedback            = 0;
    float shift = 0, shift2 = 0;
    int delayWritePos         = 0;
    int delayWritePos2        = 0;
    float delayWritePosF      = 0;
    float delayReadPos        = 0;
    float delayReadPos2       = 0;
    float delayReadFrac       = 0;
    float delayReadFrac2      = 0;
    float readPos             = 0;
    float _in_lp_a, _in_lp_b;
    float _in_lp2_a, _in_lp2_b;
    float lpF, lpF2;
     
    float low1 = 0, band1 = 0;
    float low2 = 0, band2 = 0;
    float low3 = 0, band3 = 0;
    float low4 = 0, band4 = 0;
    float low5 = 0, band5 = 0;
    float low6 = 0, band6 = 0;
    float low7 = 0, band7 = 0;
    float low8 = 0, band8 = 0;

    const int delayBufStereoSize = delayBufferSize * 0.5f;
    const float delayBufStereoSizeF = delayBufStereoSize;
    const int delayBufStereoSizeM1 = delayBufStereoSize - 1;
    const float delayBufStereoDiv4 = delayBufStereoSize * 0.25f;
    const float delayBufStereoSizeInv = 1.0f / delayBufStereoSize;

    float delaySumOut = 0;
    float delayIn = 0;

    // hp filter
    float hp_in_x0 = 0;
    float hp_in_y0 = 0;
    float hp_in_y1 = 0;
    float hp_in_x1 = 0;
    float hp_in2_x0 = 0;
    float hp_in2_y0 = 0;
    float hp_in2_y1 = 0;
    float hp_in2_x1 = 0;
    float hp_in3_x0 = 0;
    float hp_in3_y0 = 0;
    float hp_in3_y1 = 0;
    float hp_in3_x1 = 0;
    float _in_b1, _in_a0, _in_a1;
    float _in2_b1, _in2_a0, _in2_a1;
    float _in3_b1, _in3_a0, _in3_a1;

    // allpass filters
    float _ly1 = 0;
    float _lx1 = 0;
    float _ly2 = 0;
    float _lx2 = 0;
    float _ly3 = 0;
    float _lx3 = 0;
    float _ly4 = 0;
    float _lx4 = 0;
    float apcoef1, apcoef2, apcoef3, apcoef4;

    // frequency shifter 
    float hb1_x1 = 0, hb1_x2 = 0, hb1_y1 = 0, hb1_y2 = 0;
    float hb2_x1 = 0, hb2_x2 = 0, hb2_y1 = 0, hb2_y2 = 0;
    float hb3_x1 = 0, hb3_x2 = 0, hb3_y1 = 0, hb3_y2 = 0;
    float hb4_x1 = 0, hb4_x2 = 0, hb4_y1 = 0, hb4_y2 = 0;
    float hb5_x1 = 0, hb5_x2 = 0, hb5_y1 = 0, hb5_y2 = 0;
    float hb6_x1 = 0, hb6_x2 = 0, hb6_y1 = 0, hb6_y2 = 0;
    float hb7_x1 = 0, hb7_x2 = 0, hb7_y1 = 0, hb7_y2 = 0;
    float hb8_x1 = 0, hb8_x2 = 0, hb8_y1 = 0, hb8_y2 = 0;
    float phase1 = 0;
    float samplen1 = 0;
    float shifterOutMix = 0;

    // diffuser
    int inputWritePos1     = 0;
    int inputWritePos2     = 0;
    int inputWritePos3     = 0;
    int inputWritePos4     = 0;
    int inputWritePos5     = 0;
    const int inputBufferLen1 = 112;
    const int inputBufferLen2 = 210;
    const int inputBufferLen3 = 137;
    const int inputBufferLen4 = 242;
    const int inputBufferLen5 = 160;

    /** --------------end of FX conf--------------  */

    int8_t timbreNumber_;
    struct OneSynthParams params_;
    struct MixerState *mixerState_;
    float sampleBlock_[BLOCK_SIZE * 2];
    float *sbMax_;
    float numberOfVoiceInverse_;
    // numberOfVoices is to be used instead of params.engine1.numberOfVoice
    // Because numberOfVoices can be just incremented, and we're not sure the voice is ready.
    // numberOfVoices is incremented after the new voices are initialized and ready to use.
    float numberOfVoices_;

    float mixerGain_;
    Voice *voices_[MAX_NUMBER_OF_VOICES];
    bool holdPedal_;
    int8_t lastPlayedNote_;

    // 6 oscillators Max
    Osc osc1_;
    Osc osc2_;
    Osc osc3_;
    Osc osc4_;
    Osc osc5_;
    Osc osc6_;

    // And their 6 envelopes
    Env env1_;
    Env env2_;
    Env env3_;
    Env env4_;
    Env env5_;
    Env env6_;

    // Must recompute LFO steps ?
    bool recomputeNext_;
    float currentGate_;
    // Arpeggiator

    // TO REFACTOR
    float ticksPerSecond_;
    float ticksEveryNCalls_;
    int ticksEveyNCallsInteger_;

    float arpegiatorStep_;
    NoteStack note_stack_;
    EventScheduler event_scheduler_;

    uint8_t running_;
    uint8_t latch_;
    uint8_t tick_;
    uint8_t idle_ticks_;
    uint16_t bitmask_;
    int8_t current_direction_;
    int8_t current_octave_;
    int8_t current_step_;
    int8_t start_step_;
    uint8_t ignore_note_off_messages_;
    uint8_t recording_;

    // lfoUsed
    uint8_t lfoUSed_[NUMBER_OF_LFO];
    uint8_t lowerNote_;
    float lowerNoteFrequency;
    bool lowerNoteReleased_;
    // static
    static uint32_t voiceIndex_;

    // Unison phase
    static float unisonPhase[14];

    uint8_t mpeSetting_;

};

#endif /* TIMBRE_H_ */
