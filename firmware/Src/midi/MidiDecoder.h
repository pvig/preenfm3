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

#ifndef MIDIDECODER_H_
#define MIDIDECODER_H_

#include "SynthStateAware.h"
#include "Synth.h"
#include "RingBuffer.h"
#include "VisualInfo.h"
#include "Storage.h"

// number of external control change
#define NUMBER_OF_ECC 4

struct MidiEventState {
    EventState eventState;
    uint8_t numberOfBytes;
    uint16_t index;
};

struct MidiEvent {
    unsigned char channel;
    EventType eventType;
    unsigned char value[2];
};
    
enum AllControlChange {
    // ------------------------------------------------------------------------
    // Standard MIDI CCs & Channel / Performance Controls (0 - 15)
    // ------------------------------------------------------------------------
    CC_BANK_SELECT          = 0,
    CC_MODWHEEL             = 1,
    CC_BREATH               = 2,
    CC_GLIDE                = 5,
    CC_MIXER_VOLUME         = 7,
    CC_MIXER_PAN            = 10,
    CC_MIXER_SEND           = 11,
    CC_UNISON_DETUNE        = 13,
    CC_UNISON_SPREAD        = 14,

    // ------------------------------------------------------------------------
    // Engine Algorithm & FM Modulation Indexes (16 - 21, 30)
    // ------------------------------------------------------------------------
    CC_ALGO                 = 16,
    CC_IM1                  = 17,
    CC_IM2                  = 18,
    CC_IM3                  = 19,
    CC_IM4                  = 20,
    CC_IM5                  = 21,

    // ------------------------------------------------------------------------
    // Operator Mix and Pan (22 - 29, 31, 36, 37, 39)
    // Interleaved pattern: Even = Mix, Odd = Pan (-1.0 to +1.0, 64 = center)
    // ------------------------------------------------------------------------
    CC_MIX1                 = 22,
    CC_PAN1                 = 23,
    CC_MIX2                 = 24,
    CC_PAN2                 = 25,
    CC_MIX3                 = 26,
    CC_PAN3                 = 27,
    CC_MIX4                 = 28,
    CC_PAN4                 = 29,

    // Feedback modulation index (0.0 to 1.0)
    CC_IM_FEEDBACK          = 30,

    CC_MIX5                 = 31,

    // ------------------------------------------------------------------------
    // Bank Select LSB, Operator Warp (32 - 35) & Operator 5-6 Mix/Pan (36, 37, 39)
    // Note: Warp OP1-3 were moved from 97-99 to avoid NRPN protocol collisions.
    // ------------------------------------------------------------------------
    CC_BANK_SELECT_LSB      = 32,
    CC_OSC1_WARP            = 33,
    CC_OSC2_WARP            = 34,
    CC_OSC3_WARP            = 35,
    CC_PAN5                 = 36,
    CC_MIX6                 = 37,
    // Note: CC 38 is standard MIDI Data Entry LSB for NRPN protocol.
    CC_PAN6                 = 39,

    // ------------------------------------------------------------------------
    // Timbre Settings & Master FX (40 - 45)
    // ------------------------------------------------------------------------
    CC_TIMBRE_DRYWET        = 40,
    CC_TIMBRE_DYNAMICS      = 41,
    CC_ENGINE_DECIMATION    = 42,
    CC_PLAY_MODE            = 43,
    CC_MFX_RETURN           = 44,
    CC_MFX_PRESET           = 45,

    // ------------------------------------------------------------------------
    // Modulation Matrix Multipliers (46 - 49)
    // ------------------------------------------------------------------------
    CC_MATRIXROW1_MUL       = 46,
    CC_MATRIXROW2_MUL       = 47,
    CC_MATRIXROW3_MUL       = 48,
    CC_MATRIXROW4_MUL       = 49,

    // ------------------------------------------------------------------------
    // Operator Frequencies (50 - 55)
    // ------------------------------------------------------------------------
    CC_OSC1_FREQ            = 50,
    CC_OSC2_FREQ            = 51,
    CC_OSC3_FREQ            = 52,
    CC_OSC4_FREQ            = 53,
    CC_OSC5_FREQ            = 54,
    CC_OSC6_FREQ            = 55,

    // ------------------------------------------------------------------------
    // LFO Frequencies & Env Silence (56 - 59)
    // ------------------------------------------------------------------------
    CC_LFO1_FREQ            = 56,
    CC_LFO2_FREQ            = 57,
    CC_LFO3_FREQ            = 58,
    CC_LFO_ENV2_SILENCE     = 59,

    // ------------------------------------------------------------------------
    // Envelope Global Shortcuts & Hold Pedal (62 - 64)
    // ------------------------------------------------------------------------
    CC_ENV_ATK_ALL_MODULATOR = 62,
    CC_ENV_REL_ALL_MODULATOR = 63,
    CC_HOLD_PEDAL            = 64,

    // ------------------------------------------------------------------------
    // Operator 1 Attack & Filters 1 & 2 (65 - 74)
    // Note: OP1 Attack is at 65 while OP2-6 are at 75-79 to preserve legacy
    // filter and brightness (CC74) assignments between them.
    // ------------------------------------------------------------------------
    CC_ENV_ATK_OP1           = 65,
    CC_FILTER2_TYPE          = 66,
    CC_FILTER2_PARAM1        = 67,
    CC_FILTER2_PARAM2        = 68,
    CC_FILTER2_MIX           = 69,
    CC_FILTER_TYPE           = 70,
    CC_FILTER_PARAM1         = 71,
    CC_FILTER_PARAM2         = 72,
    CC_FILTER_GAIN           = 73,
    CC_MPE_SLIDE_CC74        = 74,

    // ------------------------------------------------------------------------
    // Operator 2-6 Attack (75 - 79)
    // ------------------------------------------------------------------------
    CC_ENV_ATK_OP2           = 75,
    CC_ENV_ATK_OP3           = 76,
    CC_ENV_ATK_OP4           = 77,
    CC_ENV_ATK_OP5           = 78,
    CC_ENV_ATK_OP6           = 79,

    // ------------------------------------------------------------------------
    // Envelope Carrier Shortcuts & Operator Release (80 - 87)
    // ------------------------------------------------------------------------
    CC_ENV_ATK_ALL_CARRIER   = 80,
    CC_ENV_REL_ALL_CARRIER   = 81,
    CC_ENV_REL_OP1           = 82,
    CC_ENV_REL_OP2           = 83,
    CC_ENV_REL_OP3           = 84,
    CC_ENV_REL_OP4           = 85,
    CC_ENV_REL_OP5           = 86,
    CC_ENV_REL_OP6           = 87,

    // ------------------------------------------------------------------------
    // LFO Phases, Bias & Shapes (88 - 96)
    // ------------------------------------------------------------------------
    CC_LFO1_PHASE            = 88,
    CC_LFO2_PHASE            = 89,
    CC_LFO3_PHASE            = 90,
    CC_LFO1_BIAS             = 91,
    CC_LFO2_BIAS             = 92,
    CC_LFO3_BIAS             = 93,
    CC_LFO1_SHAPE            = 94,
    CC_LFO2_SHAPE            = 95,
    CC_LFO3_SHAPE            = 96,

    // Note: CC 96 = NRPN Data Increment, CC 97 = NRPN Data Decrement,
    //       CC 98 = NRPN Parameter LSB, CC 99 = NRPN Parameter MSB.
    // These are reserved by the MIDI standard NRPN specification.

    // ------------------------------------------------------------------------
    // Arpeggiator (100 - 105)
    // ------------------------------------------------------------------------
    CC_ARP_CLOCK             = 100,
    CC_ARP_DIRECTION         = 101,
    CC_ARP_OCTAVE            = 102,
    CC_ARP_PATTERN           = 103,
    CC_ARP_DIVISION          = 104,
    CC_ARP_DURATION          = 105,

    // ------------------------------------------------------------------------
    // Sequencer (106 - 110)
    // ------------------------------------------------------------------------
    CC_SEQ_START_ALL         = 106,
    CC_SEQ_START_INST        = 107,
    CC_SEQ_RECORD_INST       = 108,
    CC_SEQ_SET_SEQUENCE      = 109,
    CC_SEQ_TRANSPOSE         = 110,

    // ------------------------------------------------------------------------
    // Operator 4-6 Warp (111 - 113)
    // ------------------------------------------------------------------------
    CC_OSC4_WARP             = 111,
    CC_OSC5_WARP             = 112,
    CC_OSC6_WARP             = 113,

    // ------------------------------------------------------------------------
    // Performance Matrix Sources & Current Instrument (115 - 119)
    // ------------------------------------------------------------------------
    CC_MATRIX_SOURCE_CC1     = 115,
    CC_MATRIX_SOURCE_CC2     = 116,
    CC_MATRIX_SOURCE_CC3     = 117,
    CC_MATRIX_SOURCE_CC4     = 118,
    CC_CURRENT_INSTRUMENT    = 119,

    // ------------------------------------------------------------------------
    // MIDI Channel Mode Messages (120 - 127)
    // ------------------------------------------------------------------------
    CC_ALL_SOUND_OFF         = 120,
    CC_ALL_NOTES_OFF         = 123,
    CC_OMNI_OFF              = 124,
    CC_OMNI_ON               = 125,
    CC_RESET                 = 127
};

struct Nrpn {
    unsigned char paramLSB;
    unsigned char paramMSB;
    unsigned char valueLSB;
    unsigned char valueMSB;
    bool readyToSend;
};


enum ActionType {
    LOAD_PRESET,
    SEND_PATCH_AS_NRPN
};

struct AsyncActionDetail {
    uint8_t actionType;
    uint8_t timbre;
    uint8_t param1;
    uint8_t param2;
    uint8_t param3;
    uint8_t param4;
    uint8_t param5;
    uint8_t param6;
};

struct AsyncAction {
    union {
        struct AsyncActionDetail action;
        uint64_t fullBytes;
    };
};

class MidiDecoder: public SynthParamListener, public SynthStateAware {
public:
    MidiDecoder();
    virtual ~MidiDecoder();
    void setStorage(Storage* storage) {
        this->storage = storage;
    }

    void newByte(unsigned char byte);
    void newMessageType(unsigned char byte);
    void newMessageData(unsigned char byte);
    void midiEventReceived(MidiEvent& midiEvent);
    void midiEventForInstrument1MPE(MidiEvent& midiEvent);
    void controlChange(int timbre, MidiEvent& midiEvent, bool allowNrpnDump = true);
    void decodeNrpn(int timbre, int midiChannel = -1, bool allowDump = true);
    void setSynth(Synth* synth);
    void setVisualInfo(VisualInfo* visualInfo);

    void newParamValueFromExternal(int timbre, int currentrow, int encoder, ParameterDisplay* param, float oldValue, float newValue);
    void newMixerValueFromExternal(int timbre, int mixerValueType, float oldValue, float newValue);
    void newParamValue(int timbre, int currentrow, int encoder, ParameterDisplay* param, float oldValue, float newValue);
    void newMixerValue(uint8_t valueType, uint8_t timbre, float oldValue, float newValue);
    void newMixerEdit(int oldButton, int newButton) {
    }

    void beforeNewParamsLoad(int timbre) {
    }
    void afterNewParamsLoad(int timbre) {
    }
    void afterNewMixerLoad() {
    }
    void showAlgo() {
    }
    void showIMInformation() {
    }

    void writeMidiCCOut(struct MidiEvent *toSend);
    void sendMidiDin5Out();
    void sendMidiUsbOut();
    void sendMidiUsbOutIfBufferFull();
    void playNote(int timbre, char note, char velocity) {
    }
    void stopNote(int timbre, char note) {
    }
    void newTimbre(int timbre) {
        currentTimbre = timbre;
    }
    void sendCurrentPatchAsNrpns(int timbre, int midiChannel = -1);

    // Firmware 2.00
    // Phase LFO1/3 added not at the right place so nrpm and params row are now
    // unlinked for compatibility reason....
    int getMidiIndexFromMemory(int index);
    int getMemoryIndexFromMidi(int index);

    // Some actions must not be called from the audio thread
    void processAsyncActions();
private:
    uint8_t analyseSysexBuffer(uint8_t *sysexBuffer, uint16_t size);

    struct MidiEventState currentEventState;
    struct MidiEvent currentEvent;
    Synth* synth;
    VisualInfo *visualInfo;
    Storage* storage;
    int currentTimbre;
    struct MidiEvent toSend;
    struct MidiEvent lastSentCC;
    struct Nrpn currentNrpn[NUMBER_OF_TIMBRES];
    bool omniOn[NUMBER_OF_TIMBRES];
    unsigned char runningStatus;

    // Midi Clock
    bool isExternalMidiClockStarted;
    int midiClockCpt;
    int songPosition;

    // int bank number
    char bankNumber[NUMBER_OF_TIMBRES];
    char bankNumberLSB[NUMBER_OF_TIMBRES];

};

#endif /* MIDIDECODER_H_ */
