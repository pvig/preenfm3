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


#include <math.h>
#include "Timbre.h"
#include "Voice.h"

#define INV127 .00787401574803149606f
#define INV16 .0625

// Regular memory
float midiNoteScale[2][NUMBER_OF_TIMBRES][128] __attribute__((section(".ram_d1")));
float Timbre::unisonPhase[14] = { .37f, .11f, .495f, .53f, .03f, .19f, .89f, 0.23f, .71f, .19f, .31f, .43f, .59f, .97f };
float Timbre::delayBuffer[NUMBER_OF_TIMBRES][delayBufferSize] __attribute__ ((section(".ram_d2b")));

struct CounterpointTimbreState {
    float   ticksEveryNCalls;
    float   arpegiatorStep;
    uint8_t activeNote;
    uint8_t prevNote;
    uint8_t clockMode;
    uint8_t running;
    uint8_t tickedThisBlock;
    uint8_t stepTriggeredThisTick;
    uint8_t curPatternStep;
};

static CounterpointTimbreState cptTimbreBus[NUMBER_OF_TIMBRES] __attribute__((section(".ram_d1")));

#define CALLED_PER_SECOND (PREENFM_FREQUENCY / 32.0f)

enum ArpeggiatorDirection {
    ARPEGGIO_DIRECTION_UP = 0,
    ARPEGGIO_DIRECTION_DOWN,
    ARPEGGIO_DIRECTION_UP_DOWN,
    ARPEGGIO_DIRECTION_PLAYED,
    ARPEGGIO_DIRECTION_RANDOM,
    ARPEGGIO_DIRECTION_CHORD,
   /*
    * ROTATE modes rotate the first note played, e.g. UP: C-E-G -> E-G-C -> G-C-E -> repeat
    */
    ARPEGGIO_DIRECTION_ROTATE_UP, ARPEGGIO_DIRECTION_ROTATE_DOWN, ARPEGGIO_DIRECTION_ROTATE_UP_DOWN,
   /*
    * SHIFT modes rotate and extend with transpose, e.g. UP: C-E-G -> E-G-C1 -> G-C1-E1 -> repeat
    */
    ARPEGGIO_DIRECTION_SHIFT_UP, ARPEGGIO_DIRECTION_SHIFT_DOWN, ARPEGGIO_DIRECTION_SHIFT_UP_DOWN,
    ARPEGGIO_DIRECTION_DOWN_UP,
    ARPEGGIO_DIRECTION_UP_DOWN_INC,
    ARPEGGIO_DIRECTION_CONVERGE,
    ARPEGGIO_DIRECTION_DIVERGE,
    ARPEGGIO_DIRECTION_THUMB,
    ARPEGGIO_DIRECTION_PINKY,
    ARPEGGIO_DIRECTION_RANDOM_WALK,
    ARPEGGIO_DIRECTION_STRIDE,
    ARPEGGIO_DIRECTION_LEAP,
    ARPEGGIO_DIRECTION_OCTAVE_NOTE,

    ARPEGGIO_DIRECTION_COUNT
};

// TODO Maybe add something like struct ArpDirectionParams { dir, can_change, use_start_step }

inline static int __getDirection( int _direction ) {
    switch( _direction ) {
    case ARPEGGIO_DIRECTION_DOWN:
    case ARPEGGIO_DIRECTION_ROTATE_DOWN:
    case ARPEGGIO_DIRECTION_SHIFT_DOWN:
    case ARPEGGIO_DIRECTION_DOWN_UP:
        return -1;
    default:
        return 1;
    }
}

inline static int __canChangeDir( int _direction ) {
    switch( _direction ) {
    case ARPEGGIO_DIRECTION_UP_DOWN:
    case ARPEGGIO_DIRECTION_ROTATE_UP_DOWN:
    case ARPEGGIO_DIRECTION_SHIFT_UP_DOWN:
    case ARPEGGIO_DIRECTION_DOWN_UP:
    case ARPEGGIO_DIRECTION_UP_DOWN_INC:
        return 1;
    default:
        return 0;
    }
}

inline static int __canTranspose( int _direction ) {
    switch( _direction ) {
    case ARPEGGIO_DIRECTION_SHIFT_UP:
    case ARPEGGIO_DIRECTION_SHIFT_DOWN:
    case ARPEGGIO_DIRECTION_SHIFT_UP_DOWN:
        return 1;
    default:
        return 0;
    }
}

arp_pattern_t lut_res_arpeggiator_patterns[ ARPEGGIATOR_PRESET_PATTERN_COUNT ]  = {
  ARP_PATTERN_ACCENT(21845, 0x1041), ARP_PATTERN_ACCENT(62965, 0x50a1), ARP_PATTERN_ACCENT(46517, 0x1081), ARP_PATTERN_ACCENT(54741, 0x4141),
  ARP_PATTERN_ACCENT(43861, 0x0841), ARP_PATTERN_ACCENT(22869, 0x1041), ARP_PATTERN_ACCENT(38293, 0x1081), ARP_PATTERN_ACCENT(2313,  0x0101),
  ARP_PATTERN_ACCENT(37449, 0x0201), ARP_PATTERN_ACCENT(21065, 0x0201), ARP_PATTERN_ACCENT(18761, 0x0101), ARP_PATTERN_ACCENT(54553, 0x4101),
  ARP_PATTERN_ACCENT(27499, 0x0841), ARP_PATTERN_ACCENT(23387, 0x1111), ARP_PATTERN_ACCENT(30583, 0x1111), ARP_PATTERN_ACCENT(28087, 0x2111),
  ARP_PATTERN_ACCENT(22359, 0x0441), ARP_PATTERN_ACCENT(28527, 0x2449), ARP_PATTERN_ACCENT(30431, 0x2449), ARP_PATTERN_ACCENT(43281, 0x0801),
  ARP_PATTERN_ACCENT(28609, 0x0901), ARP_PATTERN_ACCENT(53505, 0x4001)
};

const uint16_t lut_res_arpeggiator_ties[ ARPEGGIATOR_PRESET_PATTERN_COUNT ] = {
  0x0041, 0x0001, 0x0001, 0x0101,
  0x0041, 0x0041, 0x1001, 0x0001,
  0x0001, 0x0001, 0x0001, 0x0101,
  0x0840, 0x1010, 0x0444, 0x0100,
  0x0440, 0x0048, 0x0400, 0x0001,
  0x0801, 0x0001
};

uint16_t Timbre::getArpeggiatorPattern() const
{
  const int pattern = (int)params_.engineArp2.pattern;
  if ( pattern < ARPEGGIATOR_PRESET_PATTERN_COUNT )
    return ARP_PATTERN_GETMASK(lut_res_arpeggiator_patterns[ pattern ]);
  else
    return ARP_PATTERN_GETMASK( params_.engineArpUserPatterns.patterns[ pattern - ARPEGGIATOR_PRESET_PATTERN_COUNT ] );
}

void Timbre::updateArpPresetCache(int pattern) const
{
  int seed = getArpEuclideanSeed(&params_);
  if (arpPresetCachePattern_ != (uint8_t)pattern || arpPresetCacheSeed_ != (uint8_t)seed) {
    computeArpPresetMasks(
        ARP_PATTERN_GETMASK(lut_res_arpeggiator_patterns[pattern]),
        ARP_PATTERN_GETACCENT(lut_res_arpeggiator_patterns[pattern]),
        lut_res_arpeggiator_ties[pattern],
        seed,
        &arpPresetCacheAcc_, &arpPresetCacheTie_, &arpPresetCacheOct_, &arpPresetCacheRt_,
        &arpPresetCacheOctDown_, &arpPresetCacheDyad_, &arpPresetCacheChord_, &arpPresetCacheSlap_);
    arpPresetCachePattern_ = (uint8_t)pattern;
    arpPresetCacheSeed_ = (uint8_t)seed;
  }
}

uint16_t Timbre::getArpeggiatorAccentMask() const
{
  const int pattern = (int)params_.engineArp2.pattern;
  if (pattern < ARPEGGIATOR_PRESET_PATTERN_COUNT) {
    updateArpPresetCache(pattern);
    return arpPresetCacheAcc_;
  }
  return ARP_PATTERN_GETACCENT(params_.engineArpUserPatterns.patterns[pattern - ARPEGGIATOR_PRESET_PATTERN_COUNT]);
}

int Timbre::getArpeggiatorPatternLength() const
{
  const int pattern = (int)params_.engineArp2.pattern;
  if (pattern < ARPEGGIATOR_PRESET_PATTERN_COUNT) {
    return getArpPresetLength(getArpEuclideanSeed(&params_));
  }
  return getArpUserPatternLength(&params_, pattern - ARPEGGIATOR_PRESET_PATTERN_COUNT);
}

int Timbre::getArpeggiatorStepMode() const
{
  const int pattern = (int)params_.engineArp2.pattern;
  if (bitmask_ == 0) {
    return ARP_STEP_MODE_NORM;
  }
  if (pattern < ARPEGGIATOR_PRESET_PATTERN_COUNT) {
    updateArpPresetCache(pattern);
    if ((bitmask_ & arpPresetCacheTie_) != 0) {
      return ARP_STEP_MODE_TIE;
    }
    if ((bitmask_ & arpPresetCacheOct_) != 0) {
      return ARP_STEP_MODE_OCT_UP;
    }
    if ((bitmask_ & arpPresetCacheRt_) != 0) {
      return ARP_STEP_MODE_RATCHET2;
    }
    if ((bitmask_ & arpPresetCacheOctDown_) != 0) {
      return ARP_STEP_MODE_OCT_DOWN;
    }
    if ((bitmask_ & arpPresetCacheDyad_) != 0) {
      return ARP_STEP_MODE_DYAD;
    }
    if ((bitmask_ & arpPresetCacheChord_) != 0) {
      return ARP_STEP_MODE_CHORD;
    }
    if ((bitmask_ & arpPresetCacheSlap_) != 0) {
      return ARP_STEP_MODE_SLAP;
    }
    return ARP_STEP_MODE_NORM;
  }
  int step = 0;
  uint16_t m = bitmask_;
  while ((m & 1u) == 0 && step < 15) {
    m >>= 1;
    step++;
  }
  return getArpUserStepMode(&params_, pattern - ARPEGGIATOR_PRESET_PATTERN_COUNT, step);
}

const uint8_t midi_clock_tick_per_step[17]  = {
  192, 144, 96, 72, 64, 48, 36, 32, 24, 16, 12, 8, 6, 4, 3, 2, 1
};

extern float noise[32];

const float panTable[]  = {
        0.0000, 0.0007, 0.0020, 0.0036, 0.0055, 0.0077, 0.0101, 0.0128, 0.0156, 0.0186,
        0.0218, 0.0252, 0.0287, 0.0324, 0.0362, 0.0401, 0.0442, 0.0484, 0.0527, 0.0572,
        0.0618, 0.0665, 0.0713, 0.0762, 0.0812, 0.0863, 0.0915, 0.0969, 0.1023, 0.1078,
        0.1135, 0.1192, 0.1250, 0.1309, 0.1369, 0.1430, 0.1492, 0.1554, 0.1618, 0.1682,
        0.1747, 0.1813, 0.1880, 0.1947, 0.2015, 0.2085, 0.2154, 0.2225, 0.2296, 0.2369,
        0.2441, 0.2515, 0.2589, 0.2664, 0.2740, 0.2817, 0.2894, 0.2972, 0.3050, 0.3129,
        0.3209, 0.3290, 0.3371, 0.3453, 0.3536, 0.3619, 0.3703, 0.3787, 0.3872, 0.3958,
        0.4044, 0.4131, 0.4219, 0.4307, 0.4396, 0.4485, 0.4575, 0.4666, 0.4757, 0.4849,
        0.4941, 0.5034, 0.5128, 0.5222, 0.5316, 0.5411, 0.5507, 0.5604, 0.5700, 0.5798,
        0.5896, 0.5994, 0.6093, 0.6193, 0.6293, 0.6394, 0.6495, 0.6597, 0.6699, 0.6802,
        0.6905, 0.7009, 0.7114, 0.7218, 0.7324, 0.7430, 0.7536, 0.7643, 0.7750, 0.7858,
        0.7967, 0.8076, 0.8185, 0.8295, 0.8405, 0.8516, 0.8627, 0.8739, 0.8851, 0.8964,
        0.9077, 0.9191, 0.9305, 0.9420, 0.9535, 0.9651, 0.9767, 0.9883, 1.0000, 1.0000,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000
} ;

Timbre::Timbre() {

    initFx();
    
    recomputeNext_ = true;
    currentGate_ = 0;
    sbMax_ = &sampleBlock_[64];
    holdPedal_ = false;
    lastPlayedNote_ = 0;
    // arpegiator
    setNewBPMValue(90);
    arpegiatorStep_ = 0.0;
    idle_ticks_ = 96;
    running_ = 0;
    ignore_note_off_messages_ = 0;
    recording_ = 0;
    extClockCounter_ = 0;
    extBlocksPerTick_ = 31.25f;
    arpRatchetSubBlocks_ = 0.0f;
    arpRatchetGateBlocks_ = 0.0f;
    arpRatchetTotalBlocks_ = 0.0f;
    arpRatchetElapsedBlocks_ = 0.0f;
    arpRatchetSubElapsed_ = 0.0f;
    arpRatchetCount_ = 0;
    arpRatchetVelocity_ = 0;
    arpRatchetHitIndex_ = 0;
    arpRatchetSubDiv_ = 1;
    arpRatchetSign_ = 1;
    arpRatchetGateOpen_ = false;
    zeroFadeL_ = 0.0f;
    zeroFadeR_ = 0.0f;
    zeroFadeStepL_ = 0.0f;
    zeroFadeStepR_ = 0.0f;
    zeroFadeBlocks_ = 0;
    arpRatchetActive_ = false;
    arpStrumCount_ = 0;
    arpPendingStepMode_ = ARP_STEP_MODE_NORM;
    arpSwingStepScale_ = 1.0f;
    arpPresetCachePattern_ = 0xFF;
    arpPresetCacheSeed_ = 0xFF;
    arpPresetCacheAcc_ = 0;
    arpPresetCacheTie_ = 0;
    arpPresetCacheOct_ = 0;
    arpPresetCacheRt_ = 0;
    arpPresetCacheOctDown_ = 0;
    arpPresetCacheDyad_ = 0;
    arpPresetCacheChord_ = 0;
    arpPresetCacheSlap_ = 0;
    note_stack_.Init();
    held_note_stack_.Init();
    voice_allocator_.Init(0, VOICE_PRIORITY_LAST_KEY, voice_allocator_pool_, voice_allocator_lru_);
    event_scheduler_.Init();
    // Arpeggiator start
    Start();

}

Timbre::~Timbre() {
}

void Timbre::init(SynthState *synthState, int timbreNumber) {

    mixerState_ = &synthState->mixerState;

    env1_.init(&params_.env1Time, &params_.env1Level, 0, &params_.engine1.algo, &params_.env1Curve);
    env2_.init(&params_.env2Time, &params_.env2Level, 1, &params_.engine1.algo, &params_.env2Curve);
    env3_.init(&params_.env3Time, &params_.env3Level, 2, &params_.engine1.algo, &params_.env3Curve);
    env4_.init(&params_.env4Time, &params_.env4Level, 3, &params_.engine1.algo, &params_.env4Curve);
    env5_.init(&params_.env5Time, &params_.env5Level, 4, &params_.engine1.algo, &params_.env5Curve);
    env6_.init(&params_.env6Time, &params_.env6Level, 5, &params_.engine1.algo, &params_.env6Curve);

    struct OperatorPhaseRowParams* phaseParamsBase = &params_.phaseOp1;
    osc1_.init(synthState, &params_.osc1, phaseParamsBase, OSC1_FREQ);
    osc2_.init(synthState, &params_.osc2, phaseParamsBase, OSC2_FREQ);
    osc3_.init(synthState, &params_.osc3, phaseParamsBase, OSC3_FREQ);
    osc4_.init(synthState, &params_.osc4, phaseParamsBase, OSC4_FREQ);
    osc5_.init(synthState, &params_.osc5, phaseParamsBase, OSC5_FREQ);
    osc6_.init(synthState, &params_.osc6, phaseParamsBase, OSC6_FREQ);

    timbreNumber_ = timbreNumber;

    for (int s = 0; s < 2; s++) {
        for (int n = 0; n < 128; n++) {
            midiNoteScale[s][timbreNumber][n] = INV127 * (float) n;
        }
    }
    for (int lfo = 0; lfo < NUMBER_OF_LFO; lfo++) {
        lfoUSed_[lfo] = 0;
    }

    lowerNote_ = 64;
    lowerNoteReleased_ = true;

    /** --------------FX init--------------  */

    delayBuffer_ = delayBuffer[timbreNumber_];

    for (int s = 0; s < delayBufferSize; s++) {
        delayBuffer_[s] = 0;
    }

    CounterpointTimbreState &cpt = cptTimbreBus[timbreNumber_];
    cpt.ticksEveryNCalls = 31.25f;
    cpt.arpegiatorStep = 0.0f;
    cpt.activeNote = 0;
    cpt.prevNote = 0;
    cpt.clockMode = 0;
    cpt.running = 0;
    cpt.tickedThisBlock = 0;
    cpt.stepTriggeredThisTick = 0;
    cpt.curPatternStep = 0;

}

void Timbre::setVoiceNumber(int v, int n) {
    voiceNumber_[v] = n;
    if (n >= 0) {
        voices_[n]->setCurrentTimbre(this);
    }
}

void Timbre::initVoicePointer(int n, Voice *voice) {
    voices_[n] = voice;
}

void Timbre::noteOn(char note, char velocity) {
    if (params_.engineArp1.clock) {
        arpeggiatorNoteOn(note, velocity);
    } else {
        uint8_t midiNote = (uint8_t) note;
        if (velocity > 0 && midiNote < 128 && numberOfVoices_ > 0
                && mixerState_->instrumentState_[timbreNumber_].scaleFrequencies[midiNote] != 0.0f) {
            held_note_stack_.NoteOn(midiNote, (uint8_t) velocity);
        }
        if (getEffectiveCptPart() != CPT_PART_ALL) {
            reconcileCounterpointDirectNotes();
        } else {
            preenNoteOn(note, velocity);
        }
    }
}

void Timbre::noteOff(char note) {
    if (params_.engineArp1.clock) {
        arpeggiatorNoteOff(note);
    } else {
        uint8_t previousSize = held_note_stack_.size();
        held_note_stack_.NoteOff((uint8_t) note);
        if (getEffectiveCptPart() != CPT_PART_ALL) {
            reconcileCounterpointDirectNotes();
        } else {
            preenNoteOff(note, held_note_stack_.size() != previousSize);
        }
    }
}

void Timbre::noteOnMPE(uint8_t channel, uint8_t note, uint8_t velocity)  {
	int voiceToUse = voiceNumber_[channel -1];

    if (unlikely(voiceToUse == -1)) {
        return;
    }

    preenNoteOnUpdateMatrix(voiceToUse, note, velocity);
    float noteFrequency = mixerState_->instrumentState_[0].scaleFrequencies[(int) note];

    if (voices_[voiceToUse]->isPlaying()) {
        voices_[voiceToUse]->noteOnWithoutPop(note, noteFrequency, velocity);
    } else {
        voices_[voiceToUse]->noteOn(note, noteFrequency, velocity);
    }
}


void Timbre::noteOffMPE(uint8_t channel, uint8_t note, uint8_t velocityOff) {
    int voiceToUse = voiceNumber_[channel -1];

    if (unlikely(voiceToUse == -1)) {
		return;
	}

	if (voices_[voiceToUse]->isPlaying()) {
		voices_[voiceToUse]->noteOff();
	}
}


void Timbre::preenNoteOn(char note, char velocity, bool restored) {
    // NumberOfVoice = 0 or no mapping in scala frequencies
    if (unlikely(numberOfVoices_ == 0 || mixerState_->instrumentState_[timbreNumber_].scaleFrequencies[(int) note] == 0.0f)) {
        return;
    }

    note &= 0x7f;
    {
        uint8_t uNote = (uint8_t) note;
        uint8_t prev = cptTimbreBus[timbreNumber_].activeNote;
        if (prev > 0 && uNote != prev) {
            cptTimbreBus[timbreNumber_].prevNote = prev;
        }
        cptTimbreBus[timbreNumber_].activeNote = uNote;
    }

    // Frequency depends on the current instrument scale
    float noteFrequency = mixerState_->instrumentState_[timbreNumber_].scaleFrequencies[(int) note];

    if (params_.engine1.playMode == PLAY_MODE_POLY) {
        uint8_t unavailable[MAX_NUMBER_OF_VOICES];
        for (int k = 0; k < numberOfVoices_; ++k) {
            unavailable[k] = voices_[voiceNumber_[k]]->isNewNotePending()
                    && (voice_allocator_pool_[k] & 0x80);
        }
        uint8_t voiceIndex = voice_allocator_.NoteOn((uint8_t) note, unavailable, restored);
        if (voiceIndex == 0xff) {
            voiceIndex = voice_allocator_.NoteOn((uint8_t) note, 0, restored);
            if (voiceIndex == 0xff) {
                return;
            }
        }

        int voiceToUse = voiceNumber_[voiceIndex];
        preenNoteOnUpdateMatrix(voiceToUse, note, velocity);
        if (voices_[voiceToUse]->isPlaying()) {
            voices_[voiceToUse]->noteOnWithoutPop(note, noteFrequency, velocity);
        } else {
            voices_[voiceToUse]->noteOn(note, noteFrequency, velocity);
        }
        if (isArpeggiatorActive() && voiceNumber_[0] >= 0 && voiceToUse != voiceNumber_[0]) {
            voices_[voiceToUse]->syncStepSeqFrom(voices_[voiceNumber_[0]]);
        }

        lastPlayedNote_ = voiceToUse;
        if (note <= lowerNote_ || lowerNoteReleased_) {
            lowerNote_ = note;
            lowerNoteReleased_ = false;
            lowerNoteFrequency = voices_[voiceToUse]->getNoteRealFrequencyEstimation(noteFrequency);
        }
        return;
    }

    int voiceToUse = voiceNumber_[0];
    if (unlikely(voiceToUse == -1)) {
        return;
    }

    uint8_t priorityMode = getVoicePriorityMode();
    if (!restored && held_note_stack_.size() > 1
            && (priorityMode == VOICE_PRIORITY_HIGH_KEY || priorityMode == VOICE_PRIORITY_LOW_KEY)
            && voices_[voiceToUse]->isPlaying() && !voices_[voiceToUse]->isReleased()
            && !voices_[voiceToUse]->isHoldedByPedal()) {
        uint8_t playingNote = (uint8_t) voices_[voiceToUse]->getNote() & 0x7f;
        uint8_t incomingNote = (uint8_t) note;
        if ((priorityMode == VOICE_PRIORITY_HIGH_KEY && incomingNote < playingNote)
                || (priorityMode == VOICE_PRIORITY_LOW_KEY && incomingNote > playingNote)) {
            return;
        }
    }

    bool isPlaying = voices_[voiceToUse]->isPlaying();
    if (likely(params_.engine1.playMode != PLAY_MODE_UNISON)) {
        preenNoteOnUpdateMatrix(voiceToUse, note, velocity);
        if (isPlaying) {
            voices_[voiceToUse]->noteOnWithoutPop(note, noteFrequency, velocity);
        } else {
            voices_[voiceToUse]->noteOn(note, noteFrequency, velocity);
        }
    } else {
        float noteFrequencyUnison = noteFrequency + noteFrequency * params_.engine2.unisonDetune * .05f;
        float noteFrequencyUnisonInc = noteFrequency * params_.engine2.unisonDetune * numberOfVoiceInverse_ * .1f;
        for (int k = 0; k < numberOfVoices_; k++) {
            int n = voiceNumber_[k];
            preenNoteOnUpdateMatrix(n, note, velocity);
            if (isPlaying) {
                voices_[n]->noteOnWithoutPop(note, noteFrequencyUnison, velocity, unisonPhase[k]);
            } else {
                voices_[n]->noteOn(note, noteFrequencyUnison, velocity, unisonPhase[k]);
            }
            noteFrequencyUnison += noteFrequencyUnisonInc;
        }
    }

    lastPlayedNote_ = voiceToUse;
    lowerNote_ = note;
    lowerNoteReleased_ = false;
    lowerNoteFrequency = voices_[voiceToUse]->getNoteRealFrequencyEstimation(noteFrequency);
}

void Timbre::preenNoteOnUpdateMatrix(int voiceToUse, int note, int velocity) {
    // Update voice matrix with midi note and velocity
    voices_[voiceToUse]->matrix.setSource(MATRIX_SOURCE_NOTE1, midiNoteScale[0][timbreNumber_][note]);
    voices_[voiceToUse]->matrix.setSource(MATRIX_SOURCE_NOTE2, midiNoteScale[1][timbreNumber_][note]);
    voices_[voiceToUse]->matrix.setSource(MATRIX_SOURCE_VELOCITY, INV127 * velocity);
    voices_[voiceToUse]->matrix.setSource(MATRIX_SOURCE_RANDOM, noise[voiceToUse]);
}

void Timbre::preenNoteOff(char note, bool restoreNotes) {
    if (unlikely(numberOfVoices_ == 0)) {
        return;
    }
    bool isUnison = params_.engine1.playMode == PLAY_MODE_UNISON;

    if (note == lowerNote_) {
        lowerNoteReleased_ = true;
    }
    if (cptTimbreBus[timbreNumber_].activeNote == (uint8_t) (note & 0x7f)) {
        cptTimbreBus[timbreNumber_].prevNote = cptTimbreBus[timbreNumber_].activeNote;
        cptTimbreBus[timbreNumber_].activeNote = 0;
    }

    if (params_.engine1.playMode == PLAY_MODE_POLY) {
        if (unlikely(holdPedal_)) {
            for (int k = 0; k < numberOfVoices_; k++) {
                int n = voiceNumber_[k];
                if (n >= 0 && voices_[n]->isPlaying() && !voices_[n]->isReleased() && voices_[n]->getNote() == note) {
                    voices_[n]->setHoldedByPedal(true);
                    return;
                }
            }
            return;
        }

        uint8_t releasedVoice = voice_allocator_.NoteOff((uint8_t) note);
        if (releasedVoice != 0xff) {
            int n = voiceNumber_[releasedVoice];
            if (likely(n >= 0)) {
                bool wasActive = voices_[n]->isPlaying() && !voices_[n]->isReleased();
                voices_[n]->noteOff();
                if (restoreNotes && wasActive) {
                    restoreHeldNotes(1);
                }
            }
        } else {
            for (int k = 0; k < numberOfVoices_; k++) {
                int n = voiceNumber_[k];
                if (n >= 0 && voices_[n]->isPlaying() && !voices_[n]->isReleased() && voices_[n]->getNote() == note) {
                    voices_[n]->noteOff();
                    if (restoreNotes) {
                        restoreHeldNotes(1);
                    }
                    break;
                }
            }
        }
        return;
    }

    // Unison we must check all voices.... (different from noteOn)
    int iNov = params_.engine1.playMode == PLAY_MODE_MONO ? 1 : (int) numberOfVoices_;
    for (int k = 0; k < iNov; k++) {
        // voice number k of timbre
        int n = voiceNumber_[k];

        if (unlikely(n < 0 || !voices_[n]->isPlaying() || voices_[n]->isReleased())) {
            continue;
        }

        if (voices_[n]->getNote() == note) {
            if (unlikely(holdPedal_)) {
                voices_[n]->setHoldedByPedal(true);
                if (likely(!isUnison)) {
                    return;
                }
            } else {
                if (restoreNotes && held_note_stack_.size() > 0) {
                    const NoteEntry& previousNote = getHeldNoteForRestore();
                    preenNoteOn(previousNote.note, previousNote.velocity, true);
                    return;
                }
                voices_[n]->noteOff();
                if (likely(!isUnison)) {
                    return;
                }
            }
        }
    }
}

const NoteEntry& Timbre::getHeldNoteForRestore() const {
    uint8_t priorityMode = getVoicePriorityMode();
    if (priorityMode == VOICE_PRIORITY_HIGH_KEY) {
        return held_note_stack_.sorted_note(held_note_stack_.size() - 1);
    }
    if (priorityMode == VOICE_PRIORITY_LOW_KEY) {
        return held_note_stack_.sorted_note(0);
    }
    return held_note_stack_.most_recent_note();
}

void Timbre::restoreHeldNotes(uint8_t numberOfVoices) {
    uint8_t priorityMode = getVoicePriorityMode();
    while (numberOfVoices > 0) {
        bool restored = false;
        uint8_t heldCount = held_note_stack_.size();
        for (uint8_t k = 0; k < heldCount; k++) {
            const NoteEntry& candidate =
                (priorityMode == VOICE_PRIORITY_HIGH_KEY || priorityMode == VOICE_PRIORITY_BASS_PROTECT)
                    ? held_note_stack_.sorted_note(heldCount - k - 1)
                    : (priorityMode == VOICE_PRIORITY_LOW_KEY
                        ? held_note_stack_.sorted_note(k)
                        : held_note_stack_.recent_note(k));
            if (voice_allocator_.FindActive(candidate.note) == 0xff) {
                preenNoteOn(candidate.note, candidate.velocity, true);
                if (voice_allocator_.FindActive(candidate.note) != 0xff) {
                    --numberOfVoices;
                    restored = true;
                    break;
                }
            }
        }
        if (!restored) {
            return;
        }
    }
}


void Timbre::stopPlayingNow() {
    for (int k = 0; k < numberOfVoices_; k++) {
        // voice number k of timbre
        int n = voiceNumber_[k];
        if (n != -1) {
            voices_[n]->killNow();
        }
    }
}

void Timbre::setHoldPedal(int value) {
    if (value < 64) {
        holdPedal_ = false;
        int nVoices = numberOfVoices_;
        int releasedVoices = 0;
        for (int k = 0; k < nVoices; k++) {
            // voice number k of timbre
            int n = voiceNumber_[k];
            if (voices_[n]->isHoldedByPedal()) {
                if (params_.engine1.playMode == PLAY_MODE_POLY) {
                    voice_allocator_.NoteOff((uint8_t) voices_[n]->getNote());
                }
                voices_[n]->noteOff();
                if (voices_[n]->isReleased()) {
                    releasedVoices++;
                }
            }
        }
        if (params_.engine1.playMode == PLAY_MODE_POLY) {
            restoreHeldNotes((uint8_t) releasedVoices);
        } else if (releasedVoices > 0 && held_note_stack_.size() > 0) {
            const NoteEntry& previousNote = getHeldNoteForRestore();
            preenNoteOn(previousNote.note, previousNote.velocity, true);
        }
        arpeggiatorSetHoldPedal(0);
    } else {
        holdPedal_ = true;
        arpeggiatorSetHoldPedal(127);
    }
}

void Timbre::setNewBPMValue(float bpm) {
    float ticksPerSecond = bpm * 24.0f / 60.0f;
    ticksEveryNCalls_ = CALLED_PER_SECOND / ticksPerSecond;
}

void Timbre::setArpeggiatorClock(float clockValue) {
    if (clockValue == CLOCK_OFF) {
        FlushQueue();
        note_stack_.Clear();
    }
    if (clockValue == CLOCK_INTERNAL) {
        setNewBPMValue(params_.engineArp1.BPM);
    }
    if (clockValue == CLOCK_EXTERNAL) {
        // Let's consider we're running
        running_ = 1;
    }
}

void Timbre::updateArpegiatorInternalClock() {
    CounterpointTimbreState &myCpt = cptTimbreBus[timbreNumber_];
    myCpt.tickedThisBlock = 0;
    myCpt.stepTriggeredThisTick = 0;
    myCpt.clockMode = (uint8_t) params_.engineArp1.clock;
    myCpt.running = running_;

    // Arpeggiator clock : internal
    if (params_.engineArp1.clock == CLOCK_INTERNAL) {
        int rank = 0, total = 1, masterTimbre = timbreNumber_;
        bool syncToChannelMaster = (getCptPartMode() != CPT_PART_ALL);
        if (syncToChannelMaster) {
            getChannelStackInfo(&rank, &total, &masterTimbre);
        }
        if (syncToChannelMaster
                && masterTimbre >= 0 && masterTimbre < timbreNumber_
                && cptTimbreBus[masterTimbre].clockMode == CLOCK_INTERNAL
                && cptTimbreBus[masterTimbre].running != 0) {
            ticksEveryNCalls_ = cptTimbreBus[masterTimbre].ticksEveryNCalls;
            arpegiatorStep_ = cptTimbreBus[masterTimbre].arpegiatorStep;
            if (unlikely(cptTimbreBus[masterTimbre].tickedThisBlock != 0)) {
                if (cptTimbreBus[masterTimbre].stepTriggeredThisTick != 0
                        && cptTimbreBus[masterTimbre].curPatternStep == 0) {
                    int stepTicks = getArpStepTicks();
                    if (stepTicks > 0) {
                        tick_ = (uint8_t) (stepTicks - 1);
                    }
                }
                myCpt.tickedThisBlock = 1;
                Tick();
            }
        } else {
            arpegiatorStep_ += 1.0f;
            if (unlikely((arpegiatorStep_) > ticksEveryNCalls_)) {
                arpegiatorStep_ -= ticksEveryNCalls_;
                myCpt.tickedThisBlock = 1;
                Tick();
            }
        }
        myCpt.ticksEveryNCalls = ticksEveryNCalls_;
        myCpt.arpegiatorStep = arpegiatorStep_;
    } else if (params_.engineArp1.clock == CLOCK_EXTERNAL) {
        if (extClockCounter_ < 2000) {
            extClockCounter_++;
        }
    }
    if (unlikely(arpRatchetActive_)) {
        updateArpRatchet();
    }
    if (unlikely(arpStrumCount_ > 0)) {
        updateArpStrum();
    }
}

void Timbre::cleanNextBlock() {

    float *sp = sampleBlock_;
    while (sp < sbMax_) {
        *sp++ = 0;
        *sp++ = 0;
        *sp++ = 0;
        *sp++ = 0;
        *sp++ = 0;
        *sp++ = 0;
        *sp++ = 0;
        *sp++ = 0;
    }
}

void Timbre::prepareMatrixForNewBlock() {
    for (int k = 0; k < numberOfVoices_; k++) {
        voices_[voiceNumber_[k]]->prepareMatrixForNewBlock();
    }
    if (isArpeggiatorActive() && voiceNumber_[0] >= 0) {
        voices_[voiceNumber_[0]]->updateArpModulatorsInMatrix();
    }
}


void Timbre::glide() {
    if (unlikely(params_.engine1.playMode == PLAY_MODE_UNISON)) {
        for (int vv = 0; vv < numberOfVoices_; vv++) {
            int v = voiceNumber_[vv];
            if (v != -1 && voices_[v]->isGliding()) {
                voices_[v]->glide();
            }
        }
    } else {
        if (voiceNumber_[0] != -1 && voices_[voiceNumber_[0]]->isGliding()) {
            voices_[voiceNumber_[0]]->glide();
        }
    }
}

uint8_t Timbre::voicesNextBlock() {
    uint8_t numberOfPlayingVoices_ = 0;
    if (unlikely(params_.engine1.playMode == PLAY_MODE_UNISON)) {
        cleanNextBlock();
        // UNISON
        float pansSav[6];
        pansSav[0] = params_.engineMix1.panOsc1;
        pansSav[1] = params_.engineMix1.panOsc2;
        pansSav[2] = params_.engineMix2.panOsc3;
        pansSav[3] = params_.engineMix2.panOsc4;
        pansSav[4] = params_.engineMix3.panOsc5;
        pansSav[5] = params_.engineMix3.panOsc6;

        float* pans[] = {
            &params_.engineMix1.panOsc1,
            &params_.engineMix1.panOsc2,
            &params_.engineMix2.panOsc3,
            &params_.engineMix2.panOsc4,
            &params_.engineMix3.panOsc5,
            &params_.engineMix3.panOsc6
        };

        int currentAlgo = (int)params_.engine1.algo;
        float algoNumberOfMix = algoInformation[currentAlgo].mix;
        float numberOfCarrierOp = numberOfVoices_ * algoNumberOfMix;
        float opPan = - params_.engine2.unisonSpread;
        float opPanInc = 2.0f / numberOfCarrierOp * params_.engine2.unisonSpread;

        if (likely(voiceNumber_[0] >= 0 && voices_[voiceNumber_[0]]->isPlaying())) {
            for (int vv = 0; vv < numberOfVoices_; vv++) {
                int v = voiceNumber_[vv];
                if (unlikely(v < 0)) {
                    continue;
                }
                bool otherSide = false;
                for (int op = 0; op < 6; op ++) {
                    if (algoOpInformation[currentAlgo][op] == 1) {
                        if (otherSide) {
                            *pans[op] = -opPan;
                        } else {
                            *pans[op] = opPan;
                        }
                        opPan += opPanInc;
                        otherSide =! otherSide;
                    }
                }
                voices_[v]->nextBlock();

                if (vv > 0) {
                    // We accumulate in the first voice buffer
                    float *source = voices_[v]->getSampleBlock();
                    float *voice0SampleBlock = voices_[voiceNumber_[0]]->getSampleBlock();
                    for (int s = 0; s < BLOCK_SIZE; s++) {
                        *voice0SampleBlock++ += *source++;
                        *voice0SampleBlock++ += *source++;
                    }
                }
                numberOfPlayingVoices_++;
            }
            voices_[voiceNumber_[0]]->fxAfterBlock();

        } else if (likely(voiceNumber_[0] >= 0)) {
            voices_[voiceNumber_[0]]->emptyBuffer();
        }

        params_.engineMix1.panOsc1 = pansSav[0];
        params_.engineMix1.panOsc2 = pansSav[1];
        params_.engineMix2.panOsc3 = pansSav[2];
        params_.engineMix2.panOsc4 = pansSav[3];
        params_.engineMix3.panOsc5 = pansSav[4];
        params_.engineMix3.panOsc6 = pansSav[5];
    } else {
        for (int k = 0; k < numberOfVoices_; k++) {
            int v = voiceNumber_[k];
            if (unlikely(v < 0)) {
                continue;
            }
            if (likely(voices_[v]->isPlaying())) {
                voices_[v]->nextBlock();
                voices_[v]->fxAfterBlock();
                numberOfPlayingVoices_++;
            } else {
                voices_[v]->emptyBuffer();
            }
        }
    }

    return numberOfPlayingVoices_;
}

void Timbre::voicesToTimbre(float volumeGain) {

    int numberOfVoicesToCopy = numberOfVoices_;
    if (unlikely(numberOfVoicesToCopy <= 0 || voiceNumber_[0] < 0)) {
        cleanNextBlock();
        return;
    }

    // If unison : all voices have already been added to the first one
    if (unlikely(params_.engine1.playMode == PLAY_MODE_UNISON)) {
        numberOfVoicesToCopy = 1;
    }

    for (int k = 0; k < numberOfVoicesToCopy; k++) {
        float *timbreBlock = sampleBlock_;
        const float *voiceBlock = voices_[voiceNumber_[k]]->getSampleBlock();

        if (unlikely(k == 0)) {
            if (unlikely(numberOfVoicesToCopy == 1)) {
                for (int s = 0; s < BLOCK_SIZE; s++) {
                    *timbreBlock++ = *voiceBlock++ * volumeGain;
                    *timbreBlock++ = *voiceBlock++ * volumeGain;
                }
            } else {
                for (int s = 0; s < BLOCK_SIZE; s++) {
                    *timbreBlock++ = *voiceBlock++;
                    *timbreBlock++ = *voiceBlock++;
                }
            }
        } else if (k == numberOfVoicesToCopy - 1) {
            for (int s = 0; s < BLOCK_SIZE; s++) {
                *timbreBlock = (*timbreBlock + *voiceBlock++) * volumeGain;
                timbreBlock++;
                *timbreBlock = (*timbreBlock + *voiceBlock++) * volumeGain;
                timbreBlock++;
            }
        } else {
            for (int s = 0; s < BLOCK_SIZE; s++) {
                *timbreBlock++ += *voiceBlock++;
                *timbreBlock++ += *voiceBlock++;
            }
        }
    }
    if (unlikely(isZeroFadeActive())) {
        stepZeroFade(true);
    }
}

void Timbre::gateFx() {
    // Gate algo !!
    float gate = 0.0f;
    if (likely(params_.engine1.playMode == PLAY_MODE_POLY)) {
        for (int k = 0; k < numberOfVoices_; k++) {
            int v = voiceNumber_[k];
            if (unlikely(v < 0 || !voices_[v]->isPlaying())) {
                continue;
            }
            float voiceGate = voices_[v]->matrix.getDestination(MAIN_GATE);
            if (voiceGate > gate) {
                gate = voiceGate;
            }
        }
    } else {
        gate = voices_[lastPlayedNote_]->matrix.getDestination(MAIN_GATE);
    }
    if (unlikely(gate > 0 || currentGate_ > 0)) {
        gate *= .72547132656922730694f; // 0 < gate < 1.0
        if (gate > 1.0f) {
            gate = 1.0f;
        }
        float incGate = (gate - currentGate_) * .03125f; // ( *.03125f = / 32)
        // limit the speed.
        if (incGate > 0.002f) {
            incGate = 0.002f;
        } else if (incGate < -0.002f) {
            incGate = -0.002f;
        }

        float *sp = sampleBlock_;
        float coef;
        for (int k = 0; k < BLOCK_SIZE; k++) {
            currentGate_ += incGate;
            coef = 1.0f - currentGate_;
            *sp = *sp * coef;
            sp++;
            *sp = *sp * coef;
            sp++;
        }
        //    currentGate = gate;
    }
}

void Timbre::afterNewParamsLoad() {

    env1_.applyCurves();
    env2_.applyCurves();
    env3_.applyCurves();
    env4_.applyCurves();
    env5_.applyCurves();
    env6_.applyCurves();

    for (int k = 0; k < numberOfVoices_; k++) {
        voices_[voiceNumber_[k]]->afterNewParamsLoad();
    }

    for (int j = 0; j < NUMBER_OF_ENCODERS_PFM2; j++) {
        env1_.reloadADSR(j);
        env2_.reloadADSR(j);
        env3_.reloadADSR(j);
        env4_.reloadADSR(j);
        env5_.reloadADSR(j);
        env6_.reloadADSR(j);
    }

    resetArpeggiator();
    for (int k = 0; k < NUMBER_OF_ENCODERS_PFM2; k++) {
        setNewEffecParam(k);
    }
    // Update midi note scale
    updateMidiNoteScale(0);
    updateMidiNoteScale(1);
    updateVoicePriorityMode();
}

void Timbre::resetArpeggiator() {
    // Reset Arpeggiator
    FlushQueue();
    note_stack_.Clear();
    setArpeggiatorClock(params_.engineArp1.clock);
    setLatchMode(params_.engineArp2.latche);
}

void Timbre::setNewValue(int index, struct ParameterDisplay *param, float newValue) {
    if (newValue > param->maxValue) {
        // in v2, matrix target were removed so some values are > to max value but we need to accept it
        bool mustConstraint = true;
        if (param->valueNameOrder != NULL) {
            for (int v = 0; v < param->numberOfValues; v++) {
                if ((int) param->valueNameOrder[v] == (int) (newValue + .01)) {
                    mustConstraint = false;
                }
            }
        }
        if (mustConstraint) {
            newValue = param->maxValue;
        }
    } else if (newValue < param->minValue) {
        newValue = param->minValue;
    }
    ((float*) &params_)[index] = newValue;
}

int Timbre::getSeqStepValue(int whichStepSeq, int step) {

    if (whichStepSeq == 0) {
        return params_.lfoSteps1.steps[step] & 0x0F;
    } else {
        return params_.lfoSteps2.steps[step] & 0x0F;
    }
}

void Timbre::setSeqStepValue(int whichStepSeq, int step, int value) {

    if (whichStepSeq == 0) {
        params_.lfoSteps1.steps[step] = (char) ((params_.lfoSteps1.steps[step] & 0xF0) | (value & 0x0F));
    } else {
        params_.lfoSteps2.steps[step] = (char) ((params_.lfoSteps2.steps[step] & 0xF0) | (value & 0x0F));
    }
}

int Timbre::getSeqStepMode(int whichStepSeq, int step) {
    uint8_t raw = (whichStepSeq == 0) ? (uint8_t) params_.lfoSteps1.steps[step] : (uint8_t) params_.lfoSteps2.steps[step];
    int mode = (raw >> 4) & 0x0F;
    return (mode < STEPSEQ_MODE_COUNT) ? mode : STEPSEQ_MODE_NORM;
}

void Timbre::setSeqStepMode(int whichStepSeq, int step, int mode) {
    if (mode < 0 || mode >= STEPSEQ_MODE_COUNT) {
        return;
    }
    if (whichStepSeq == 0) {
        params_.lfoSteps1.steps[step] = (char) (((mode & 0x0F) << 4) | (params_.lfoSteps1.steps[step] & 0x0F));
    } else {
        params_.lfoSteps2.steps[step] = (char) (((mode & 0x0F) << 4) | (params_.lfoSteps2.steps[step] & 0x0F));
    }
}

void Timbre::setNewEffecParam(int encoder) {

    for (int k = 0; k < numberOfVoices_; k++) {
        voices_[voiceNumber_[k]]->setNewEffectParam(encoder);
    }

}

// Code bellowed have been adapted by Xavier Hosxe for PreenFM2
// It come from Muteable Instrument midiPAL

/////////////////////////////////////////////////////////////////////////////////
// Copyright 2011 Emilie Gillet.
//
// Author: Emilie Gillet (emilie.o.gillet@gmail.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// -----------------------------------------------------------------------------
//
// Arpeggiator app.

void Timbre::arpeggiatorNoteOn(char note, char velocity) {
    bool isNewLatchedChord = (latch_ && !recording_);
    bool wasEmpty = (note_stack_.size() == 0) || isNewLatchedChord;
    // CLOCK_MODE_INTERNAL
    if (params_.engineArp1.clock == CLOCK_INTERNAL) {
        if (wasEmpty || idle_ticks_ >= 96 || !running_) {
            Start();
        }
        idle_ticks_ = 0;
    } else if (wasEmpty) {
        resetArpPatternState();
    }

    if (isNewLatchedChord) {
        note_stack_.Clear();
        recording_ = 1;
    }
    note_stack_.NoteOn(note, velocity);
}

void Timbre::arpeggiatorNoteOff(char note) {
    if (ignore_note_off_messages_) {
        return;
    }
    if (!latch_) {
        note_stack_.NoteOff(note);
        if (note_stack_.size() == 0) {
            resetArpPatternState();
        }
    } else {
        if (note == note_stack_.most_recent_note().note) {
            recording_ = 0;
        }
    }
}

void Timbre::OnMidiContinue() {
    if (params_.engineArp1.clock == CLOCK_EXTERNAL) {
        running_ = 1;
    }
}

void Timbre::OnMidiStart() {
    if (params_.engineArp1.clock == CLOCK_EXTERNAL) {
        Start();
    }
}

void Timbre::OnMidiStop() {
    if (params_.engineArp1.clock == CLOCK_EXTERNAL) {
        running_ = 0;
        FlushQueue();
    }
}

void Timbre::OnMidiClock() {
    if (params_.engineArp1.clock == CLOCK_EXTERNAL) {
        if (extClockCounter_ > 2 && extClockCounter_ < 2000) {
            extBlocksPerTick_ = (float) extClockCounter_;
        }
        extClockCounter_ = 0;
        if (running_) {
            Tick();
        }
    }
}

void Timbre::updateArpStrum() {
    uint8_t remaining = 0;
    bool stoppedPrevRatchet = false;
    for (uint8_t i = 0; i < arpStrumCount_; i++) {
        if (arpStrumDelay_[i] > 0) {
            arpStrumDelay_[i]--;
        }
        if (arpStrumDelay_[i] == 0) {
            if (!stoppedPrevRatchet) {
                stopArpRatchet();
                stoppedPrevRatchet = true;
            }
            SendNote(arpStrumNotes_[i], arpStrumVelocities_[i], (int) arpPendingStepMode_);
        } else {
            if (remaining != i) {
                arpStrumNotes_[remaining] = arpStrumNotes_[i];
                arpStrumVelocities_[remaining] = arpStrumVelocities_[i];
                arpStrumDelay_[remaining] = arpStrumDelay_[i];
            }
            remaining++;
        }
    }
    arpStrumCount_ = remaining;
}

void Timbre::stopArpRatchet() {
    if (arpRatchetActive_) {
        if (arpRatchetGateOpen_) {
            for (uint8_t i = 0; i < arpRatchetCount_; i++) {
                preenNoteOff(arpRatchetNotes_[i]);
            }
            arpRatchetGateOpen_ = false;
        }
        arpRatchetActive_ = false;
        arpRatchetCount_ = 0;
    }
}

static inline float getRatchetDensityScale(float subBlocks) {
    if (subBlocks >= 90.0f) {
        return 1.0f;
    }
    float scale = 0.70f + 0.30f * ((subBlocks - 42.0f) / 48.0f);
    return (scale < 0.70f) ? 0.70f : scale;
}

void Timbre::updateArpRatchet() {
    arpRatchetElapsedBlocks_ += 1.0f;
    arpRatchetSubElapsed_ += 1.0f;

    if (arpRatchetGateOpen_ && arpRatchetSubElapsed_ >= arpRatchetGateBlocks_) {
        for (uint8_t i = 0; i < arpRatchetCount_; i++) {
            preenNoteOff(arpRatchetNotes_[i]);
        }
        arpRatchetGateOpen_ = false;
    }

    if (arpRatchetElapsedBlocks_ >= (arpRatchetTotalBlocks_ - 0.5f)) {
        stopArpRatchet();
        return;
    }

    if (arpRatchetSubElapsed_ >= arpRatchetSubBlocks_) {
        arpRatchetSubElapsed_ -= arpRatchetSubBlocks_;
        if ((arpRatchetTotalBlocks_ - arpRatchetElapsedBlocks_) >= (arpRatchetSubBlocks_ * 0.4f)) {
            arpRatchetHitIndex_++;
            float progress = (arpRatchetSubDiv_ > 1) ? ((float) arpRatchetHitIndex_ / (float) (arpRatchetSubDiv_ - 1)) : 1.0f;
            if (progress > 1.0f) {
                progress = 1.0f;
            }
            float rampScale = (arpRatchetSign_ >= 0) ? (0.25f + 0.75f * progress) : (1.0f - 0.75f * progress);
            float densityScale = getRatchetDensityScale(arpRatchetSubBlocks_);

            int effVel = (int) arpRatchetVelocity_;
            if (voiceNumber_[0] >= 0) {
                effVel += (int) roundf(voices_[voiceNumber_[0]]->matrix.getDestination(MTX_DEST_ARP_VEL) * 64.0f);
            }
            effVel = (int) roundf((float) effVel * rampScale * densityScale);
            if (effVel < 1) {
                effVel = 1;
            } else if (effVel > 127) {
                effVel = 127;
            }
            uint8_t hitVelocity = (uint8_t) effVel;
            for (uint8_t i = 0; i < arpRatchetCount_; i++) {
                preenNoteOn(arpRatchetNotes_[i], hitVelocity);
            }
            arpRatchetGateOpen_ = true;
        } else {
            stopArpRatchet();
        }
    }
}

int Timbre::getArpStepTicks() const {
    int divIdx = (int) params_.engineArp2.division;
    if (voiceNumber_[0] >= 0) {
        divIdx += (int) roundf(voices_[voiceNumber_[0]]->matrix.getDestination(MTX_DEST_ARP_DIV) * 12.0f);
    }
    if (divIdx < 0) {
        divIdx = 0;
    } else if (divIdx > 16) {
        divIdx = 16;
    }
    return midi_clock_tick_per_step[divIdx];
}

void Timbre::SendNote(uint8_t note, uint8_t velocity, int stepMode) {
    Matrix *arpMatrix = (voiceNumber_[0] >= 0) ? &voices_[voiceNumber_[0]]->matrix : NULL;
    uint8_t baseVelocity = velocity;
    int ratchetMode = 0;
    int8_t ratchetSign = 1;
    if (arpMatrix != NULL) {
        int effNote = (int) note + (int) roundf(arpMatrix->getDestination(MTX_DEST_ARP_NOTE) * 12.0f);
        if (effNote < 0) {
            effNote = 0;
        } else if (effNote > 127) {
            effNote = 127;
        }
        note = (uint8_t) effNote;

        int effVel = (int) velocity + (int) roundf(arpMatrix->getDestination(MTX_DEST_ARP_VEL) * 64.0f);
        if (effVel < 1) {
            effVel = 1;
        } else if (effVel > 127) {
            effVel = 127;
        }
        velocity = (uint8_t) effVel;

        float ratchetVal = arpMatrix->getDestination(MTX_DEST_ARP_RATCHET);
        ratchetSign = (ratchetVal < 0.0f) ? -1 : 1;
        float absRatchet = fabsf(ratchetVal);
        if (absRatchet >= 0.04f) {
            float curved = 0.25f * absRatchet + 0.75f * absRatchet * absRatchet * absRatchet;
            ratchetMode = (int) roundf(curved * 7.0f);
            if (ratchetMode < 1) {
                ratchetMode = 1;
            } else if (ratchetMode > 11) {
                ratchetMode = 11;
            }
        }
    }

    // If there are some Note Off messages for the note about to be triggeered
    // remove them from the queue and process them now.
    if (event_scheduler_.Remove(note, 0)) {
        preenNoteOff(note);
    }

    // Send a note on and schedule a note off later.
    if (stepMode < 0) {
        stepMode = getArpeggiatorStepMode();
    }
    if (stepMode == ARP_STEP_MODE_RATCHET2 && ratchetMode == 0) {
        ratchetMode = 1;
    } else if (stepMode == ARP_STEP_MODE_SLAP) {
        int slapVel = (int) velocity + ((int) velocity >> 2);
        velocity = (uint8_t) ((slapVel > 127) ? 127 : slapVel);
    }

    int rawDurTicks = (int) midi_clock_tick_per_step[(int) params_.engineArp2.duration];
    int stepTicks = getArpStepTicks();
    if (arpMatrix != NULL) {
        rawDurTicks += (int) roundf(arpMatrix->getDestination(MTX_DEST_ARP_GATE) * (float) stepTicks);
    }
    if (stepMode == ARP_STEP_MODE_TIE && rawDurTicks < (stepTicks + 1)) {
        rawDurTicks = stepTicks + 1;
    } else if (stepMode == ARP_STEP_MODE_SLAP) {
        rawDurTicks = (int) roundf((float) rawDurTicks * 0.25f);
        if (rawDurTicks < 2) {
            rawDurTicks = 2;
        } else if (rawDurTicks > 7) {
            rawDurTicks = 7;
        }
    }

    int subDiv = 1;
    float totalBlocks = 0.0f;
    float subBlocks = 0.0f;
    float densityScale = 1.0f;
    if (ratchetMode > 0) {
        // Keep ratchet spacing >= 42 audio blocks (~28 ms at 1500 Hz) so transients
        // do not stack into audio-rate buzz.
        const float minSubBlocks = 42.0f;
        float blocksPerTick = (params_.engineArp1.clock == CLOCK_EXTERNAL)
            ? extBlocksPerTick_
            : ((ticksEveryNCalls_ > 2.0f) ? ticksEveryNCalls_ : 31.25f);
        totalBlocks = (float) stepTicks * blocksPerTick * arpSwingStepScale_;
        int maxSubDiv = (int) floorf(totalBlocks / minSubBlocks);
        subDiv = ratchetMode + 1;
        if (subDiv > maxSubDiv) {
            subDiv = maxSubDiv;
        }
        if (subDiv < 2) {
            ratchetMode = 0;
        } else {
            subBlocks = totalBlocks / (float) subDiv;
            densityScale = getRatchetDensityScale(subBlocks);
        }
    }

    uint8_t firstHitVel = velocity;
    if (ratchetMode > 0) {
        float firstScale = (ratchetSign > 0) ? (0.25f * densityScale) : densityScale;
        int scaled = (int) roundf((float) velocity * firstScale);
        firstHitVel = (uint8_t) ((scaled < 1) ? 1 : scaled);
    }
    preenNoteOn(note, firstHitVel);

    if (ratchetMode > 0) {
        // Ratchet subdivisions of the current arpeggiator step (stepTicks):
        // 1..7 -> x2..x8 subdivisions (and up to x12 when Mult > 1.0), bounded by minSubBlocks.
        arpRatchetTotalBlocks_ = totalBlocks;
        arpRatchetSubBlocks_ = subBlocks;
        float gateRatio = (stepTicks > 0) ? ((float) rawDurTicks / (float) stepTicks) : 0.75f;
        if (gateRatio < 0.35f) {
            gateRatio = 0.35f;
        } else if (gateRatio > 0.85f) {
            gateRatio = 0.85f;
        }
        arpRatchetGateBlocks_ = arpRatchetSubBlocks_ * gateRatio;
        if (arpRatchetGateBlocks_ < 10.0f) {
            arpRatchetGateBlocks_ = 10.0f;
        }
        if (arpRatchetGateBlocks_ > (arpRatchetSubBlocks_ - 4.0f)) {
            arpRatchetGateBlocks_ = arpRatchetSubBlocks_ - 4.0f;
        }
        arpRatchetElapsedBlocks_ = 0.0f;
        arpRatchetSubElapsed_ = 0.0f;
        if (arpRatchetCount_ < 8) {
            arpRatchetNotes_[arpRatchetCount_++] = note;
        }
        arpRatchetVelocity_ = baseVelocity;
        arpRatchetHitIndex_ = 0;
        arpRatchetSubDiv_ = (uint8_t) subDiv;
        arpRatchetSign_ = ratchetSign;
        arpRatchetGateOpen_ = true;
        arpRatchetActive_ = true;
        return;
    }

    int durTicks = rawDurTicks - 1;
    if (arpSwingStepScale_ < 1.0f && durTicks > 1) {
        durTicks = (int) roundf((float) durTicks * arpSwingStepScale_);
    }
    if (durTicks < 1) {
        durTicks = 1;
    } else if (durTicks > 192) {
        durTicks = 192;
    }
    event_scheduler_.Schedule(note, 0, (uint8_t) durTicks, 0);
}

void Timbre::SendLater(uint8_t note, uint8_t velocity, uint8_t when, uint8_t tag) {
    event_scheduler_.Schedule(note, velocity, when, tag);
}

void Timbre::SendScheduledNotes() {
    uint8_t current = event_scheduler_.root();
    while (current) {
        const SchedulerEntry &entry = event_scheduler_.entry(current);
        if (entry.when) {
            break;
        }
        if (entry.note != kZombieSlot) {
            if (entry.velocity == 0) {
                preenNoteOff(entry.note);
            } else {
                preenNoteOn(entry.note, entry.velocity);
            }
        }
        current = entry.next;
    }
    event_scheduler_.Tick();
}

void Timbre::FlushQueue() {
    stopArpRatchet();
    arpStrumCount_ = 0;
    arpPendingStepMode_ = ARP_STEP_MODE_NORM;
    arpSwingStepScale_ = 1.0f;
    while (event_scheduler_.size()) {
        SendScheduledNotes();
    }
}

void Timbre::Tick() {
    ++tick_;

    if (note_stack_.size()) {
        idle_ticks_ = 0;
    }
    ++idle_ticks_;
    if (idle_ticks_ >= 96) {
        idle_ticks_ = 96;
        if (params_.engineArp1.clock == CLOCK_INTERNAL) {
            running_ = 0;
            FlushQueue();
        }
    }

    SendScheduledNotes();

    int stepTicks = getArpStepTicks();
    if (tick_ >= stepTicks) {
        tick_ = 0;
        arpStrumCount_ = 0;

        uint16_t pattern = getArpeggiatorPattern();
        uint16_t accentMask = getArpeggiatorAccentMask();
        int patternLen = getArpeggiatorPatternLength();
        uint16_t endMask = (uint16_t) (1u << (patternLen - 1));
        if (bitmask_ == 0 || bitmask_ > endMask) {
            bitmask_ = 1;
        }
        int curPatternStep = 0;
        {
            uint16_t m = bitmask_;
            while ((m & 1u) == 0 && curPatternStep < 15) {
                m >>= 1;
                curPatternStep++;
            }
        }
        int nextPatternStep = (curPatternStep + 1 < patternLen) ? (curPatternStep + 1) : 0;

        int grooveIdx = 0, humanize = 0;
        getArpSwingParams(&params_, &grooveIdx, &humanize);
        bool curSwung = isArpSwingStep(curPatternStep, grooveIdx);
        bool nextSwung = isArpSwingStep(nextPatternStep, grooveIdx);

        float swingNorm = params_.midiNote1Curve.unused1 * 0.01f;
        if (swingNorm < 0.0f) {
            swingNorm = 0.0f;
        } else if (swingNorm > 0.50f) {
            swingNorm = 0.50f;
        }
        if (swingNorm > 0.0f) {
            float curOffset = curSwung ? swingNorm : 0.0f;
            float nextOffset = nextSwung ? swingNorm : 0.0f;
            arpSwingStepScale_ = 1.0f + nextOffset - curOffset;
        } else {
            arpSwingStepScale_ = 1.0f;
        }
        float swingVelScale = (curSwung && swingNorm > 0.0f)
            ? (1.0f - swingNorm * getArpSwingVelAttenCoeff(grooveIdx))
            : 1.0f;

        uint16_t swingDelay = 0;
        float blocksPerTick = (params_.engineArp1.clock == CLOCK_EXTERNAL)
            ? extBlocksPerTick_
            : ((ticksEveryNCalls_ > 2.0f) ? ticksEveryNCalls_ : 31.25f);
        float stepBlocks = (float) stepTicks * blocksPerTick;
        if (curSwung && swingNorm > 0.0f) {
            int d = (int) roundf(swingNorm * stepBlocks);
            swingDelay = (uint16_t) ((d < 1) ? 1 : d);
        }
        if (humanize > 0 && curPatternStep > 0) {
            int hd = (int) roundf(fabsf(noise[0]) * stepBlocks * ((float) humanize * 0.0030f));
            if (hd > 0) {
                int maxTotal = (int) stepBlocks - 4;
                int totalDelay = (int) swingDelay + hd;
                if (maxTotal > 1 && totalDelay > maxTotal) {
                    totalDelay = maxTotal;
                }
                swingDelay = (uint16_t) totalDelay;
            }
        }
        int humnVel = (humanize > 0) ? (int) roundf(noise[1] * ((float) humanize * 0.40f)) : 0;
        const uint8_t effectivePart = getEffectiveCptPart();
        uint8_t has_arpeggiator_note = (bitmask_ & pattern) ? 255 : 0;
        bool has_accent = (bitmask_ & accentMask) != 0;
        int stepMode = getArpeggiatorStepMode();

        uint8_t poolIndices[16];
        const int num_notes = buildCounterpointPool(note_stack_, poolIndices);
        const int direction = params_.engineArp1.direction;

        CounterpointTimbreState &myCpt = cptTimbreBus[timbreNumber_];
        myCpt.stepTriggeredThisTick = 1;
        myCpt.curPatternStep = (uint8_t) curPatternStep;

        if (num_notes == 0) {
            stopArpRatchet();
            if (event_scheduler_.size() == 0 && cptTimbreBus[timbreNumber_].activeNote != 0) {
                preenNoteOff((char) cptTimbreBus[timbreNumber_].activeNote, false);
            }
        }
        if (num_notes && has_arpeggiator_note) {
            if (swingDelay == 0) {
                stopArpRatchet();
            }
            arpPendingStepMode_ = (uint8_t) stepMode;
            if (voiceNumber_[0] >= 0) {
                voices_[voiceNumber_[0]]->arpStepSeqNoteOn();
            }
            if (ARPEGGIO_DIRECTION_CHORD != direction && ARPEGGIO_DIRECTION_STRIDE != direction && stepMode != ARP_STEP_MODE_CHORD) {
                StepArpeggio((uint8_t) num_notes);
                int step, transpose = 0;
                if (direction == ARPEGGIO_DIRECTION_CONVERGE) {
                    int pos = current_step_;
                    step = (pos & 1) ? ((num_notes - 1) - (pos >> 1)) : (pos >> 1);
                } else if (direction == ARPEGGIO_DIRECTION_DIVERGE) {
                    int revPos = (num_notes - 1) - current_step_;
                    step = (revPos & 1) ? ((num_notes - 1) - (revPos >> 1)) : (revPos >> 1);
                } else if (direction == ARPEGGIO_DIRECTION_THUMB && num_notes > 2) {
                    step = (current_step_ & 1) ? (1 + (current_step_ >> 1)) : 0;
                } else if (direction == ARPEGGIO_DIRECTION_PINKY && num_notes > 2) {
                    step = (current_step_ & 1) ? (num_notes - 1) : (current_step_ >> 1);
                } else if (direction == ARPEGGIO_DIRECTION_LEAP) {
                    if (num_notes > 2) {
                        int half = (num_notes + 1) >> 1;
                        step = (current_step_ & 1) ? (half + (current_step_ >> 1)) : (current_step_ >> 1);
                        if (step >= num_notes) {
                            step = num_notes - 1;
                        }
                    } else {
                        step = current_step_;
                    }
                } else if (current_direction_ > 0) {
                    step = start_step_ + current_step_;
                    if (step >= num_notes) {
                        step -= num_notes;
                        transpose = 12;
                    }
                } else {
                    step = (num_notes - 1) - (start_step_ + current_step_);
                    if (step < 0) {
                        step += num_notes;
                        transpose = -12;
                    }
                }
                if (step < 0) {
                    step = 0;
                } else if (step >= num_notes) {
                    step = num_notes - 1;
                }

                const NoteEntry &noteEntry = (ARPEGGIO_DIRECTION_PLAYED == direction && effectivePart == CPT_PART_ALL)
                    ? note_stack_.played_note(poolIndices[step])
                    : note_stack_.sorted_note(poolIndices[step]);

                uint8_t note = noteEntry.note;
                int effVel = (int) roundf((float) noteEntry.velocity * swingVelScale) + (has_accent ? 32 : 0) + humnVel;
                if (effVel < 1) {
                    effVel = 1;
                } else if (effVel > 127) {
                    effVel = 127;
                }
                uint8_t velocity = (uint8_t) effVel;
                int totalTranspose = 12 * (int) current_octave_;
                if (__canTranspose(direction)) {
                    totalTranspose += transpose;
                }
                if (stepMode == ARP_STEP_MODE_OCT_UP) {
                    totalTranspose += 12;
                } else if (stepMode == ARP_STEP_MODE_OCT_DOWN && (int) note + totalTranspose >= 12) {
                    totalTranspose -= 12;
                }
                int rawPitch = (int) note + totalTranspose;
                while (rawPitch > 127) {
                    rawPitch -= 12;
                }
                if (rawPitch < 0) {
                    rawPitch = 0;
                }
                note = (uint8_t) rawPitch;

                uint8_t compNote = 0;
                bool hasDyad = false;
                if (stepMode == ARP_STEP_MODE_DYAD) {
                    if (num_notes >= 2) {
                        int compIdx = (num_notes >= 3) ? ((step + 2) % num_notes) : ((step + 1) & 1);
                        const NoteEntry &compEntry = (ARPEGGIO_DIRECTION_PLAYED == direction && effectivePart == CPT_PART_ALL)
                            ? note_stack_.played_note(poolIndices[compIdx])
                            : note_stack_.sorted_note(poolIndices[compIdx]);
                        int cn = (int) compEntry.note + 12 * current_octave_ + (__canTranspose(direction) ? transpose : 0);
                        while (cn > 127) {
                            cn -= 12;
                        }
                        if (cn >= 0 && cn != note) {
                            compNote = (uint8_t) cn;
                            hasDyad = true;
                        }
                    } else {
                        compNote = (note <= 120) ? (uint8_t) (note + 7) : (uint8_t) (note - 5);
                        hasDyad = true;
                    }
                }

                if (swingDelay == 0) {
                    SendNote(note, velocity, stepMode);
                    if (hasDyad) {
                        uint8_t compVel = (uint8_t) ((velocity * 7) >> 3);
                        if (compVel < 1) {
                            compVel = 1;
                        }
                        SendNote(compNote, compVel, stepMode);
                    }
                } else {
                    arpStrumNotes_[0] = note;
                    arpStrumVelocities_[0] = velocity;
                    arpStrumDelay_[0] = swingDelay;
                    arpStrumCount_ = 1;
                    if (hasDyad) {
                        uint8_t compVel = (uint8_t) ((velocity * 7) >> 3);
                        if (compVel < 1) {
                            compVel = 1;
                        }
                        arpStrumNotes_[1] = compNote;
                        arpStrumVelocities_[1] = compVel;
                        arpStrumDelay_[1] = (uint16_t) (swingDelay + 3);
                        arpStrumCount_ = 2;
                    }
                }
            } else {
                uint32_t seed = (uint32_t) chord_cycle_ * 2654435761u;
                bool strumDown = has_accent
                    ? (((seed >> 8) % 5u) == 0u)
                    : (((seed >> 8) & 1u) != 0u);
                bool highPhase = (chord_cycle_ & 1u) != 0u;
                chord_cycle_++;

                int chordIndices[16];
                int chordCount = 0;
                if (direction == ARPEGGIO_DIRECTION_STRIDE && stepMode != ARP_STEP_MODE_CHORD) {
                    if (num_notes == 1) {
                        chordIndices[chordCount++] = poolIndices[0];
                    } else if (num_notes == 2) {
                        chordIndices[chordCount++] = poolIndices[highPhase ? 1 : 0];
                    } else if (!highPhase) {
                        chordIndices[chordCount++] = poolIndices[0];
                    } else if (num_notes == 3) {
                        chordIndices[chordCount++] = poolIndices[1];
                        chordIndices[chordCount++] = poolIndices[2];
                    } else {
                        chordIndices[chordCount++] = poolIndices[num_notes - 2];
                        chordIndices[chordCount++] = poolIndices[num_notes - 1];
                    }
                } else if (num_notes <= 3) {
                    for (int idx = 0; idx < num_notes && idx < 16; ++idx) {
                        chordIndices[chordCount++] = poolIndices[idx];
                    }
                } else if (!highPhase) {
                    chordIndices[chordCount++] = poolIndices[0];
                    chordIndices[chordCount++] = poolIndices[1];
                } else {
                    chordIndices[chordCount++] = poolIndices[num_notes - 2];
                    chordIndices[chordCount++] = poolIndices[num_notes - 1];
                }

                bool hasRatchet = (stepMode == ARP_STEP_MODE_RATCHET2)
                    || (voiceNumber_[0] >= 0 && fabsf(voices_[voiceNumber_[0]]->matrix.getDestination(MTX_DEST_ARP_RATCHET)) >= 0.04f);
                // Tight, organic strum spacing: ~2-5 ms per note (3-7 blocks at 1500 Hz).
                int baseStrum = has_accent
                    ? (3 + (int) ((seed >> 2) % 3u))
                    : (4 + (int) ((seed >> 2) % 4u));
                int accumDelay = 0;

                for (int idx = 0; idx < chordCount; ++idx) {
                    int i = strumDown ? chordIndices[chordCount - 1 - idx] : chordIndices[idx];
                    const NoteEntry &noteEntry = note_stack_.sorted_note(i);
                    int effVel = (int) roundf((float) noteEntry.velocity * swingVelScale) + (has_accent ? 24 : 0) + humnVel;
                    if (chordCount > 1) {
                        int velJitter = ((int) ((seed >> (12 + (idx & 3) * 3)) & 7u) - 3) * 2;
                        effVel += (-4 * idx) + velJitter;
                    }
                    if (effVel < 32) {
                        effVel = 32;
                    } else if (effVel > 127) {
                        effVel = 127;
                    }
                    uint8_t chordNote = noteEntry.note;
                    if (stepMode == ARP_STEP_MODE_OCT_UP && chordNote <= 115) {
                        chordNote += 12;
                    } else if (stepMode == ARP_STEP_MODE_OCT_DOWN && chordNote >= 12) {
                        chordNote -= 12;
                    }

                    if ((swingDelay == 0 && (idx == 0 || hasRatchet)) || arpStrumCount_ >= 8) {
                        SendNote(chordNote, (uint8_t) effVel, stepMode);
                    } else {
                        uint16_t noteDelay = swingDelay;
                        if (idx > 0 && !hasRatchet) {
                            int jitter = (int) ((seed >> (4 + (idx & 3) * 2)) & 3u) - 1;
                            int stepDelay = baseStrum + jitter;
                            if (stepDelay < 2) {
                                stepDelay = 2;
                            }
                            accumDelay += stepDelay;
                            if (accumDelay > 48) {
                                accumDelay = 48;
                            }
                            noteDelay = (uint16_t) (swingDelay + accumDelay);
                        }
                        arpStrumNotes_[arpStrumCount_] = chordNote;
                        arpStrumVelocities_[arpStrumCount_] = (uint8_t) effVel;
                        arpStrumDelay_[arpStrumCount_] = noteDelay;
                        arpStrumCount_++;
                    }
                }
            }
        }
        if (bitmask_ >= endMask) {
            bitmask_ = 1;
        } else {
            bitmask_ <<= 1;
        }
    }
}

void Timbre::StepArpeggio(uint8_t effectiveNumNotes) {

    if (current_octave_ == 127) {
        StartArpeggio();
        return;
    }

    int direction = params_.engineArp1.direction;
    uint8_t num_notes = (effectiveNumNotes > 0) ? effectiveNumNotes : note_stack_.size();
    if (num_notes == 0) {
        return;
    }
    if (direction == ARPEGGIO_DIRECTION_RANDOM) {
        uint8_t random_byte = *(uint8_t*) noise;
        current_octave_ = random_byte & 0xf;
        current_step_ = (random_byte & 0xf0) >> 4;
        while (current_octave_ >= params_.engineArp1.octave) {
            current_octave_ -= params_.engineArp1.octave;
        }
        while (current_step_ >= num_notes) {
            current_step_ -= num_notes;
        }
    } else if (direction == ARPEGGIO_DIRECTION_RANDOM_WALK) {
        uint8_t random_byte = *(uint8_t*) noise;
        int delta = (random_byte % 3) - 1;
        int totalNotes = num_notes * (int) params_.engineArp1.octave;
        int flat = (int) current_octave_ * num_notes + (int) current_step_ + delta;
        if (flat < 0) {
            flat = (totalNotes > 1) ? 1 : 0;
        } else if (flat >= totalNotes) {
            flat = (totalNotes > 1) ? (totalNotes - 2) : 0;
        }
        current_octave_ = (uint8_t) (flat / num_notes);
        current_step_ = (uint8_t) (flat % num_notes);
    } else if (direction == ARPEGGIO_DIRECTION_OCTAVE_NOTE) {
        if (++current_octave_ >= params_.engineArp1.octave) {
            current_octave_ = 0;
            ++current_step_;
        }
        if (current_step_ >= num_notes) {
            current_step_ = 0;
        }
    } else {
        // NOTE: We always count [0 - cycleLen) here; the actual handling of direction is in Tick()
        uint8_t cycleLen = num_notes;
        if ((direction == ARPEGGIO_DIRECTION_THUMB || direction == ARPEGGIO_DIRECTION_PINKY) && num_notes > 2) {
            cycleLen = (uint8_t) (2 * (num_notes - 1));
        }

        uint8_t trigger_change = 0;
        if (++current_step_ >= cycleLen) {
            current_step_ = 0;
            trigger_change = 1;
        }

        // special case the 'ROTATE' and 'SHIFT' modes, they might not change the octave until the cycle is through
        if (trigger_change && (direction >= ARPEGGIO_DIRECTION_ROTATE_UP && direction <= ARPEGGIO_DIRECTION_SHIFT_UP_DOWN)) {
            if (++start_step_ >= num_notes)
                start_step_ = 0;
            else
                trigger_change = 0;
        }

        if (trigger_change) {
            current_octave_ += current_direction_;
            if (current_octave_ >= params_.engineArp1.octave || current_octave_ < 0) {
                if (__canChangeDir(direction)) {
                    current_direction_ = -current_direction_;
                    StartArpeggio();
                    if (direction != ARPEGGIO_DIRECTION_UP_DOWN_INC && (num_notes > 1 || params_.engineArp1.octave > 1)) {
                        StepArpeggio(num_notes);
                    }
                } else {
                    StartArpeggio();
                }
            }
        }
    }
}

void Timbre::StartArpeggio() {

    current_step_ = 0;
    start_step_ = 0;
    if (current_direction_ == 1) {
        current_octave_ = 0;
    } else {
        current_octave_ = params_.engineArp1.octave - 1;
    }
}

void Timbre::resetArpPatternState() {
    bitmask_ = 1;
    chord_cycle_ = 0;
    arpStrumCount_ = 0;
    arpPendingStepMode_ = ARP_STEP_MODE_NORM;
    arpSwingStepScale_ = 1.0f;
    cptTimbreBus[timbreNumber_].curPatternStep = 0;
    for (int k = 0; k < numberOfVoices_; k++) {
        if (voiceNumber_[k] >= 0) {
            voices_[voiceNumber_[k]]->arpStartStepSeq();
        }
    }
    current_octave_ = 127;
    current_direction_ = __getDirection(params_.engineArp1.direction);
}

void Timbre::getChannelStackInfo(int *outRank, int *outTotal, int *outMasterTimbre) const {
    *outRank = 0;
    *outTotal = 1;
    *outMasterTimbre = timbreNumber_;
    if (unlikely(mixerState_ == 0)) {
        return;
    }
    uint8_t myChan = mixerState_->instrumentState_[timbreNumber_].midiChannel;
    int rank = 0;
    int total = 0;
    int master = -1;
    for (int t = 0; t < NUMBER_OF_TIMBRES; t++) {
        if (mixerState_->instrumentState_[t].numberOfVoices == 0) {
            continue;
        }
        uint8_t tChan = mixerState_->instrumentState_[t].midiChannel;
        if (tChan == myChan || tChan == 0 || myChan == 0) {
            if (master < 0) {
                master = t;
            }
            if (t == timbreNumber_) {
                rank = total;
            }
            total++;
        }
    }
    *outRank = rank;
    *outTotal = (total > 0) ? total : 1;
    *outMasterTimbre = (master >= 0) ? master : timbreNumber_;
}

uint8_t Timbre::getEffectiveCptPart() const {
    uint8_t part = getCptPartMode();
    if (part != CPT_PART_DIVISI) {
        return part;
    }
    int rank = 0, total = 1, master = timbreNumber_;
    getChannelStackInfo(&rank, &total, &master);
    if (total <= 1) {
        return CPT_PART_ALL;
    }
    if (total == 2) {
        return (rank == 0) ? CPT_PART_BASS : CPT_PART_INNER;
    }
    if (total == 3) {
        return (rank == 0) ? CPT_PART_BASS : ((rank == 1) ? CPT_PART_INNER : CPT_PART_SOPRANO);
    }
    if (rank == 0) {
        return CPT_PART_BASS;
    }
    if (rank == 1) {
        return CPT_PART_CANTUS;
    }
    if (rank == 2) {
        return CPT_PART_ALTO;
    }
    return CPT_PART_SOPRANO;
}

int Timbre::buildCounterpointPool(const NoteStack &stack, uint8_t *outIndices) const {
    int n = (int) stack.size();
    if (n <= 0) {
        return 0;
    }
    uint8_t rawPart = getCptPartMode();
    bool arpOn = params_.engineArp1.clock > 0.0f;
    int count = 0;

    if (rawPart == CPT_PART_DIVISI) {
        int rank = 0, total = 1, master = timbreNumber_;
        getChannelStackInfo(&rank, &total, &master);
        if (total <= 1) {
            for (int i = 0; i < n && count < 16; ++i) {
                outIndices[count++] = (uint8_t) i;
            }
            return count;
        }
        // Strict Divisi: timbre of rank r only enters when at least (r + 1) notes are held.
        if (n <= rank) {
            return 0;
        }
        if (rank < total - 1) {
            outIndices[count++] = (uint8_t) rank;
        } else {
            for (int i = rank; i < n && count < 16; ++i) {
                outIndices[count++] = (uint8_t) i;
            }
        }
        return count;
    }

    switch (rawPart) {
    case CPT_PART_BASS:
        outIndices[count++] = 0;
        break;
    case CPT_PART_CANTUS:
        // Strict Divisi: Tenor enters from the 2nd held note.
        if (n < 2) {
            return 0;
        }
        if (n <= 3) {
            outIndices[count++] = 1;
        } else if (arpOn) {
            outIndices[count++] = 1;
            outIndices[count++] = 2;
        } else {
            uint8_t prev = cptTimbreBus[timbreNumber_].activeNote;
            if (prev == 0) {
                prev = cptTimbreBus[timbreNumber_].prevNote;
            }
            if (prev > 0) {
                int d1 = (int) stack.sorted_note(1).note - (int) prev;
                int d2 = (int) stack.sorted_note(2).note - (int) prev;
                if (d1 < 0) d1 = -d1;
                if (d2 < 0) d2 = -d2;
                outIndices[count++] = (d2 < d1) ? 2 : 1;
            } else {
                outIndices[count++] = 1;
            }
        }
        break;
    case CPT_PART_ALTO:
        // Strict Divisi: Alto enters from the 3rd held note.
        if (n < 3) {
            return 0;
        }
        if (n == 3 || !arpOn) {
            outIndices[count++] = (uint8_t) (n - 2);
        } else {
            outIndices[count++] = (uint8_t) (n - 3);
            outIndices[count++] = (uint8_t) (n - 2);
        }
        break;
    case CPT_PART_SOPRANO: {
        int rank = 0, total = 1, master = timbreNumber_;
        getChannelStackInfo(&rank, &total, &master);
        // Strict Divisi: if stacked with other timbres on the channel, do not double a single bass note.
        if (total > 1 && n < 2) {
            return 0;
        }
        outIndices[count++] = (uint8_t) (n - 1);
        break;
    }
    case CPT_PART_INNER:
        // Strict Divisi: Inner voices only exist when at least 3 notes are held (between Bass and Soprano).
        if (n < 3) {
            return 0;
        }
        for (int i = 1; i <= n - 2 && count < 16; ++i) {
            outIndices[count++] = (uint8_t) i;
        }
        break;
    case CPT_PART_OUTER:
        outIndices[count++] = 0;
        if (n >= 2) {
            outIndices[count++] = (uint8_t) (n - 1);
        }
        break;
    case CPT_PART_UPPER:
        // Strict Divisi: everything except the lowest note (Bass).
        if (n < 2) {
            return 0;
        }
        for (int i = 1; i < n && count < 16; ++i) {
            outIndices[count++] = (uint8_t) i;
        }
        break;
    case CPT_PART_LOWER:
        // Everything except the highest note (Soprano) when >= 2 notes are held.
        if (n == 1) {
            outIndices[count++] = 0;
        } else {
            for (int i = 0; i <= n - 2 && count < 16; ++i) {
                outIndices[count++] = (uint8_t) i;
            }
        }
        break;
    case CPT_PART_EVEN:
        for (int i = 0; i < n && count < 16; i += 2) {
            outIndices[count++] = (uint8_t) i;
        }
        break;
    case CPT_PART_ODD:
        // Strict Divisi: interleaved odd indices enter from the 2nd held note.
        if (n < 2) {
            return 0;
        }
        for (int i = 1; i < n && count < 16; i += 2) {
            outIndices[count++] = (uint8_t) i;
        }
        break;
    case CPT_PART_LAST: {
        uint8_t recentNote = stack.most_recent_note().note;
        uint8_t foundIdx = (uint8_t) (n - 1);
        for (int i = 0; i < n; ++i) {
            if (stack.sorted_note(i).note == recentNote) {
                foundIdx = (uint8_t) i;
                break;
            }
        }
        outIndices[count++] = foundIdx;
        break;
    }
    default:
        for (int i = 0; i < n && count < 16; ++i) {
            outIndices[count++] = (uint8_t) i;
        }
        break;
    }
    return count;
}

void Timbre::reconcileCounterpointDirectNotes() {
    if (unlikely(numberOfVoices_ == 0)) {
        return;
    }
    int heldCount = (int) held_note_stack_.size();
    uint8_t poolIndices[16];
    int poolSize = (heldCount > 0) ? buildCounterpointPool(held_note_stack_, poolIndices) : 0;
    if (poolSize <= 0) {
        for (int k = 0; k < numberOfVoices_; k++) {
            int v = voiceNumber_[k];
            if (v >= 0 && voices_[v]->isPlaying() && !voices_[v]->isReleased() && !voices_[v]->isHoldedByPedal()) {
                preenNoteOff(voices_[v]->getNote(), false);
            }
        }
        return;
    }

    if (params_.engine1.playMode != PLAY_MODE_POLY || numberOfVoices_ == 1) {
        uint8_t prio = getVoicePriorityMode();
        int chosenPoolIdx = (prio == VOICE_PRIORITY_HIGH_KEY || prio == VOICE_PRIORITY_BASS_PROTECT)
            ? (poolSize - 1)
            : 0;
        const NoteEntry &entry = held_note_stack_.sorted_note(poolIndices[chosenPoolIdx]);
        uint8_t targetNote = entry.note;
        int v0 = voiceNumber_[0];
        if (v0 >= 0 && voices_[v0]->isPlaying() && !voices_[v0]->isReleased()
                && ((uint8_t) voices_[v0]->getNote() & 0x7f) == targetNote) {
            return;
        }
        preenNoteOn((char) targetNote, (char) entry.velocity, true);
        return;
    }

    uint8_t targetNotes[16];
    uint8_t targetVels[16];
    for (int i = 0; i < poolSize; ++i) {
        const NoteEntry &e = held_note_stack_.sorted_note(poolIndices[i]);
        targetNotes[i] = e.note;
        targetVels[i] = e.velocity;
    }

    for (int k = 0; k < numberOfVoices_; k++) {
        int v = voiceNumber_[k];
        if (v < 0 || !voices_[v]->isPlaying() || voices_[v]->isReleased() || voices_[v]->isHoldedByPedal()) {
            continue;
        }
        uint8_t playingNote = (uint8_t) voices_[v]->getNote() & 0x7f;
        bool inTarget = false;
        for (int i = 0; i < poolSize; ++i) {
            if (targetNotes[i] == playingNote) {
                inTarget = true;
                break;
            }
        }
        if (!inTarget) {
            preenNoteOff((char) playingNote, false);
        }
    }

    for (int i = 0; i < poolSize; ++i) {
        if (voice_allocator_.FindActive(targetNotes[i]) == 0xff) {
            preenNoteOn((char) targetNotes[i], (char) targetVels[i], true);
        }
    }
}

void Timbre::updateCounterpointMode() {
    if (params_.engineArp1.clock == 0.0f && held_note_stack_.size() > 0) {
        reconcileCounterpointDirectNotes();
    }
}

void Timbre::Start() {
    resetArpPatternState();
    recording_ = 0;
    running_ = 1;
    tick_ = getArpStepTicks() - 1;
}

void Timbre::arpeggiatorSetHoldPedal(uint8_t value) {
    if (ignore_note_off_messages_ && !value) {
        // Pedal was released, kill all pending arpeggios.
        note_stack_.Clear();
    }
    ignore_note_off_messages_ = value;
}

void Timbre::setLatchMode(uint8_t value) {
    // When disabling latch mode, clear the note stack.
    latch_ = value;
    if (value == 0) {
        note_stack_.Clear();
        recording_ = 0;
    }
}

void Timbre::setDirection(uint8_t value) {
    // When changing the arpeggio direction, reset the pattern.
    current_direction_ = __getDirection(value);
    StartArpeggio();
}

void Timbre::lfoValueChange(int currentRow, int encoder, float newValue) {
    for (int k = 0; k < numberOfVoices_; k++) {
        voices_[voiceNumber_[k]]->lfoValueChange(currentRow, encoder, newValue);
    }
}

void Timbre::updateMidiNoteScale(int scale) {

    int intBreakNote;
    int curveBefore;
    int curveAfter;
    if (scale == 0) {
        intBreakNote = params_.midiNote1Curve.breakNote;
        curveBefore = params_.midiNote1Curve.curveBefore;
        curveAfter = params_.midiNote1Curve.curveAfter;
    } else {
        intBreakNote = params_.midiNote2Curve.breakNote;
        curveBefore = params_.midiNote2Curve.curveBefore;
        curveAfter = params_.midiNote2Curve.curveAfter;
    }
    float floatBreakNote = intBreakNote;
    float multiplier = 1.0f;

    switch (curveBefore) {
        case MIDI_NOTE_CURVE_FLAT:
            for (int n = 0; n < intBreakNote; n++) {
                midiNoteScale[scale][timbreNumber_][n] = 0;
            }
            break;
        case MIDI_NOTE_CURVE_M_LINEAR:
            multiplier = -1.0f;
            goto linearBefore;
        case MIDI_NOTE_CURVE_M_LINEAR2:
            multiplier = -8.0f;
            goto linearBefore;
        case MIDI_NOTE_CURVE_LINEAR2:
            multiplier = 8.0f;
            goto linearBefore;
        case MIDI_NOTE_CURVE_LINEAR:
linearBefore: for (int n = 0; n < intBreakNote; n++) {
                float fn = (floatBreakNote - n);
                midiNoteScale[scale][timbreNumber_][n] = fn * INV127 * multiplier;
            }
            break;
        case MIDI_NOTE_CURVE_M_EXP:
            multiplier = -1.0f;
        case MIDI_NOTE_CURVE_EXP:
            if (floatBreakNote == 0) {
                floatBreakNote = 1;
            }
            for (int n = 0; n < intBreakNote; n++) {
                float fn = (floatBreakNote - n);
                fn = fn * fn / floatBreakNote;
                midiNoteScale[scale][timbreNumber_][n] = fn * INV16 * multiplier;
            }
            break;
    }

    // BREAK NOTE = 0;
    midiNoteScale[scale][timbreNumber_][intBreakNote] = 0;

    switch (curveAfter) {
        case MIDI_NOTE_CURVE_FLAT:
            for (int n = intBreakNote + 1; n < 128; n++) {
                midiNoteScale[scale][timbreNumber_][n] = 0;
            }
            break;
        case MIDI_NOTE_CURVE_M_LINEAR:
            multiplier = -1.0f;
            goto linearAfter;
        case MIDI_NOTE_CURVE_M_LINEAR2:
            multiplier = -8.0f;
            goto linearAfter;
        case MIDI_NOTE_CURVE_LINEAR2:
            multiplier = 8.0f;
            goto linearAfter;
        case MIDI_NOTE_CURVE_LINEAR:
linearAfter: for (int n = intBreakNote + 1; n < 128; n++) {
                float fn = n - floatBreakNote;
                midiNoteScale[scale][timbreNumber_][n] = fn * INV127 * multiplier;
            }
            break;
        case MIDI_NOTE_CURVE_M_EXP:
            multiplier = -1.0f;
        case MIDI_NOTE_CURVE_EXP:
            if (floatBreakNote == 0) {
                floatBreakNote = 1;
            }
            for (int n = intBreakNote + 1; n < 128; n++) {
                float fn = n - floatBreakNote;
                fn = fn * fn / floatBreakNote;
                midiNoteScale[scale][timbreNumber_][n] = fn * INV16 * multiplier;
            }
            break;
    }
    /*
     lcd.setCursor(0,0);
     lcd.print((int)(midiNoteScale[timbreNumber][25] * 127.0f));
     lcd.print(" ");
     lcd.setCursor(10,0);
     lcd.print((int)(midiNoteScale[timbreNumber][intBreakNote - 5] * 127.0f));
     lcd.print(" ");
     lcd.setCursor(0,1);
     lcd.print((int)(midiNoteScale[timbreNumber][intBreakNote + 5] * 127.0f));
     lcd.print(" ");
     lcd.setCursor(10,1);
     lcd.print((int)(midiNoteScale[timbreNumber][102] * 127.0f));
     lcd.print(" ");
     */

}

void Timbre::midiClockContinue(int songPosition) {

    for (int k = 0; k < numberOfVoices_; k++) {
        voices_[voiceNumber_[k]]->midiClockContinue(songPosition);
    }

    recomputeNext_ = ((songPosition & 0x1) == 0);
    OnMidiContinue();
}

void Timbre::midiClockStart() {

    for (int k = 0; k < numberOfVoices_; k++) {
        voices_[voiceNumber_[k]]->midiClockStart();
    }

    recomputeNext_ = true;
    OnMidiStart();
}

void Timbre::midiClockSongPositionStep(int songPosition) {

    for (int k = 0; k < numberOfVoices_; k++) {
        voices_[voiceNumber_[k]]->midiClockSongPositionStep(songPosition, recomputeNext_);
    }

    if ((songPosition & 0x1) == 0) {
        recomputeNext_ = true;
    }
}

void Timbre::resetMatrixDestination(float oldValue) {
    for (int k = 0; k < numberOfVoices_; k++) {
        voices_[voiceNumber_[k]]->matrix.resetDestination(oldValue);
    }
}

void Timbre::setMatrixSource(enum SourceEnum source, float newValue) {
    for (int k = 0; k < numberOfVoices_; k++) {
        voices_[voiceNumber_[k]]->matrix.setSource(source, newValue);
    }
}

void Timbre::setMatrixSourceMPE(uint8_t channel, enum SourceEnum source, float newValue) {
    int voiceToUse = voiceNumber_[channel -1];

    if (unlikely(voiceToUse == -1)) {
        return;
    }

	voices_[voiceToUse]->matrix.setSource(source, newValue);
}


void Timbre::setMatrixPolyAfterTouch(uint8_t note, float newValue) {
    for (int k = 0; k < numberOfVoices_; k++) {
        int voiceIndex = voiceNumber_[k];
        if (unlikely(voices_[voiceIndex]->isPlaying() && voices_[voiceIndex]->getNote() == note)) {
            voices_[voiceNumber_[k]]->matrix.setSource(MATRIX_SOURCE_POLYPHONIC_AFTERTOUCH, newValue);
            break;
        }
    }
}



void Timbre::verifyLfoUsed(int encoder, float oldValue, float newValue) {
    // No need to recompute
    if (numberOfVoices_ == 0.0f) {
        return;
    }
    // We used it and still use it
    if (encoder == ENCODER_MATRIX_MUL && oldValue != 0.0f && newValue != 0.0f) {
        return;
    }

    for (int lfo = 0; lfo < NUMBER_OF_LFO; lfo++) {
        lfoUSed_[lfo] = 0;
    }

    MatrixRowParams *matrixRows = &params_.matrixRowState1;

    for (int r = 0; r < MATRIX_SIZE; r++) {
        if ((matrixRows[r].source >= MATRIX_SOURCE_LFO1 && matrixRows[r].source <= MATRIX_SOURCE_LFOSEQ2) && matrixRows[r].mul != 0.0f
            && (matrixRows[r].dest1 != 0.0f || matrixRows[r].dest2 != 0.0f)) {
            lfoUSed_[(int) matrixRows[r].source - MATRIX_SOURCE_LFO1]++;
        }

        // Check if we have a Mtx* that would require LFO even if mul is 0
        // http://ixox.fr/forum/index.php?topic=69220.0
        if (matrixRows[r].mul != 0.0f && matrixRows[r].source != 0.0f) {
            if (matrixRows[r].dest1 >= MTX1_MUL && matrixRows[r].dest1 <= MTX4_MUL) {
                int index = matrixRows[r].dest1 - MTX1_MUL;
                if (matrixRows[index].source >= MATRIX_SOURCE_LFO1 && matrixRows[index].source <= MATRIX_SOURCE_LFOSEQ2) {
                    lfoUSed_[(int) matrixRows[index].source - MATRIX_SOURCE_LFO1]++;
                }
            }
            // same test for dest2
            if (matrixRows[r].dest2 >= MTX1_MUL && matrixRows[r].dest2 <= MTX4_MUL) {
                int index = matrixRows[r].dest2 - MTX1_MUL;
                if (matrixRows[index].source >= MATRIX_SOURCE_LFO1 && matrixRows[index].source <= MATRIX_SOURCE_LFOSEQ2) {
                    lfoUSed_[(int) matrixRows[index].source - MATRIX_SOURCE_LFO1]++;
                }
            }
        }
    }

    // Keep LFO 1..3 running when referenced as a step trig/ratchet gate in an active StepSeq.
    for (int s = 0; s < 2; s++) {
        if (lfoUSed_[5 + s] > 0) {
            const StepSequencerSteps *seqSteps = (s == 0) ? &params_.lfoSteps1 : &params_.lfoSteps2;
            for (int k = 0; k < 16; k++) {
                int mode = ((uint8_t) seqSteps->steps[k] >> 4) & 0x0F;
                if (mode >= STEPSEQ_MODE_LFO1 && mode <= STEPSEQ_MODE_LFO3) {
                    lfoUSed_[mode - STEPSEQ_MODE_LFO1]++;
                }
            }
        }
    }

    /*
     lcd.setCursor(11, 1);
     lcd.print('>');
     for (int lfo=0; lfo < NUMBER_OF_LFO; lfo++) {
     lcd.print((int)lfoUSed[lfo]);
     }
     lcd.print('<');
     */

}
