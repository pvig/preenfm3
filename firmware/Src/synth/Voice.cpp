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

#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wsequence-point"
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#include "Voice.h"
#include "Timbre.h"

#if defined(__GNUC__) && !defined(DEBUG)
#pragma GCC push_options
#pragma GCC optimize("Ofast", "fast-math")
#endif

float Voice::glidePhaseInc[13];
float Voice::mpeBitchBend[6];
float Voice::envEndScratch_[6] __attribute__((section(".ram_d1")));
float Voice::envIncScratch_[6] __attribute__((section(".ram_d1")));

static float voiceDecimHistory[MAX_NUMBER_OF_VOICES][2][6] __attribute__((section(".ram_d1")));
static float osFirstHalfScratch[BLOCK_SIZE * 2] __attribute__((section(".ram_d1")));
static uint8_t nextVoiceDecimIndex = 0;

Voice::Voice(void)
{

    if (glidePhaseInc[0] != .2f)
    {
        float tmp[] = {2.0f, 5.0f, 9.0f, 15.0f, 22.0f, 35.0f, 50.0f, 90.0f, 140.0f, 200.0f, 500.0f, 1200.0f, 2700.0f};
        for (int k = 0; k < 13; k++)
        {
            glidePhaseInc[k] = 1.0f / tmp[k];
        }
    }

    // We multiply by 10 the number of semitone to have the direct index in exp2_harm
    if (mpeBitchBend[1] != 480.0f)
    {
        float tmp[] = {0.0f, 480.0f, 360.0f, 240.0f, 120.0, 0.0f};
        for (int k = 0; k < 6; k++)
        {
            mpeBitchBend[k] = tmp[k];
        }
    }

    freqAi = freqAo = 0.0f;
    freqBi = freqBo = 0.0f;
    freqHarm = 1.0f;
    baseFreqHarm = 1.0f;
    targetFreqHarm = 0.0f;
    noteAddScaleRatio = 1.0f;
    lastNoteAdd = 0;
    lastDecimationMode_ = -1;
    oversample2x_ = 0;
    fxNoteTrigger_ = 0;
    decimHistory_ = voiceDecimHistory[nextVoiceDecimIndex % MAX_NUMBER_OF_VOICES];
    nextVoiceDecimIndex++;
    for (int k = 0; k < 6; k++)
    {
        decimHistory_[0][k] = 0.0f;
        decimHistory_[1][k] = 0.0f;
    }

    // Init FX variables
    v0L = v1L = v2L = v3L = v4L = v5L = v6L = v7L = v8L = v0R = v1R = v2R = v3R = v4R = v5R = v6R = v7R = v8R = 0.0f;
}

Voice::~Voice(void)
{
}

void Voice::init()
{
    this->currentTimbre = 0;
    this->playing = false;
    this->isFullOfZero = false;
    this->newNotePending = false;
    this->note = 0;
    this->midiVelocity = 0;
    this->holdedByPedal = false;
    this->newNotePlayed = false;
    this->noteAlreadyFinished = 0;
    this->lastDecimationMode_ = -1;
    this->oversample2x_ = 0;
    this->nextMainFrequency = 0.0f;
    for (int k = 0; k < 6; k++)
    {
        this->decimHistory_[0][k] = 0.0f;
        this->decimHistory_[1][k] = 0.0f;
    }
}

void Voice::glideToNote(float newNoteFrequency)
{
    // Must glide...
    this->gliding = true;
    this->glidePhase = 0.0f;

    currentTimbre->osc1_.glideToNote(&oscState1_, newNoteFrequency);
    currentTimbre->osc2_.glideToNote(&oscState2_, newNoteFrequency);
    currentTimbre->osc3_.glideToNote(&oscState3_, newNoteFrequency);
    currentTimbre->osc4_.glideToNote(&oscState4_, newNoteFrequency);
    currentTimbre->osc5_.glideToNote(&oscState5_, newNoteFrequency);
    currentTimbre->osc6_.glideToNote(&oscState6_, newNoteFrequency);
}

void Voice::noteOnWithoutPop(short newNote, float newNoteFrequency, short velocity, float phase)
{
    // We can glide in mono and unison
    if (!this->released && !this->newNotePending && currentTimbre->params_.engine1.playMode != PLAY_MODE_POLY && currentTimbre->params_.engine2.glideType != GLIDE_TYPE_OFF)
    {
        // Keep phase continuity on overlap glide: only pitch glides to the new note.
        glideToNote(newNoteFrequency);
        this->note = newNote;
        this->midiVelocity = velocity;
        int zeroVelo = (16 - currentTimbre->params_.engine1.velocity) * 8;
        int newVelocity = zeroVelo + ((velocity * (128 - zeroVelo)) >> 7);
        this->velocity = newVelocity * .0078125f;
        this->fxNoteTrigger_ = 1;
        this->noteFrequency = newNoteFrequency;
        this->nextMainFrequency = newNoteFrequency;
        this->holdedByPedal = false;
        this->noteAlreadyFinished = 0;
        this->lastNoteAdd = 0x7FFF;
    }
    else
    {
        // update note now so that the noteOff is triggered by the new note
        this->note = newNote;
        this->midiVelocity = velocity;
        this->noteFrequency = newNoteFrequency;
        // Quick dead !
        bool alreadyPending = this->newNotePending;
        this->newNotePending = true;
        this->pendingNoteVelocity = velocity;
        this->pendingNote = newNote;
        this->pendingNoteFrequency = newNoteFrequency;
        this->phase_ = phase;

        // Not release anymore... not available for new notes...
        this->released = false;
        this->holdedByPedal = false;
        this->noteAlreadyFinished = 0;

        if (!alreadyPending)
        {
            this->currentTimbre->env1_.noteOffQuick(&envState1_);
            this->currentTimbre->env2_.noteOffQuick(&envState2_);
            this->currentTimbre->env3_.noteOffQuick(&envState3_);
            this->currentTimbre->env4_.noteOffQuick(&envState4_);
            this->currentTimbre->env5_.noteOffQuick(&envState5_);
            this->currentTimbre->env6_.noteOffQuick(&envState6_);
        }
        // We don't want noteOnAfterMatrixCompute to restart the env or => Hanging note !!!
        this->newNotePlayed = false;
    }
}

void Voice::glide()
{
    this->glidePhase += glidePhaseInc[(int)(currentTimbre->params_.engine1.glideSpeed + .05f)];
    if (glidePhase < 1.0f)
    {

        currentTimbre->osc1_.glideStep(&oscState1_, this->glidePhase);
        currentTimbre->osc2_.glideStep(&oscState2_, this->glidePhase);
        currentTimbre->osc3_.glideStep(&oscState3_, this->glidePhase);
        currentTimbre->osc4_.glideStep(&oscState4_, this->glidePhase);
        currentTimbre->osc5_.glideStep(&oscState5_, this->glidePhase);
        currentTimbre->osc6_.glideStep(&oscState6_, this->glidePhase);
    }
    else
    {
        // last with phase set to 1 to have exact frequencry
        currentTimbre->osc1_.glideStep(&oscState1_, 1);
        currentTimbre->osc2_.glideStep(&oscState2_, 1);
        currentTimbre->osc3_.glideStep(&oscState3_, 1);
        currentTimbre->osc4_.glideStep(&oscState4_, 1);
        currentTimbre->osc5_.glideStep(&oscState5_, 1);
        currentTimbre->osc6_.glideStep(&oscState6_, 1);
        this->gliding = false;
    }
}

// This is the frequency used for oscilloscope sync
// It's based on oscState1_ characteristics

float Voice::getNoteRealFrequencyEstimation(float newNoteFrequency)
{
    return currentTimbre->osc1_.getNoteRealFrequencyEstimation(&oscState1_, newNoteFrequency);
}

static inline float phaseDegreeToNormalized(float phaseDegree)
{
    // Convert UI phase in degrees to wavetable phase in [0..1].
    if (phaseDegree < 0.0f)
    {
        phaseDegree = 0.0f;
    }
    else if (phaseDegree > 360.0f)
    {
        phaseDegree = 360.0f;
    }
    return phaseDegree * (1.0f / 360.0f);
}

static inline float phaseDegreeToNormalizedWrapped(float phaseDegree)
{
    while (phaseDegree < 0.0f)
    {
        phaseDegree += 360.0f;
    }
    while (phaseDegree > 360.0f)
    {
        phaseDegree -= 360.0f;
    }
    return phaseDegree * (1.0f / 360.0f);
}

void Voice::applyOperatorStartPhases(float mainFrequency)
{
    float oscPhase1 = phaseDegreeToNormalizedWrapped(currentTimbre->params_.phaseOp1.phase + matrix.getDestination(OSC1_PHASE));
    float oscPhase2 = phaseDegreeToNormalizedWrapped(currentTimbre->params_.phaseOp2.phase + matrix.getDestination(OSC2_PHASE));
    float oscPhase3 = phaseDegreeToNormalizedWrapped(currentTimbre->params_.phaseOp3.phase + matrix.getDestination(OSC3_PHASE));
    float oscPhase4 = phaseDegreeToNormalizedWrapped(currentTimbre->params_.phaseOp4.phase + matrix.getDestination(OSC4_PHASE));
    float oscPhase5 = phaseDegreeToNormalizedWrapped(currentTimbre->params_.phaseOp5.phase + matrix.getDestination(OSC5_PHASE));
    float oscPhase6 = phaseDegreeToNormalizedWrapped(currentTimbre->params_.phaseOp6.phase + matrix.getDestination(OSC6_PHASE));

    // Preserve legacy unison behavior for negative detune.
    if (unlikely(currentTimbre->params_.engine2.unisonDetune < 0.0f))
    {
        oscPhase1 = 0.25f;
        oscPhase2 = 0.25f;
        oscPhase3 = 0.25f;
        oscPhase4 = 0.25f;
        oscPhase5 = 0.25f;
        oscPhase6 = 0.25f;
    }

    currentTimbre->osc1_.newNote(&oscState1_, mainFrequency, oscPhase1);
    currentTimbre->osc2_.newNote(&oscState2_, mainFrequency, oscPhase2);
    currentTimbre->osc3_.newNote(&oscState3_, mainFrequency, oscPhase3);
    currentTimbre->osc4_.newNote(&oscState4_, mainFrequency, oscPhase4);
    currentTimbre->osc5_.newNote(&oscState5_, mainFrequency, oscPhase5);
    currentTimbre->osc6_.newNote(&oscState6_, mainFrequency, oscPhase6);
    this->lastDecimationMode_ = -1;
}

void Voice::noteOn(short newNote, float newNoteFrequency, short velocity, float phase)
{

    (void)phase;

    bool retriggerLfo = true;

    // On noteOn we can only glide in mono and unison with glideType ALWAYS
    if (unlikely(currentTimbre->params_.engine1.playMode != PLAY_MODE_POLY && currentTimbre->params_.engine2.glideType == GLIDE_TYPE_ALWAYS))
    {
        bool wasPlaying = this->playing;
        if (unlikely(this->nextMainFrequency != 0.0f))
        {
            this->noteFrequency = this->nextMainFrequency;
        }
        else
        {
            this->noteFrequency = newNoteFrequency;
        }
        // Keep phase continuity in glide-always mode: glide changes pitch only.
        glideToNote(newNoteFrequency);
        this->note = newNote;
        this->nextMainFrequency = newNoteFrequency;
        // In glide mode, overlapping notes are slides and should not retrigger one-shot LFO.
        retriggerLfo = !wasPlaying;
    }
    else
    {
        this->note = newNote;
        this->noteFrequency = newNoteFrequency;
        // Absolute per-operator start phase at note-on.
        applyOperatorStartPhases(newNoteFrequency);
    }

    this->midiVelocity = velocity;
    this->fxNoteTrigger_ = 1;
    this->released = false;
    this->playing = true;
    this->currentTimbre->newNotePlayed = true;
    this->isFullOfZero = false;
    this->pendingNote = 0;
    this->newNotePending = false;
    this->holdedByPedal = false;
    this->noteAlreadyFinished = 0;
    this->lastNoteAdd = 0x7FFF;
    this->fdbLastValue[0] = 0.0f;
    this->fdbLastValue[1] = 0.0f;
    this->fdbLastValue[2] = 0.0f;
    for (int k = 0; k < 6; k++)
    {
        this->decimHistory_[0][k] = 0.0f;
        this->decimHistory_[1][k] = 0.0f;
    }

    this->velIm1 = currentTimbre->params_.engineIm1.modulationIndexVelo1 * (float)velocity * .0078125f;
    this->velIm2 = currentTimbre->params_.engineIm1.modulationIndexVelo2 * (float)velocity * .0078125f;
    this->velIm3 = currentTimbre->params_.engineIm2.modulationIndexVelo3 * (float)velocity * .0078125f;
    this->velIm4 = currentTimbre->params_.engineIm2.modulationIndexVelo4 * (float)velocity * .0078125f;
    this->velIm5 = currentTimbre->params_.engineIm3.modulationIndexVelo5 * (float)velocity * .0078125f;
    this->velIm6 = currentTimbre->params_.engineIm3.modulationIndexVelo6 * (float)velocity * .0078125f;

    int zeroVelo = (16 - currentTimbre->params_.engine1.velocity) * 8;
    int newVelocity = zeroVelo + ((velocity * (128 - zeroVelo)) >> 7);
    this->velocity = newVelocity * .0078125f; // divide by 127

    // Tell nextBlock() to init Env...
    this->newNotePlayed = true;

    if (retriggerLfo)
    {
        lfoNoteOn();
    }
}

void Voice::endNoteOrBeginNextOne()
{
    if (this->newNotePending)
    {
        // pendingNote can be 255 : see Voice:noteOff
        if (pendingNote <= 127)
        {
            noteOn(pendingNote, pendingNoteFrequency, pendingNoteVelocity, phase_);
        }
        else
        {
            noteOn(pendingNote - 128, pendingNoteFrequency, pendingNoteVelocity, phase_);
            this->noteAlreadyFinished = 1;
        }
    }
    else
    {
        this->playing = false;
        this->released = false;
    }
    this->freqAi = 0.0f;
    this->freqAo = 0.0f;
    this->freqBi = 0.0f;
    this->freqBo = 0.0f;
}

void Voice::noteOff()
{

    if (unlikely(!this->playing))
    {
        return;
    }
    if (likely(!this->newNotePending))
    {
        if (this->newNotePlayed)
        {
            // Note hasn't played
            // Let's give it a chance to create one block then noteOff
            this->noteAlreadyFinished = 1;
            return;
        }
        this->released = true;
        this->gliding = false;
        this->holdedByPedal = false;

        currentTimbre->env1_.noteOff(&envState1_, &matrix);
        currentTimbre->env2_.noteOff(&envState2_, &matrix);
        currentTimbre->env3_.noteOff(&envState3_, &matrix);
        currentTimbre->env4_.noteOff(&envState4_, &matrix);
        currentTimbre->env5_.noteOff(&envState5_, &matrix);
        currentTimbre->env6_.noteOff(&envState6_, &matrix);

        lfoNoteOff();
    }
    else
    {
        // We receive a note off before the note actually started, let's flag it
        if (this->pendingNote <= 127)
        {
            this->pendingNote += 128;
        }
    }
}

// Called when loading new preset
void Voice::noteOffQuick()
{

    if (unlikely(!this->playing))
    {
        return;
    }

    if (this->newNotePlayed)
    {
        // Note hasn't played
        killNow();
        return;
    }
    this->released = true;
    this->gliding = false;
    this->holdedByPedal = false;

    currentTimbre->env1_.noteOffQuick(&envState1_);
    currentTimbre->env2_.noteOffQuick(&envState2_);
    currentTimbre->env3_.noteOffQuick(&envState3_);
    currentTimbre->env4_.noteOffQuick(&envState4_);
    currentTimbre->env5_.noteOffQuick(&envState5_);
    currentTimbre->env6_.noteOffQuick(&envState6_);

    lfoNoteOff();
}

void Voice::killNow()
{
    this->newNotePlayed = false;
    this->playing = false;
    this->newNotePending = false;
    this->pendingNote = 0;
    this->gliding = false;
    this->released = false;
    this->holdedByPedal = false;
    this->noteAlreadyFinished = 0;
    this->env1ValueMem = 0;
    this->env2ValueMem = 0;
    this->env3ValueMem = 0;
    this->env4ValueMem = 0;
    this->env5ValueMem = 0;
    this->env6ValueMem = 0;
}

void Voice::emptyBuffer()
{
    if (!isFullOfZero)
    {
        isFullOfZero = true;
        float *sample = sampleBlock;

        for (int k = 0; k < BLOCK_SIZE; k++)
        {
            *sample++ = 0.0f;
            *sample++ = 0.0f;
        }
    }
}

void Voice::nextBlock()
{
    // After matrix
    if (unlikely(this->newNotePlayed))
    {
        this->newNotePlayed = false;
        this->env1ValueMem = currentTimbre->env1_.noteOnAfterMatrixCompute(&envState1_, &matrix);
        this->env2ValueMem = currentTimbre->env2_.noteOnAfterMatrixCompute(&envState2_, &matrix);
        this->env3ValueMem = currentTimbre->env3_.noteOnAfterMatrixCompute(&envState3_, &matrix);
        this->env4ValueMem = currentTimbre->env4_.noteOnAfterMatrixCompute(&envState4_, &matrix);
        this->env5ValueMem = currentTimbre->env5_.noteOnAfterMatrixCompute(&envState5_, &matrix);
        this->env6ValueMem = currentTimbre->env6_.noteOnAfterMatrixCompute(&envState6_, &matrix);
    }

    updateAllMixOscsAndPans();

    // * 20 so that we have a full tone for full pitchbend
    bool updateFreqHarm = false;
    float newTargetFreqHarm = matrix.getDestination(ALL_OSC_FREQ_HARM) * 20.0f;
    if (unlikely(currentTimbre->getMPESetting() > 0))
    {
        newTargetFreqHarm += (matrix.getSource(MATRIX_SOURCE_PITCHBEND_MPE) * mpeBitchBend[currentTimbre->getMPESetting()]);
    }
    if (unlikely(newTargetFreqHarm != targetFreqHarm))
    {
        targetFreqHarm = newTargetFreqHarm;

        float findex = 512.0f + targetFreqHarm;
        int index = findex;
        // Max = 1024
        index &= 0x3ff;
        float fp = findex - index;
        baseFreqHarm = (exp2_harm[index] * (1.0f - fp) + exp2_harm[index + 1] * fp);
        updateFreqHarm = true;
    }

    int noteAdd = (int)roundf(matrix.getDestination(MTX_DEST_NOTE_ADD) * 12.0f);
    if (unlikely(noteAdd != lastNoteAdd))
    {
        lastNoteAdd = noteAdd;
        if (noteAdd == 0)
        {
            noteAddScaleRatio = 1.0f;
        }
        else
        {
            int baseNote = (int)((uint8_t)this->note & 0x7F);
            int effectiveNote = baseNote + noteAdd;
            if (effectiveNote < 0)
            {
                effectiveNote = 0;
            }
            else if (effectiveNote > 127)
            {
                effectiveNote = 127;
            }
            float baseScaleFreq = currentTimbre->getScaleFrequency(baseNote);
            float targetScaleFreq = currentTimbre->getScaleFrequency(effectiveNote);
            if (baseScaleFreq > 0.0f && targetScaleFreq > 0.0f)
            {
                noteAddScaleRatio = targetScaleFreq / baseScaleFreq;
            }
        }
        updateFreqHarm = true;
    }

    if (unlikely(updateFreqHarm))
    {
        freqHarm = baseFreqHarm * noteAddScaleRatio;
    }

    if (unlikely(oversample2x_))
    {
        runAlgoBlock(1, 0.015625f, 0.997625f, false);
        for (int k = 0; k < BLOCK_SIZE * 2; k++)
        {
            osFirstHalfScratch[k] = sampleBlock[k];
        }
        runAlgoBlock(2, 0.015625f, 0.997625f, true);
        decimateHalfBand2x(osFirstHalfScratch, sampleBlock);
    }
    else
    {
        runAlgoBlock(0, 0.03125f, 0.99525f, true);
    }

    if (unlikely(this->noteAlreadyFinished > 0))
    {
        if (noteAlreadyFinished == 2)
        {
            this->noteAlreadyFinished = 0;
            noteOff();
        }
        else
        {
            noteAlreadyFinished++;
        }
    }
}

#if defined(__GNUC__)
__attribute__((optimize("Ofast", "fast-math")))
#endif
void Voice::decimateHalfBand2x(const float *firstHalf, float *secondHalfAndOut)
{
    // 7-tap Kaiser half-band FIR (M=3, beta=3.2): 3 MACs/sample/channel, zero scratch copy.
    static constexpr float h0 = 0.497951501f;
    static constexpr float h1 = 0.284185170f;
    static constexpr float h3 = -0.033160920f;

    const float *a = firstHalf;
    float *b = secondHalfAndOut;
    float *hL = decimHistory_[0];
    float *hR = decimHistory_[1];

    float nextL[6];
    float nextR[6];
    for (int k = 0; k < 6; k++)
    {
        nextL[k] = b[52 + k * 2];
        nextR[k] = b[53 + k * 2];
    }

    for (int k = 31; k >= 19; k--)
    {
        int i = 2 * k - 35;
        b[2 * k] = h0 * b[2 * i] + h1 * (b[2 * (i - 1)] + b[2 * (i + 1)]) + h3 * (b[2 * (i - 3)] + b[2 * (i + 3)]);
        b[2 * k + 1] = h0 * b[2 * i + 1] + h1 * (b[2 * (i - 1) + 1] + b[2 * (i + 1) + 1]) + h3 * (b[2 * (i - 3) + 1] + b[2 * (i + 3) + 1]);
    }

    b[36] = h0 * b[2] + h1 * (b[0] + b[4]) + h3 * (a[60] + b[8]);
    b[37] = h0 * b[3] + h1 * (b[1] + b[5]) + h3 * (a[61] + b[9]);

    b[34] = h0 * a[62] + h1 * (a[60] + b[0]) + h3 * (a[56] + b[4]);
    b[35] = h0 * a[63] + h1 * (a[61] + b[1]) + h3 * (a[57] + b[5]);

    b[32] = h0 * a[58] + h1 * (a[56] + a[60]) + h3 * (a[52] + b[0]);
    b[33] = h0 * a[59] + h1 * (a[57] + a[61]) + h3 * (a[53] + b[1]);

    for (int k = 15; k >= 3; k--)
    {
        int c = 2 * k - 3;
        b[2 * k] = h0 * a[2 * c] + h1 * (a[2 * (c - 1)] + a[2 * (c + 1)]) + h3 * (a[2 * (c - 3)] + a[2 * (c + 3)]);
        b[2 * k + 1] = h0 * a[2 * c + 1] + h1 * (a[2 * (c - 1) + 1] + a[2 * (c + 1) + 1]) + h3 * (a[2 * (c - 3) + 1] + a[2 * (c + 3) + 1]);
    }

    b[4] = h0 * a[2] + h1 * (a[0] + a[4]) + h3 * (hL[4] + a[8]);
    b[5] = h0 * a[3] + h1 * (a[1] + a[5]) + h3 * (hR[4] + a[9]);

    b[2] = h0 * hL[5] + h1 * (hL[4] + a[0]) + h3 * (hL[2] + a[4]);
    b[3] = h0 * hR[5] + h1 * (hR[4] + a[1]) + h3 * (hR[2] + a[5]);

    b[0] = h0 * hL[3] + h1 * (hL[2] + hL[4]) + h3 * (hL[0] + a[0]);
    b[1] = h0 * hR[3] + h1 * (hR[2] + hR[4]) + h3 * (hR[0] + a[1]);

    for (int k = 0; k < 6; k++)
    {
        hL[k] = nextL[k];
        hR[k] = nextR[k];
    }
}

void Voice::runAlgoBlock(int osMode, float invStep, float dcPole, bool allowEndNote)
{
    float env1Value;
    float env2Value;
    float env3Value;
    float env4Value;
    float env5Value;
    float env6Value;

    float env1Inc;
    float env2Inc;
    float env3Inc;
    float env4Inc;
    float env5Inc;
    float env6Inc;

    float *sample = sampleBlock;

    switch ((int)currentTimbre->params_.engine1.algo)
    {

    case ALGO1:
        /*
         IM3
         <----
         .---.  .---.
         | 2 |  | 3*|
         '---'  '---'
         \IM1   /IM2
         .---.
         | 1 |
         '---'
         */

        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlockWithFeedbackAndEnveloppe(&oscState3_, feedbackModulation, env3Value, env3Inc,
                                                                                          oscState3_.frequency, fdbLastValue);

            float f2x;
            float f2xm1 = freqAi;
            float freq2 = freqAo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {
                // getNextBlockWithFeedbackAndEnveloppe already applied envelope and multiplier
                float freq3 = osc3Values[k];

                oscState2_.frequency = freq3 * voiceIm3 + oscState2_.mainFrequencyPlusMatrix;
                f2x = currentTimbre->osc2_.getNextSample(&oscState2_) * env2Value * oscState2_.frequency;
                freq2 = f2x - f2xm1 + dcPole * freq2;
                f2xm1 = f2x;

                oscState1_.frequency = oscState1_.mainFrequencyPlusMatrix + voiceIm1 * freq2 + voiceIm2 * freq3;
                float currentSample = currentTimbre->osc1_.getNextSample(&oscState1_) * this->velocity * env1Value * mix1;

                *sample++ = currentSample * pan1Right;
                *sample++ = currentSample * pan1Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
            }

            freqAi = f2xm1;
            freqAo = freq2;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }

    case ALGO2:
        /*
         .---.
         | 3 |
         '---'
         |
         .------.
         |IM1   |IM2
         .---.  .---.
         | 1 |  | 2 |
         '---'  '---'
         |Mix1  |Mix2
         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlockWithFeedbackAndEnveloppe(&oscState3_, feedbackModulation, env3Value, env3Inc,
                                                                                          oscState3_.frequency, fdbLastValue);

            float div2TimesVelocity = this->velocity * .5f;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {
                // getNextBlockWithFeedbackAndEnveloppe already applied envelope and multiplier
                float freq3 = osc3Values[k];

                oscState2_.frequency = freq3 * voiceIm2 + oscState2_.mainFrequencyPlusMatrix;
                float car2 = currentTimbre->osc2_.getNextSample(&oscState2_) * env2Value * mix2 * div2TimesVelocity;

                oscState1_.frequency = freq3 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div2TimesVelocity;

                *sample++ = car1 * pan1Right + car2 * pan2Right;
                *sample++ = car1 * pan1Left + car2 * pan2Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
            }

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env2_.isDead(&envState2_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
        }
        break;

    case ALGO3:
        /*

         IM4
         ---->
         .---.  .---.     .---.
         | 2 |  | 3 |     | 4*|
         '---'  '---'     '---'
         \IM1  |IM2    /IM3
         .---.
         | 1 |
         '---'
         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlockWithFeedbackAndEnveloppe(&oscState3_, feedbackModulation, env3Value, env3Inc,
                                                                                          oscState3_.frequency, fdbLastValue);

            float f4x;
            float f4xm1 = freqAi;
            float freq4 = freqAo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {
                float freq2 = osc2Values[k] * env2Value * oscState2_.frequency;
                float freq3 = osc3Values[k];

                oscState4_.frequency = freq3 * voiceIm4 + oscState4_.mainFrequencyPlusMatrix;

                f4x = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * oscState4_.frequency;
                freq4 = f4x - f4xm1 + dcPole * freq4;
                f4xm1 = f4x;

                oscState1_.frequency = freq2 * voiceIm1 + freq3 * voiceIm2 + freq4 * voiceIm3 + oscState1_.mainFrequencyPlusMatrix;
                float currentSample = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * this->velocity;

                *sample++ = currentSample * pan1Right;
                *sample++ = currentSample * pan1Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env4Value += env4Inc;
            }

            freqAi = f4xm1;
            freqAo = freq4;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }
    case ALGO4:
        /*           IM4
         .---. <----   .---.
         | 3 |         | 4 |
         '---'         '---'
         |IM1 \IM3     |IM2
         .---.          .---.
         | 1 |          | 2 |
         '---'          '---'
         |Mix1          |Mix2

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlockWithFeedbackAndEnveloppe(&oscState4_, feedbackModulation, env4Value, env4Inc,
                                                                                          oscState4_.frequency, fdbLastValue);

            float div2TimesVelocity = this->velocity * .5f;

            float f3x;
            float f3xm1 = freqAi;
            float freq3 = freqAo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {
                float freq4 = osc4Values[k];

                oscState3_.frequency = freq4 * voiceIm4 + oscState3_.mainFrequencyPlusMatrix;
                f3x = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * oscState3_.frequency;
                freq3 = f3x - f3xm1 + dcPole * freq3;
                f3xm1 = f3x;

                oscState2_.frequency = freq4 * voiceIm2 + freq3 * voiceIm3 + oscState2_.mainFrequencyPlusMatrix;
                float car2 = currentTimbre->osc2_.getNextSample(&oscState2_) * env2Value * mix2 * div2TimesVelocity;

                oscState1_.frequency = freq3 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div2TimesVelocity;

                *sample++ = car1 * pan1Right + car2 * pan2Right;
                *sample++ = car1 * pan1Left + car2 * pan2Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
            }

            freqAi = f3xm1;
            freqAo = freq3;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env2_.isDead(&envState2_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
            break;
        }
    case ALGO5:
        /*
         .---.
         | 4 |
         '---'  \
               |IM3  |
         .---.   |
         | 3 |   | IM4
         '---'   |
         |IM2  |
         .---.  /
         | 2 |
         '---'
         |IM1
         .---.
         | 1 |
         '---'

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlockWithFeedbackAndEnveloppe(&oscState4_, feedbackModulation, env4Value, env4Inc,
                                                                                          oscState4_.frequency, fdbLastValue);

            float f2x;
            float f2xm1 = freqAi;
            float freq2 = freqAo;
            float f3x;
            float f3xm1 = freqBi;
            float freq3 = freqBo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq4 = osc4Values[k];

                oscState3_.frequency = freq4 * voiceIm3 + oscState3_.mainFrequencyPlusMatrix;
                f3x = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * oscState3_.frequency;
                freq3 = f3x - f3xm1 + 0.99535f * freq3;
                f3xm1 = f3x;

                oscState2_.frequency = freq3 * voiceIm2 + freq4 * voiceIm4 + oscState2_.mainFrequencyPlusMatrix;
                f2x = currentTimbre->osc2_.getNextSample(&oscState2_) * env2Value * oscState2_.frequency;
                freq2 = f2x - f2xm1 + dcPole * freq2;
                f2xm1 = f2x;

                oscState1_.frequency = freq2 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * this->velocity * mix1;

                *sample++ = car1 * pan1Right;
                *sample++ = car1 * pan1Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
            }
            freqAi = f2xm1;
            freqAo = freq2;
            freqBi = f3xm1;
            freqBo = freq3;
            if (unlikely(currentTimbre->env1_.isDead(&envState1_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
            break;
        }
    case ALGO6:
        /*
         .---.
         | 4 |
         '---'
         /IM1 |IM2 \IM3
         .---.  .---.  .---.
         | 1 |  | 2 |  | 3 |
         '---'  '---'  '---'
         |Mix1  |Mix2  | Mix3

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlockWithFeedbackAndEnveloppe(&oscState4_, feedbackModulation, env4Value, env4Inc,
                                                                                          oscState4_.frequency, fdbLastValue);

            float div3TimesVelocity = .33f * this->velocity;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq4 = osc4Values[k];

                oscState3_.frequency = freq4 * voiceIm3 + oscState3_.mainFrequencyPlusMatrix;
                float car3 = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * mix3 * div3TimesVelocity;

                oscState2_.frequency = freq4 * voiceIm2 + oscState2_.mainFrequencyPlusMatrix;
                float car2 = currentTimbre->osc2_.getNextSample(&oscState2_) * env2Value * mix2 * div3TimesVelocity;

                oscState1_.frequency = freq4 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div3TimesVelocity;

                *sample++ = car1 * pan1Right + car2 * pan2Right + car3 * pan3Right;
                *sample++ = car1 * pan1Left + car2 * pan2Left + car3 * pan3Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
            }

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env2_.isDead(&envState2_) && currentTimbre->env3_.isDead(&envState3_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }

    case ALGO7:
        /*
         IM4
         ---->
         .---.  .---.  .---.
         | 2 |  | 4 |  | 6 |
         '---'  '---'  '---'
         |IM1   |IM2   |IM3
         .---.  .---.  .---.
         | 1 |  | 3 |  | 5 |
         '---'  '---'  '---'
         |Mix1  |Mix2  |Mix3

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlockWithFeedbackAndEnveloppe(&oscState4_, feedbackModulation, env4Value, env4Inc,
                                                                                          oscState4_.frequency, fdbLastValue);

            float div3TimesVelocity = .33f * this->velocity;

            // F2 will be for OP6
            float f6x;
            float f6xm1 = freqAi;
            float freq6 = freqAo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq4 = osc4Values[k];
                float freq2 = osc2Values[k] * env2Value * oscState2_.frequency;

                oscState6_.frequency = freq4 * voiceIm4 + oscState6_.mainFrequencyPlusMatrix;
                f6x = currentTimbre->osc6_.getNextSample(&oscState6_) * env6Value * oscState6_.frequency;
                freq6 = f6x - f6xm1 + dcPole * freq6;
                f6xm1 = f6x;

                oscState1_.frequency = freq2 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div3TimesVelocity;

                oscState3_.frequency = freq4 * voiceIm2 + oscState3_.mainFrequencyPlusMatrix;
                float car3 = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * mix2 * div3TimesVelocity;

                oscState5_.frequency = freq6 * voiceIm3 + oscState5_.mainFrequencyPlusMatrix;
                float car5 = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * mix3 * div3TimesVelocity;

                *sample++ = car1 * pan1Right + car3 * pan2Right + car5 * pan3Right;
                *sample++ = car1 * pan1Left + car3 * pan2Left + car5 * pan3Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env5Value += env5Inc;
                env6Value += env6Inc;
            }

            freqAi = f6xm1;
            freqAo = freq6;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env3_.isDead(&envState3_) && currentTimbre->env5_.isDead(&envState5_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
        }
        break;
    case ALGO8:
        /*
         .---.  .---.  .---.       .---.
         | 2 |  | 3 |  | 4 |       | 6 |
         '---'  '---'  '---'       '---'
         \IM1  |IM2  /IM3         | IM4
         .---.              .---.
         | 1 |              | 5 |
         '---'              '---'
         |Mix1              | Mix2

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlock(&oscState3_);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlockWithFeedbackAndEnveloppe(&oscState4_, feedbackModulation, env4Value, env4Inc,
                                                                                          oscState4_.frequency, fdbLastValue);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlock(&oscState6_);

            float div2TimesVelocity = .5f * this->velocity;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq6 = osc6Values[k] * env6Value;
                freq6 *= oscState6_.frequency;

                float freq4 = osc4Values[k];

                float freq3 = osc3Values[k] * env3Value;
                freq3 *= oscState3_.frequency;

                float freq2 = osc2Values[k] * env2Value;
                freq2 *= oscState2_.frequency;

                oscState1_.frequency = freq2 * voiceIm1 + freq3 * voiceIm2 + freq4 * voiceIm3 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div2TimesVelocity;

                oscState5_.frequency = freq6 * voiceIm4 + oscState5_.mainFrequencyPlusMatrix;
                float car5 = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * mix2 * div2TimesVelocity;

                *sample++ = car1 * pan1Right + car5 * pan2Right;
                *sample++ = car1 * pan1Left + car5 * pan2Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env5Value += env5Inc;
                env6Value += env6Inc;
            }

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env5_.isDead(&envState5_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }

    case ALGO9:
        /*
         .---.
         | 6 |
         '---'
         |IM4
         .---.      .---.       .---.
         | 2 |      | 3 |       | 5 |
         '---'      '---'       '---'
         \IM1    /IM2         | IM3
         .---.            .---.
         | 1 |            | 4 |
         '---'            '---'
         |Mix1            | Mix2

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlock(&oscState3_);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlockWithFeedbackAndEnveloppe(&oscState6_, feedbackModulation, env6Value, env6Inc,
                                                                                          oscState6_.frequency, fdbLastValue);

            float div2TimesVelocity = .5f * this->velocity;

            // F2 will be for OP5
            float f5x;
            float f5xm1 = freqAi;
            float freq5 = freqAo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq6 = osc6Values[k];

                float freq3 = osc3Values[k] * env3Value;
                freq3 *= oscState3_.frequency;

                float freq2 = osc2Values[k] * env2Value;
                freq2 *= oscState2_.frequency;

                oscState1_.frequency = freq2 * voiceIm1 + freq3 * voiceIm2 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div2TimesVelocity;

                oscState5_.frequency = freq6 * voiceIm4 + oscState5_.mainFrequencyPlusMatrix;
                f5x = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * oscState5_.frequency;
                freq5 = f5x - f5xm1 + dcPole * freq5;
                f5xm1 = f5x;

                oscState4_.frequency = freq5 * voiceIm3 + oscState4_.mainFrequencyPlusMatrix;
                float car4 = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * mix2 * div2TimesVelocity;

                *sample++ = car4 * pan2Right + car1 * pan1Right;
                *sample++ = car4 * pan2Left + car1 * pan1Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
            }
            freqAi = f5xm1;
            freqAo = freq5;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env4_.isDead(&envState4_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }

        // =========== BELOW THIS LINE DX7 ALGORITHMS !!! =================

    case ALG10:
        /* DX7 Algo 1 & 2
         .---.
         | 6 |
         '---'
         |IM4
         .---.
         | 5 |
         '---'
         |IM3
         .---.            .---.
         | 2 |            | 4 |
         '---'            '---'
         |IM1             | IM2
         .---.            .---.
         | 1 |            | 3 |
         '---'            '---'
         |Mix1            | Mix2

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlockWithFeedbackAndEnveloppe(&oscState2_, feedbackModulation, env2Value, env2Inc,
                                                                                          oscState2_.frequency, fdbLastValue);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlock(&oscState6_);

            float div2TimesVelocity = .5f * this->velocity;

            float f4x;
            float f4xm1 = freqAi;
            float freq4 = freqAo;

            float f5x;
            float f5xm1 = freqBi;
            float freq5 = freqBo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq2 = osc2Values[k];

                oscState1_.frequency = freq2 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div2TimesVelocity;

                float freq6 = osc6Values[k] * env6Value;
                freq6 *= oscState6_.frequency;

                oscState5_.frequency = freq6 * voiceIm4 + oscState5_.mainFrequencyPlusMatrix;
                f5x = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * oscState5_.frequency;
                freq5 = f5x - f5xm1 + dcPole * freq5;
                f5xm1 = f5x;

                oscState4_.frequency = freq5 * voiceIm3 + oscState4_.mainFrequencyPlusMatrix;
                f4x = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * oscState4_.frequency;
                freq4 = f4x - f4xm1 + dcPole * freq4;
                f4xm1 = f4x;

                oscState3_.frequency = freq4 * voiceIm2 + oscState3_.mainFrequencyPlusMatrix;
                float car3 = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * mix2 * div2TimesVelocity;

                *sample++ = car3 * pan2Right + car1 * pan1Right;
                *sample++ = car3 * pan2Left + car1 * pan1Left;

                env1Value += env1Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
                env6Value += env6Inc;
            }

            freqBi = f5xm1;
            freqBo = freq5;

            freqAi = f4xm1;
            freqAo = freq4;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env3_.isDead(&envState3_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }
    case ALG11:
        /*
         * DX7 algo 3 & 4
         .---.            .---.
         | 3 |            | 6 |
         '---'            '---'
         |IM2             |IM4
         .---.            .---.
         | 2 |            | 5 |
         '---'            '---'
         |IM1            | IM3
         .---.            .---.
         | 1 |            | 4 |
         '---'            '---'
         |Mix1            | Mix2

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlock(&oscState3_);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlockWithFeedbackAndEnveloppe(&oscState6_, feedbackModulation, env6Value, env6Inc,
                                                                                          oscState6_.frequency, fdbLastValue);

            float div2TimesVelocity = .5f * this->velocity;

            float f2x;
            float f2xm1 = freqAi;
            float freq2 = freqAo;

            float f5x;
            float f5xm1 = freqBi;
            float freq5 = freqBo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq3 = osc3Values[k] * env3Value;
                freq3 *= oscState3_.frequency;

                oscState2_.frequency = freq3 * voiceIm2 + oscState2_.mainFrequencyPlusMatrix;
                f2x = currentTimbre->osc2_.getNextSample(&oscState2_) * env2Value * oscState2_.frequency;
                freq2 = f2x - f2xm1 + dcPole * freq2;
                f2xm1 = f2x;

                oscState1_.frequency = freq2 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div2TimesVelocity;

                float freq6 = osc6Values[k];

                oscState5_.frequency = freq6 * voiceIm4 + oscState5_.mainFrequencyPlusMatrix;
                f5x = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * oscState5_.frequency;
                freq5 = f5x - f5xm1 + dcPole * freq5;
                f5xm1 = f5x;

                oscState4_.frequency = freq5 * voiceIm3 + oscState4_.mainFrequencyPlusMatrix;
                float car4 = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * mix2 * div2TimesVelocity;

                *sample++ = car4 * pan2Right + car1 * pan1Right;
                *sample++ = car4 * pan2Left + car1 * pan1Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
            }

            freqAi = f2xm1;
            freqAo = freq2;

            freqBi = f5xm1;
            freqBo = freq5;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env4_.isDead(&envState4_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }
    case ALG12:
        /*
         * DX7 algo 5 & 6

         .---.  .---.  .---.
         | 2 |  | 4 |  | 6 |
         '---'  '---'  '---'
         |IM1   |IM2   |IM3
         .---.  .---.  .---.
         | 1 |  | 3 |  | 5 |
         '---'  '---'  '---'
         |Mix1  |Mix2  |Mix3

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlock(&oscState4_);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlockWithFeedbackAndEnveloppe(&oscState6_, feedbackModulation, env6Value, env6Inc,
                                                                                          oscState6_.frequency, fdbLastValue);

            float div3TimesVelocity = .33f * this->velocity;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq4 = osc4Values[k] * env4Value;
                freq4 *= oscState4_.frequency;

                float freq2 = osc2Values[k] * env2Value;
                freq2 *= oscState2_.frequency;

                float freq6 = osc6Values[k];

                oscState1_.frequency = freq2 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div3TimesVelocity;

                oscState3_.frequency = freq4 * voiceIm2 + oscState3_.mainFrequencyPlusMatrix;
                float car3 = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * mix2 * div3TimesVelocity;

                oscState5_.frequency = freq6 * voiceIm3 + oscState5_.mainFrequencyPlusMatrix;
                float car5 = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * mix3 * div3TimesVelocity;

                *sample++ = car1 * pan1Right + car3 * pan2Right + car5 * pan3Right;
                *sample++ = car1 * pan1Left + car3 * pan2Left + car5 * pan3Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
            }
            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env3_.isDead(&envState3_) && currentTimbre->env5_.isDead(&envState5_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
        }
        break;
    case ALG13:
        /* DX7 Algo 7, 8, 9

         .---.
         | 6 |
         '---'
         |IM4
         .---.     .---.  .---.
         | 2 |     | 4*|  | 5 |
         '---'     '---'  '---'
         |IM1      |IM2 /IM3
         .---.     .---.
         | 1 |     | 3 |
         '---'     '---'
         |Mix1     |Mix2
         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlockWithFeedbackAndEnveloppe(&oscState4_, feedbackModulation, env4Value, env4Inc,
                                                                                          oscState4_.frequency, fdbLastValue);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlock(&oscState6_);

            float div2TimesVelocity = this->velocity * .5f;

            float f5x;
            float f5xm1 = freqAi;
            float freq5 = freqAo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq6 = osc6Values[k] * env6Value;
                freq6 *= oscState6_.frequency;

                float freq4 = osc4Values[k];

                float freq2 = osc2Values[k] * env2Value;
                freq2 *= oscState2_.frequency;

                oscState1_.frequency = freq2 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div2TimesVelocity;

                oscState5_.frequency = freq6 * voiceIm4 + oscState5_.mainFrequencyPlusMatrix;
                f5x = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * oscState5_.frequency;
                freq5 = f5x - f5xm1 + dcPole * freq5;
                f5xm1 = f5x;

                oscState3_.frequency = freq4 * voiceIm2 + freq5 * voiceIm3 + oscState3_.mainFrequencyPlusMatrix;
                float car3 = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * mix2 * div2TimesVelocity;

                *sample++ = car1 * pan1Right + car3 * pan2Right;
                *sample++ = car1 * pan1Left + car3 * pan2Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env5Value += env5Inc;
                env6Value += env6Inc;
            }

            freqAi = f5xm1;
            freqAo = freq5;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env3_.isDead(&envState3_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }
    case ALG14:
        /* DX Algo 10 & 11
         .---.
         | 3 |
         '---'
         |IM2
         .---.            .---.   .---.
         | 2 |            | 5 |   | 6*|
         '---'            '---'   '---'
         |IM1            | IM3 / IM4
         .---.            .---.
         | 1 |            | 4 |
         '---'            '---'
         |Mix1            | Mix2

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlock(&oscState3_);

            oscState5_.frequency = oscState5_.mainFrequencyPlusMatrix;
            float *osc5Values = currentTimbre->osc5_.getNextBlock(&oscState5_);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlockWithFeedbackAndEnveloppe(&oscState6_, feedbackModulation, env6Value, env6Inc,
                                                                                          oscState6_.frequency, fdbLastValue);

            float div2TimesVelocity = .5f * this->velocity;

            float f2x;
            float f2xm1 = freqAi;
            float freq2 = freqAo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq3 = osc3Values[k] * env3Value;
                freq3 *= oscState3_.frequency;

                float freq5 = osc5Values[k] * env5Value;
                freq5 *= oscState5_.frequency;

                float freq6 = osc6Values[k];

                oscState2_.frequency = freq3 * voiceIm2 + oscState2_.mainFrequencyPlusMatrix;
                f2x = currentTimbre->osc2_.getNextSample(&oscState2_) * env2Value * oscState2_.frequency;
                freq2 = f2x - f2xm1 + dcPole * freq2;
                f2xm1 = f2x;

                oscState1_.frequency = freq2 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div2TimesVelocity;

                oscState4_.frequency = freq5 * voiceIm3 + freq6 * voiceIm4 + oscState4_.mainFrequencyPlusMatrix;
                float car4 = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * mix2 * div2TimesVelocity;

                *sample++ = car4 * pan2Right + car1 * pan1Right;
                *sample++ = car4 * pan2Left + car1 * pan1Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
            }

            freqAi = f2xm1;
            freqAo = freq2;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env4_.isDead(&envState4_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }
    case ALG15:
        /*
         * DX Algo 12 & 13


         .---.     .---.   .---.   .---.
         | 2 |     | 4 |   | 5 |   | 6*|
         '---'     '---'   '---'   '---'
         |IM1        \IM2  | IM3 / IM4
         .---.             .---.
         | 1 |             | 3 |
         '---'             '---'
         |Mix1             | Mix2

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlock(&oscState4_);

            oscState5_.frequency = oscState5_.mainFrequencyPlusMatrix;
            float *osc5Values = currentTimbre->osc5_.getNextBlock(&oscState5_);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlockWithFeedbackAndEnveloppe(&oscState6_, feedbackModulation, env6Value, env6Inc,
                                                                                          oscState6_.frequency, fdbLastValue);

            float div2TimesVelocity = .5f * this->velocity;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq2 = osc2Values[k] * env2Value;
                freq2 *= oscState2_.frequency;

                oscState1_.frequency = freq2 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div2TimesVelocity;

                float freq4 = osc4Values[k] * env4Value;
                freq4 *= oscState4_.frequency;
                float freq5 = osc5Values[k] * env5Value;
                freq5 *= oscState5_.frequency;
                float freq6 = osc6Values[k];

                oscState3_.frequency = freq4 * voiceIm2 + freq5 * voiceIm3 + freq6 * voiceIm4 + oscState3_.mainFrequencyPlusMatrix;
                float car3 = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * mix2 * div2TimesVelocity;

                *sample++ = car3 * pan2Right + car1 * pan1Right;
                *sample++ = car3 * pan2Left + car1 * pan1Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
            }

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env3_.isDead(&envState3_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }
    case ALG16:
        /* DX7 Algo 14 and 15
         .---.  .---.
         | 5 |  | 6 |
         '---'  '---'
         |IM3/IM4
         .---.            .---.
         | 2*|            | 4 |
         '---'            '---'
         |IM1             | IM2
         .---.            .---.
         | 1 |            | 3 |
         '---'            '---'
         |Mix1            | Mix2

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlockWithFeedbackAndEnveloppe(&oscState2_, feedbackModulation, env2Value, env2Inc,
                                                                                          oscState2_.frequency, fdbLastValue);

            oscState5_.frequency = oscState5_.mainFrequencyPlusMatrix;
            float *osc5Values = currentTimbre->osc5_.getNextBlock(&oscState5_);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlock(&oscState6_);

            float div2TimesVelocity = .5f * this->velocity;

            float f4x;
            float f4xm1 = freqAi;
            float freq4 = freqAo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq2 = osc2Values[k];

                oscState1_.frequency = freq2 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div2TimesVelocity;

                float freq6 = osc6Values[k] * env6Value;
                freq6 *= oscState6_.frequency;

                float freq5 = osc5Values[k] * env5Value;
                freq5 *= oscState5_.frequency;

                oscState4_.frequency = freq5 * voiceIm3 + freq6 * voiceIm4 + oscState4_.mainFrequencyPlusMatrix;
                f4x = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * oscState4_.frequency;
                freq4 = f4x - f4xm1 + dcPole * freq4;
                f4xm1 = f4x;

                oscState3_.frequency = freq4 * voiceIm2 + oscState3_.mainFrequencyPlusMatrix;
                float car3 = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * mix2 * div2TimesVelocity;

                *sample++ = car3 * pan2Right + car1 * pan1Right;
                *sample++ = car3 * pan2Left + car1 * pan1Left;

                env1Value += env1Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
                env6Value += env6Inc;
            }

            freqAi = f4xm1;
            freqAo = freq4;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env3_.isDead(&envState3_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }
    case ALG17:
        /*
         * DX ALGO 16 & 17
         .---.     .---.
         | 4 |     | 6 |
         '---'     '---'
         |IM3      |IM5
         .---.     .---.     .---.
         | 2*|     | 3 |     | 5 |
         '---'     '---'     '---'
         \IM1    |IM2   / IM4
         .---.
         | 1 |
         '---'
         |Mix1

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;
            float voiceIm5 = modulationIndex5;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlockWithFeedbackAndEnveloppe(&oscState2_, feedbackModulation, env2Value, env2Inc,
                                                                                          oscState2_.frequency, fdbLastValue);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlock(&oscState4_);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlock(&oscState6_);

            float f3x;
            float f3xm1 = freqAi;
            float freq3 = freqAo;

            float f5x;
            float f5xm1 = freqBi;
            float freq5 = freqBo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq6 = osc6Values[k] * env6Value;
                freq6 *= oscState6_.frequency;

                float freq4 = osc4Values[k] * env4Value;
                freq4 *= oscState4_.frequency;

                float freq2 = osc2Values[k];

                oscState3_.frequency = freq4 * voiceIm3 + oscState3_.mainFrequencyPlusMatrix;
                f3x = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * oscState3_.frequency;
                freq3 = f3x - f3xm1 + dcPole * freq3;
                f3xm1 = f3x;

                oscState5_.frequency = freq6 * voiceIm5 + oscState5_.mainFrequencyPlusMatrix;
                f5x = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * oscState5_.frequency;
                freq5 = f5x - f5xm1 + dcPole * freq5;
                f5xm1 = f5x;

                oscState1_.frequency = freq2 * voiceIm1 + freq3 * voiceIm2 + freq5 * voiceIm4 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * velocity;

                *sample++ = car1 * pan1Right;
                *sample++ = car1 * pan1Left;

                env1Value += env1Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
                env6Value += env6Inc;
            }

            freqAi = f3xm1;
            freqAo = freq3;

            freqBi = f5xm1;
            freqBo = freq5;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }
    case ALG18:
        /*
         * DX ALGO 18
         .---.
         | 6 |
         '---'
         |IM5
         .---.
         | 5 |
         '---'
         |IM4
         .---.     .---.     .---.
         | 2 |     | 3*|     | 4 |
         '---'     '---'     '---'
         \IM1    |IM2   / IM3
         .---.
         | 1 |
         '---'
         |Mix1
         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;
            float voiceIm5 = modulationIndex5;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlockWithFeedbackAndEnveloppe(&oscState3_, feedbackModulation, env3Value, env3Inc,
                                                                                          oscState3_.frequency, fdbLastValue);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlock(&oscState6_);

            float f4x;
            float f4xm1 = freqAi;
            float freq4 = freqAo;

            float f5x;
            float f5xm1 = freqBi;
            float freq5 = freqBo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq6 = osc6Values[k] * env6Value;
                freq6 *= oscState6_.frequency;

                float freq3 = osc3Values[k];

                float freq2 = osc2Values[k] * env2Value;
                freq2 *= oscState2_.frequency;

                oscState5_.frequency = freq6 * voiceIm5 + oscState5_.mainFrequencyPlusMatrix;
                f5x = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * oscState5_.frequency;
                freq5 = f5x - f5xm1 + dcPole * freq5;
                f5xm1 = f5x;

                oscState4_.frequency = freq5 * voiceIm4 + oscState4_.mainFrequencyPlusMatrix;
                f4x = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * oscState4_.frequency;
                freq4 = f4x - f4xm1 + dcPole * freq4;
                f4xm1 = f4x;

                oscState1_.frequency = freq2 * voiceIm1 + freq3 * voiceIm2 + freq4 * voiceIm3 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * velocity;

                *sample++ = car1 * pan1Right;
                *sample++ = car1 * pan1Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
                env6Value += env6Inc;
            }

            freqAi = f4xm1;
            freqAo = freq4;

            freqBi = f5xm1;
            freqBo = freq5;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }

    case ALG19:
        /*
         * DX7 algo 19

         .---.
         | 3 |
         '---'
         |IM2
         .---.        .---.
         | 2 |        | 6*|
         '---'        '---'
         |IM1    IM3/   \IM4
         .---.  .---.     .---.
         | 1 |  | 4 |     | 5 |
         '---'  '---'     '---'
         |Mix1  |Mix2     |Mix3

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlock(&oscState3_);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlockWithFeedbackAndEnveloppe(&oscState6_, feedbackModulation, env6Value, env6Inc,
                                                                                          oscState6_.frequency, fdbLastValue);

            float div3TimesVelocity = .33f * this->velocity;

            float f2x;
            float f2xm1 = freqAi;
            float freq2 = freqAo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq3 = osc3Values[k] * env3Value;
                freq3 *= oscState3_.frequency;

                float freq6 = osc6Values[k];

                oscState2_.frequency = freq3 * voiceIm2 + oscState2_.mainFrequencyPlusMatrix;
                f2x = currentTimbre->osc2_.getNextSample(&oscState2_) * env2Value * oscState2_.frequency;
                freq2 = f2x - f2xm1 + dcPole * freq2;
                f2xm1 = f2x;

                oscState1_.frequency = freq2 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div3TimesVelocity;

                oscState4_.frequency = freq6 * voiceIm3 + oscState4_.mainFrequencyPlusMatrix;
                float car4 = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * mix2 * div3TimesVelocity;

                oscState5_.frequency = freq6 * voiceIm4 + oscState5_.mainFrequencyPlusMatrix;
                float car5 = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * mix3 * div3TimesVelocity;

                *sample++ = car1 * pan1Right + car4 * pan2Right + car5 * pan3Right;
                *sample++ = car1 * pan1Left + car4 * pan2Left + car5 * pan3Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
            }

            freqAi = f2xm1;
            freqAo = freq2;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env4_.isDead(&envState4_) && currentTimbre->env5_.isDead(&envState5_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
        }
        break;
    case ALG20:
        /*
         * DX7 algo 20, 26 & 27

         .---.     .---.  .---.
         | 3*|     | 5 |  | 6 |
         '---'     '---'  '---'
         /IM1 \IM2     \IM3 /IM4
         .---.  .---.    .---.
         | 1 |  | 2 |    | 4 |
         '---'  '---'    '---'
         |Mix1  |Mix2    |Mix3
         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlockWithFeedbackAndEnveloppe(&oscState3_, feedbackModulation, env3Value, env3Inc,
                                                                                          oscState3_.frequency, fdbLastValue);

            oscState5_.frequency = oscState5_.mainFrequencyPlusMatrix;
            float *osc5Values = currentTimbre->osc5_.getNextBlock(&oscState5_);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlock(&oscState6_);

            float div3TimesVelocity = .33f * this->velocity;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq3 = osc3Values[k];

                float freq5 = osc5Values[k] * env5Value;
                freq5 *= oscState5_.frequency;

                float freq6 = osc6Values[k] * env6Value;
                freq6 *= oscState6_.frequency;

                oscState1_.frequency = freq3 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div3TimesVelocity;

                oscState2_.frequency = freq3 * voiceIm2 + oscState2_.mainFrequencyPlusMatrix;
                float car2 = currentTimbre->osc2_.getNextSample(&oscState2_) * env2Value * mix2 * div3TimesVelocity;

                oscState4_.frequency = freq5 * voiceIm3 + freq6 * voiceIm4 + oscState4_.mainFrequencyPlusMatrix;
                float car4 = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * mix3 * div3TimesVelocity;

                *sample++ = car1 * pan1Right + car2 * pan2Right + car4 * pan3Right;
                *sample++ = car1 * pan1Left + car2 * pan2Left + car4 * pan3Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
                env6Value += env6Inc;
            }

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env2_.isDead(&envState2_) && currentTimbre->env4_.isDead(&envState4_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
        }
        break;
    case ALG21:
        /* DX ALGO 21 & 23
         .---.         .---.
         | 3*|         | 6 |
         '---'         '---'
         /IM1  \IM2   /IM3  \IM4
         .---.  .---.  .---.  .---.
         | 1 |  | 2 |  | 4 |  | 5 |
         '---'  '---'  '---'  '---'
         |Mix1  |Mix2  |Mix3  |Mix4

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlockWithFeedbackAndEnveloppe(&oscState3_, feedbackModulation, env3Value, env3Inc,
                                                                                          oscState3_.frequency, fdbLastValue);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlock(&oscState6_);

            float div4TimesVelocity = .25f * this->velocity;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq6 = osc6Values[k] * env6Value;
                freq6 *= oscState6_.frequency;

                float freq3 = osc3Values[k];

                oscState1_.frequency = freq3 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div4TimesVelocity;

                oscState2_.frequency = freq3 * voiceIm2 + oscState2_.mainFrequencyPlusMatrix;
                float car2 = currentTimbre->osc2_.getNextSample(&oscState2_) * env2Value * mix2 * div4TimesVelocity;

                oscState4_.frequency = freq6 * voiceIm3 + oscState4_.mainFrequencyPlusMatrix;
                float car4 = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * mix3 * div4TimesVelocity;

                oscState5_.frequency = freq6 * voiceIm4 + oscState5_.mainFrequencyPlusMatrix;
                float car5 = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * mix4 * div4TimesVelocity;

                *sample++ = car1 * pan1Right + car2 * pan2Right + car4 * pan3Right + car5 * pan4Right;
                *sample++ = car1 * pan1Left + car2 * pan2Left + car4 * pan3Left + car5 * pan4Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
                env6Value += env6Inc;
            }
            if (unlikely(
                    currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env2_.isDead(&envState2_) && currentTimbre->env4_.isDead(&envState4_) && currentTimbre->env5_.isDead(&envState5_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
        }
        break;
    case ALG22:
        /* DX7 algo 22
         .---.         .---.
         | 2 |         | 6*|
         '---'         '---'
         |IM1      /IM2 |IM3 \IM4
         .---.  .---.  .---.  .---.
         | 1 |  | 3 |  | 4 |  | 5 |
         '---'  '---'  '---'  '---'
         |Mix1  |Mix2  |Mix3  |Mix4

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlockWithFeedbackAndEnveloppe(&oscState6_, feedbackModulation, env6Value, env6Inc,
                                                                                          oscState6_.frequency, fdbLastValue);

            float div4TimesVelocity = .25f * this->velocity;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq6 = osc6Values[k];

                float freq2 = osc2Values[k] * env2Value;
                freq2 *= oscState2_.frequency;

                oscState1_.frequency = freq2 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div4TimesVelocity;

                oscState3_.frequency = freq6 * voiceIm2 + oscState3_.mainFrequencyPlusMatrix;
                float car3 = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * mix2 * div4TimesVelocity;

                oscState4_.frequency = freq6 * voiceIm3 + oscState4_.mainFrequencyPlusMatrix;
                float car4 = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * mix3 * div4TimesVelocity;

                oscState5_.frequency = freq6 * voiceIm4 + oscState5_.mainFrequencyPlusMatrix;
                float car5 = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * mix4 * div4TimesVelocity;

                *sample++ = car1 * pan1Right + car3 * pan2Right + car4 * pan3Right + car5 * pan4Right;
                *sample++ = car1 * pan1Left + car3 * pan2Left + car4 * pan3Left + car5 * pan4Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
            }
            if (unlikely(
                    currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env3_.isDead(&envState3_) && currentTimbre->env4_.isDead(&envState4_) && currentTimbre->env5_.isDead(&envState5_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
        }
        break;
    case ALG23:
        /* DX ALGO 24, 25 & 31

         .---.
         | 6 |
         '---'
         /IM1|IM2\IM3
         .---.  .---.  .---.  .---.   .---.
         | 1 |  | 2 |  | 3 |  | 4 |   | 5 |
         '---'  '---'  '---'  '---'   '---'
         |Mix1  |Mix2  |Mix3  |Mix4   |Mix5

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState1_.frequency = oscState1_.mainFrequencyPlusMatrix;
            float *osc1Values = currentTimbre->osc1_.getNextBlock(&oscState1_);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlockWithFeedbackAndEnveloppe(&oscState6_, feedbackModulation, env6Value, env6Inc,
                                                                                          oscState6_.frequency, fdbLastValue);

            float div5TimesVelocity = .2f * this->velocity;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq6 = osc6Values[k];

                float car1 = osc1Values[k] * env1Value * mix1 * div5TimesVelocity;

                float car2 = osc2Values[k] * env2Value * mix2 * div5TimesVelocity;

                oscState3_.frequency = freq6 * voiceIm1 + oscState3_.mainFrequencyPlusMatrix;
                float car3 = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * mix3 * div5TimesVelocity;

                oscState4_.frequency = freq6 * voiceIm2 + oscState4_.mainFrequencyPlusMatrix;
                float car4 = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * mix4 * div5TimesVelocity;

                oscState5_.frequency = freq6 * voiceIm3 + oscState5_.mainFrequencyPlusMatrix;
                float car5 = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * mix5 * div5TimesVelocity;

                *sample++ = car1 * pan1Right + car2 * pan2Right + car3 * pan3Right + car4 * pan4Right + car5 * pan5Right;
                *sample++ = car1 * pan1Left + car2 * pan2Left + car3 * pan3Left + car4 * pan4Left + car5 * pan5Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
            }
            if (unlikely(
                    currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env2_.isDead(&envState2_) && currentTimbre->env4_.isDead(&envState4_) && currentTimbre->env5_.isDead(&envState5_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
        }
        break;
    case ALG24:
        /*
         * DX ALGO 28

         .---.
         | 5*|
         '---'
         |IM3
         .---.  .---.
         | 2 |  | 4 |
         '---'  '---'
         |IM1   |IM2
         .---.  .---.  .---.
         | 1 |  | 3 |  | 6 |
         '---'  '---'  '---'
         |Mix1  |Mix2  |Mix3


         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState5_.frequency = oscState5_.mainFrequencyPlusMatrix;
            float *osc5Values = currentTimbre->osc5_.getNextBlockWithFeedbackAndEnveloppe(&oscState5_, feedbackModulation, env5Value, env5Inc,
                                                                                          oscState5_.frequency, fdbLastValue);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlock(&oscState6_);

            float div3TimesVelocity = .33f * this->velocity;

            float f4x;
            float f4xm1 = freqAi;
            float freq4 = freqAo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq2 = osc2Values[k] * env2Value;
                freq2 *= oscState2_.frequency;

                float freq5 = osc5Values[k];

                oscState4_.frequency = freq5 * voiceIm3 + oscState4_.mainFrequencyPlusMatrix;
                f4x = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * oscState4_.frequency;
                freq4 = f4x - f4xm1 + dcPole * freq4;
                f4xm1 = f4x;

                oscState1_.frequency = freq2 * voiceIm1 + oscState1_.mainFrequencyPlusMatrix;
                float car1 = currentTimbre->osc1_.getNextSample(&oscState1_) * env1Value * mix1 * div3TimesVelocity;

                oscState3_.frequency = freq4 * voiceIm2 + oscState3_.mainFrequencyPlusMatrix;
                float car3 = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * mix2 * div3TimesVelocity;

                float car6 = osc6Values[k] * env6Value * mix3 * div3TimesVelocity;

                *sample++ = car1 * pan1Right + car3 * pan2Right + car6 * pan3Right;
                *sample++ = car1 * pan1Left + car3 * pan2Left + car6 * pan3Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env6Value += env6Inc;
            }

            freqAi = f4xm1;
            freqAo = freq4;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env3_.isDead(&envState3_) && currentTimbre->env6_.isDead(&envState6_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
        }
        break;
    case ALG25:
        /*
         * DX ALGO 29
         *

         .---.  .---.
         | 4 |  | 6*|
         '---'  '---'
         |IM1   |IM2
         .---.  .---.  .---.  .---.
         | 1 |  | 2 |  | 3 |  | 5 |
         '---'  '---'  '---'  '---'
         |Mix1  |Mix2  |Mix3  |Mix4

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState1_.frequency = oscState1_.mainFrequencyPlusMatrix;
            float *osc1Values = currentTimbre->osc1_.getNextBlock(&oscState1_);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlock(&oscState4_);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlockWithFeedbackAndEnveloppe(&oscState6_, feedbackModulation, env6Value, env6Inc,
                                                                                          oscState6_.frequency, fdbLastValue);

            float div4TimesVelocity = .25f * this->velocity;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq6 = osc6Values[k];

                float freq4 = osc4Values[k] * env4Value;
                freq4 *= oscState4_.frequency;

                float car1 = osc1Values[k] * env1Value * mix1 * div4TimesVelocity;

                float car2 = osc2Values[k] * env2Value * mix2 * div4TimesVelocity;

                oscState3_.frequency = freq4 * voiceIm1 + oscState3_.mainFrequencyPlusMatrix;
                float car3 = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * mix3 * div4TimesVelocity;

                oscState5_.frequency = freq6 * voiceIm2 + oscState5_.mainFrequencyPlusMatrix;
                float car5 = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * mix4 * div4TimesVelocity;

                *sample++ = car1 * pan1Right + car2 * pan2Right + car3 * pan3Right + car5 * pan4Right;
                *sample++ = car1 * pan1Left + car2 * pan2Left + car3 * pan3Left + car5 * pan4Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
            }
            if (unlikely(
                    currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env2_.isDead(&envState2_) && currentTimbre->env3_.isDead(&envState3_) && currentTimbre->env5_.isDead(&envState5_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
        }
        break;
    case ALG26:
        /*
         * DX Algo 30
         .---.
         | 5*|
         '---'
         |IM2
         .---.
         | 4 |
         '---'
         |IM1
         .---.  .---.  .---.  .---.
         | 1 |  | 2 |  | 3 |  | 6 |
         '---'  '---'  '---'  '---'
         |Mix1  |Mix2  |Mix3  |Mix4


         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState1_.frequency = oscState1_.mainFrequencyPlusMatrix;
            float *osc1Values = currentTimbre->osc1_.getNextBlock(&oscState1_);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState5_.frequency = oscState5_.mainFrequencyPlusMatrix;
            float *osc5Values = currentTimbre->osc5_.getNextBlockWithFeedbackAndEnveloppe(&oscState5_, feedbackModulation, env5Value, env5Inc,
                                                                                          oscState5_.frequency, fdbLastValue);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlock(&oscState6_);

            float div4TimesVelocity = .25f * this->velocity;

            float f4x;
            float f4xm1 = freqAi;
            float freq4 = freqAo;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float freq5 = osc5Values[k];

                float car1 = osc1Values[k] * env1Value * mix1 * div4TimesVelocity;

                float car2 = osc2Values[k] * env2Value * mix2 * div4TimesVelocity;

                oscState4_.frequency = freq5 * voiceIm2 + oscState4_.mainFrequencyPlusMatrix;
                f4x = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * oscState4_.frequency;
                freq4 = f4x - f4xm1 + dcPole * freq4;
                f4xm1 = f4x;

                oscState3_.frequency = freq4 * voiceIm1 + oscState3_.mainFrequencyPlusMatrix;
                float car3 = currentTimbre->osc3_.getNextSample(&oscState3_) * env3Value * mix3 * div4TimesVelocity;

                float car6 = osc6Values[k] * env6Value * mix4 * div4TimesVelocity;

                *sample++ = car1 * pan1Right + car2 * pan2Right + car3 * pan3Right + car6 * pan4Right;
                *sample++ = car1 * pan1Left + car2 * pan2Left + car3 * pan3Left + car6 * pan4Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env6Value += env6Inc;
            }

            freqAi = f4xm1;
            freqAo = freq4;

            if (unlikely(
                    currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env2_.isDead(&envState2_) && currentTimbre->env3_.isDead(&envState3_) && currentTimbre->env6_.isDead(&envState6_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
        }
        break;

    case ALG27:
        /*
         * DX ALGO 32
         *
         .---.  .---.  .---.  .---.   .---.   .---.
         | 1 |  | 2 |  | 3 |  | 4 |   | 5 |   | 6*|
         '---'  '---'  '---'  '---'   '---'   '---'
         |Mix1  |Mix2  |Mix3  |Mix4   |Mix5   |Mix6

         */
        {
            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState1_.frequency = oscState1_.mainFrequencyPlusMatrix;
            float *osc1Values = currentTimbre->osc1_.getNextBlock(&oscState1_);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlock(&oscState3_);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlock(&oscState4_);

            // Not enough float arrays for a fifth optimization
            oscState5_.frequency = oscState5_.mainFrequencyPlusMatrix;
            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;

            // Multiplier is 1 because op6 is a carrier
            float *osc6Values = currentTimbre->osc6_.getNextBlockWithFeedbackAndEnveloppe(&oscState6_, feedbackModulation, env6Value, env6Inc, 1.0f,
                                                                                          fdbLastValue);

            float div6TimesVelocity = .16f * this->velocity;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float car1 = osc1Values[k] * env1Value * mix1;
                float car2 = osc2Values[k] * env2Value * mix2;
                float car3 = osc3Values[k] * env3Value * mix3;
                float car4 = osc4Values[k] * env4Value * mix4;
                // getNextBlockWithFeedbackAndEnveloppe already applied envelope
                float car6 = osc6Values[k] * mix6;

                float car5 = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * mix5;

                *sample++ = (car1 * pan1Right + car2 * pan2Right + car3 * pan3Right + car4 * pan4Right + car5 * pan5Right + car6 * pan6Right) * div6TimesVelocity;
                *sample++ = (car1 * pan1Left + car2 * pan2Left + car3 * pan3Left + car4 * pan4Left + car5 * pan5Left + car6 * pan6Left) * div6TimesVelocity;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
            }
            if (unlikely(
                    currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env2_.isDead(&envState2_) && currentTimbre->env3_.isDead(&envState3_) && currentTimbre->env4_.isDead(&envState4_) && currentTimbre->env5_.isDead(&envState5_) && currentTimbre->env6_.isDead(&envState6_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
        }
        break;
    case ALG28:
        /*
         * DX ALGO 31

         .---.
         | 6*|
         '---'
         |IM1
         .---.  .---.  .---.  .---.   .---.
         | 1 |  | 2 |  | 3 |  | 4 |   | 5 |
         '---'  '---'  '---'  '---'   '---'
         |Mix1  |Mix2  |Mix3  |Mix4   |Mix5

         */
        {
            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);
            currentTimbre->osc5_.calculateFrequencyWithMatrix(&oscState5_, &matrix, freqHarm);
            currentTimbre->osc6_.calculateFrequencyWithMatrix(&oscState6_, &matrix, freqHarm);

            float voiceIm1 = modulationIndex1;

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            prepareEnvBlock(currentTimbre->env5_, envState5_, this->env5ValueMem, env5Value, env5Inc, 4, osMode, invStep);

            prepareEnvBlock(currentTimbre->env6_, envState6_, this->env6ValueMem, env6Value, env6Inc, 5, osMode, invStep);

            oscState1_.frequency = oscState1_.mainFrequencyPlusMatrix;
            float *osc1Values = currentTimbre->osc1_.getNextBlock(&oscState1_);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlock(&oscState3_);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlock(&oscState4_);

            oscState6_.frequency = oscState6_.mainFrequencyPlusMatrix;
            float *osc6Values = currentTimbre->osc6_.getNextBlockWithFeedbackAndEnveloppe(&oscState6_, feedbackModulation, env6Value, env6Inc,
                                                                                          oscState6_.frequency, fdbLastValue);

            float div5TimesVelocity = .2f * this->velocity;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                float car1 = osc1Values[k] * env1Value * mix1;
                float car2 = osc2Values[k] * env2Value * mix2;
                float car3 = osc3Values[k] * env3Value * mix3;
                float car4 = osc4Values[k] * env4Value * mix4;

                float freq6 = osc6Values[k];

                oscState5_.frequency = freq6 * voiceIm1 + oscState5_.mainFrequencyPlusMatrix;
                float car5 = currentTimbre->osc5_.getNextSample(&oscState5_) * env5Value * mix5;

                *sample++ = (car1 * pan1Right + car2 * pan2Right + car3 * pan3Right + car4 * pan4Right + car5 * pan5Right) * div5TimesVelocity;
                *sample++ = (car1 * pan1Left + car2 * pan2Left + car3 * pan3Left + car4 * pan4Left + car5 * pan5Left) * div5TimesVelocity;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
                env5Value += env5Inc;
            }
            if (unlikely(
                    currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env2_.isDead(&envState2_) && currentTimbre->env3_.isDead(&envState3_) && currentTimbre->env4_.isDead(&envState4_) && currentTimbre->env5_.isDead(&envState5_) && currentTimbre->env6_.isDead(&envState6_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
        }
        break;

    case ALG29:
        /* Windowed sync AM, 1 synced by 3
        inspiration from https://electro-music.com/nm_classic/015_workshops/Clavia/NordModularWorkshops&Threads/WerkMap/WorkShops/Hordijk1999-2000/VOSIM.html
                         https://scsynth.org/t/fake-resonance-aka-windowed-sync/2720

                 IM3
                <----
         .---.         .---.
         | 3*|         | 4 |
         '---'         '---'
         |IM1           |IM2
         .---.          .---.
         |*1 |          | 2 |
         '---'          '---'
         |Mix1          |Mix2

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlockWithFeedbackAndEnveloppe(&oscState4_, feedbackModulation, env4Value, env4Inc,
                                                                                          oscState4_.frequency, fdbLastValue);

            float div2TimesVelocity = this->velocity * .5f;

            float mix1V = mix1 * div2TimesVelocity;
            float mix2V = mix2 * div2TimesVelocity;

            float osc1Matrix = oscState1_.mainFrequencyPlusMatrix - oscState1_.mainFrequency;
            float osc1FrequencyPlusMatrix = oscState1_.mainFrequency + osc1Matrix * fabsf(osc1Matrix) * 0.25f;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {
                float freq4 = osc4Values[k];

                oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix + voiceIm3 * freq4;
                float osc3 = currentTimbre->osc3_.getNextSample(&oscState3_);
                float phase = currentTimbre->osc3_.getPhase(&oscState3_);
                float freq3 = osc3 * env3Value * oscState3_.frequency;

                // sync check
                bool isSync = prevPhase > phase;
                prevPhase = phase;
                float window = 1.f - phase;

                oscState1_.frequency = freq3 * voiceIm1 + osc1FrequencyPlusMatrix;

                // sync slave osc (sub-sample accurate)
                if (unlikely(isSync))
                {
                    currentTimbre->osc1_.hardSyncToMaster(&oscState1_, &currentTimbre->osc3_, &oscState3_, phase);
                }

                float carSample1 = currentTimbre->osc1_.getNextSample(&oscState1_) * window;
                carSample1 *= fabsf(carSample1); // vosim shape
                float car1 = carSample1 * env1Value * mix1V;

                oscState2_.frequency = freq4 * voiceIm2 + oscState2_.mainFrequencyPlusMatrix;
                float car2 = currentTimbre->osc2_.getNextSample(&oscState2_) * env2Value * mix2V;

                *sample++ = car1 * pan1Right + car2 * pan2Right;
                *sample++ = car1 * pan1Left + car2 * pan2Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
            }

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env2_.isDead(&envState2_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
            break;
        }
    case ALG30:
        /* Windowed sync AM, 1 & 2 synced by 3

         .---.         .---.
         | 3*|         | 4 |
         '---'         '---'
         |IM1 \IM2 IM3/ |IM4
         .---.          .---.
         |*1 |          |*2 |
         '---'          '---'
         |Mix1          |Mix2

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
            float *osc4Values = currentTimbre->osc4_.getNextBlockWithFeedbackAndEnveloppe(&oscState4_, feedbackModulation, env4Value, env4Inc,
                                                                                          oscState4_.frequency, fdbLastValue);

            float div2TimesVelocity = this->velocity * .5f;

            float mix1V = mix1 * div2TimesVelocity;
            float mix2V = mix2 * div2TimesVelocity;

            float osc1Matrix = oscState1_.mainFrequencyPlusMatrix - oscState1_.mainFrequency;
            float osc1FrequencyPlusMatrix = oscState1_.mainFrequency + osc1Matrix * fabsf(osc1Matrix) * 0.25f;

            float osc2Matrix = oscState2_.mainFrequencyPlusMatrix - oscState2_.mainFrequency;
            float osc2FrequencyPlusMatrix = oscState2_.mainFrequency + osc2Matrix * fabsf(osc2Matrix) * 0.25f;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {
                float freq4 = osc4Values[k];

                oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
                float osc3 = currentTimbre->osc3_.getNextSample(&oscState3_);
                float phase = currentTimbre->osc3_.getPhase(&oscState3_);
                float freq3 = osc3 * env3Value * oscState3_.frequency;

                // sync check
                bool isSync = prevPhase > phase;
                prevPhase = phase;
                float window = 1.f - phase;

                oscState1_.frequency = freq3 * voiceIm1 + osc1FrequencyPlusMatrix + freq4 * voiceIm3;
                oscState2_.frequency = freq3 * voiceIm2 + osc2FrequencyPlusMatrix + freq4 * voiceIm4;

                // sync slave osc (sub-sample accurate)
                if (unlikely(isSync))
                {
                    currentTimbre->osc1_.hardSyncToMaster(&oscState1_, &currentTimbre->osc3_, &oscState3_, phase);
                    currentTimbre->osc2_.hardSyncToMaster(&oscState2_, &currentTimbre->osc3_, &oscState3_, phase);
                }

                float carSample1 = currentTimbre->osc1_.getNextSample(&oscState1_) * window;
                carSample1 *= fabsf(carSample1); // vosim shape
                float car1 = carSample1 * env1Value * mix1V;

                float carSample2 = currentTimbre->osc2_.getNextSample(&oscState2_) * window;
                carSample2 *= fabsf(carSample2); // vosim shape
                float car2 = carSample2 * env2Value * mix2V;

                *sample++ = car1 * pan1Right + car2 * pan2Right;
                *sample++ = car1 * pan1Left + car2 * pan2Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
            }

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env2_.isDead(&envState2_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }
            break;
        }
    case ALG31:
        /* Windowed sync AM, 1 2 3 synced by 4

         .---.
         | 4*|
         '---'
         /IM1 |IM2 \IM3
         .---.  .---.  .---.
         |*1 |  |*2 |  |*3 |
         '---'  '---'  '---'
         |Mix1  |Mix2  | Mix3

         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            float div3TimesVelocity = .33f * this->velocity;

            float mix1V = mix1 * div3TimesVelocity;
            float mix2V = mix2 * div3TimesVelocity;
            float mix3V = mix3 * div3TimesVelocity;

            float osc1Matrix = oscState1_.mainFrequencyPlusMatrix - oscState1_.mainFrequency;
            float osc1FrequencyPlusMatrix = oscState1_.mainFrequency + osc1Matrix * fabsf(osc1Matrix) * 0.25f;

            float osc2Matrix = oscState2_.mainFrequencyPlusMatrix - oscState2_.mainFrequency;
            float osc2FrequencyPlusMatrix = oscState2_.mainFrequency + osc2Matrix * fabsf(osc2Matrix) * 0.25f;

            float osc3Matrix = oscState3_.mainFrequencyPlusMatrix - oscState3_.mainFrequency;
            float osc3FrequencyPlusMatrix = oscState3_.mainFrequency + osc3Matrix * fabsf(osc3Matrix) * 0.25f;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {

                oscState4_.frequency = oscState4_.mainFrequencyPlusMatrix;
                float osc4 = currentTimbre->osc4_.getNextSample(&oscState4_);
                float phase = currentTimbre->osc4_.getPhase(&oscState4_);
                float freq4 = osc4 * env4Value * oscState4_.frequency;

                bool isSync = prevPhase > phase;
                prevPhase = phase;
                float window = 1.f - phase;

                oscState3_.frequency = freq4 * voiceIm3 + osc3FrequencyPlusMatrix;
                oscState2_.frequency = freq4 * voiceIm2 + osc2FrequencyPlusMatrix;
                oscState1_.frequency = freq4 * voiceIm1 + osc1FrequencyPlusMatrix;

                // sync slave osc (sub-sample accurate)
                if (unlikely(isSync))
                {
                    currentTimbre->osc1_.hardSyncToMaster(&oscState1_, &currentTimbre->osc4_, &oscState4_, phase);
                    currentTimbre->osc2_.hardSyncToMaster(&oscState2_, &currentTimbre->osc4_, &oscState4_, phase);
                    currentTimbre->osc3_.hardSyncToMaster(&oscState3_, &currentTimbre->osc4_, &oscState4_, phase);
                }

                float carSample3 = currentTimbre->osc3_.getNextSample(&oscState3_) * window;
                carSample3 *= fabsf(carSample3); // vosim shape
                float car3 = carSample3 * env3Value * mix3V;

                float carSample2 = currentTimbre->osc2_.getNextSample(&oscState2_) * window;
                carSample2 *= fabsf(carSample2); // vosim shape
                float car2 = carSample2 * env2Value * mix2V;

                float carSample1 = currentTimbre->osc1_.getNextSample(&oscState1_) * window;
                carSample1 *= fabsf(carSample1); // vosim shape
                float car1 = carSample1 * env1Value * mix1V;

                *sample++ = car1 * pan1Right + car2 * pan2Right + car3 * pan3Right;
                *sample++ = car1 * pan1Left + car2 * pan2Left + car3 * pan3Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env3Value += env3Inc;
                env4Value += env4Inc;
            }

            if (unlikely(currentTimbre->env1_.isDead(&envState1_) && currentTimbre->env2_.isDead(&envState2_) && currentTimbre->env3_.isDead(&envState3_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }

    case ALG32:
        /* Windowed sync AM, 1 synced by 4

                      IM4
                     ---->
         .---.  .---.     .---.
         | 2 |  | 3 |     | 4*|
         '---'  '---'     '---'
         \IM1  |IM2    /IM3
         .---.
         |*1 |
         '---'
         */
        {
            float voiceIm1 = modulationIndex1;
            float voiceIm2 = modulationIndex2;
            float voiceIm3 = modulationIndex3;
            float voiceIm4 = modulationIndex4;

            currentTimbre->osc1_.calculateFrequencyWithMatrix(&oscState1_, &matrix, freqHarm);
            currentTimbre->osc2_.calculateFrequencyWithMatrix(&oscState2_, &matrix, freqHarm);
            currentTimbre->osc3_.calculateFrequencyWithMatrix(&oscState3_, &matrix, freqHarm);
            currentTimbre->osc4_.calculateFrequencyWithMatrix(&oscState4_, &matrix, freqHarm);

            prepareEnvBlock(currentTimbre->env1_, envState1_, this->env1ValueMem, env1Value, env1Inc, 0, osMode, invStep);

            prepareEnvBlock(currentTimbre->env2_, envState2_, this->env2ValueMem, env2Value, env2Inc, 1, osMode, invStep);

            prepareEnvBlock(currentTimbre->env3_, envState3_, this->env3ValueMem, env3Value, env3Inc, 2, osMode, invStep);

            prepareEnvBlock(currentTimbre->env4_, envState4_, this->env4ValueMem, env4Value, env4Inc, 3, osMode, invStep);

            oscState2_.frequency = oscState2_.mainFrequencyPlusMatrix;
            float *osc2Values = currentTimbre->osc2_.getNextBlock(&oscState2_);

            oscState3_.frequency = oscState3_.mainFrequencyPlusMatrix;
            float *osc3Values = currentTimbre->osc3_.getNextBlockWithFeedbackAndEnveloppe(&oscState3_, feedbackModulation, env3Value, env3Inc,
                                                                                          oscState3_.frequency, fdbLastValue);

            float f4x;
            float f4xm1 = freqAi;
            float freq4 = freqAo;
            float mix1V = mix1 * this->velocity;

            float osc1Matrix = oscState1_.mainFrequencyPlusMatrix - oscState1_.mainFrequency;
            float osc1FrequencyPlusMatrix = oscState1_.mainFrequency + osc1Matrix * fabsf(osc1Matrix) * 0.25f;

            for (int k = 0; k < BLOCK_SIZE; k++)
            {
                float freq2 = osc2Values[k] * env2Value * oscState2_.frequency;
                float freq3 = osc3Values[k];

                oscState4_.frequency = freq3 * voiceIm4 + oscState4_.mainFrequencyPlusMatrix;

                f4x = currentTimbre->osc4_.getNextSample(&oscState4_) * env4Value * oscState4_.frequency;
                float phase = currentTimbre->osc4_.getPhase(&oscState4_);
                freq4 = f4x - f4xm1 + dcPole * freq4;
                f4xm1 = f4x;

                bool isSync = prevPhase > phase;
                prevPhase = phase;
                float window = 1.f - phase;

                oscState1_.frequency = freq2 * voiceIm1 + freq3 * voiceIm2 + freq4 * voiceIm3 + osc1FrequencyPlusMatrix;

                // sync slave osc (sub-sample accurate)
                if (unlikely(isSync))
                {
                    currentTimbre->osc1_.hardSyncToMaster(&oscState1_, &currentTimbre->osc4_, &oscState4_, phase);
                }

                float carSample1 = currentTimbre->osc1_.getNextSample(&oscState1_) * window;
                carSample1 *= fabsf(carSample1); // vosim shape
                float currentSample = carSample1 * env1Value * mix1V;

                *sample++ = currentSample * pan1Right;
                *sample++ = currentSample * pan1Left;

                env1Value += env1Inc;
                env2Value += env2Inc;
                env4Value += env4Inc;
            }

            freqAi = f4xm1;
            freqAo = freq4;

            if (unlikely(currentTimbre->env1_.isDead(&envState1_)))
            {
                if (allowEndNote)
                {
                    endNoteOrBeginNextOne();
                }
            }

            break;
        }

    } // End switch
}

void Voice::setCurrentTimbre(Timbre *timbre)
{
    if (timbre == 0)
    {
        this->currentTimbre = 0;
        return;
    }

    struct LfoParams *lfoParams[] = {&timbre->getParamRaw()->lfoOsc1, &timbre->getParamRaw()->lfoOsc2, &timbre->getParamRaw()->lfoOsc3};
    float *lfoSyncModes = &timbre->getParamRaw()->lfoSyncModes.lfo1;
    struct StepSequencerParams *stepseqparams[] = {&timbre->getParamRaw()->lfoSeq1, &timbre->getParamRaw()->lfoSeq2};
    struct StepSequencerSteps *stepseqs[] = {&timbre->getParamRaw()->lfoSteps1, &timbre->getParamRaw()->lfoSteps2};

    this->currentTimbre = timbre;

    matrix.init(&timbre->getParamRaw()->matrixRowState1);

    // OSC
    for (int k = 0; k < NUMBER_OF_LFO_OSC; k++)
    {
        float *phase = &((float *)&timbre->getParamRaw()->lfoPhases.phaseLfo1)[k];
        lfoOsc[k].init(lfoParams[k], &lfoSyncModes[k], phase, &this->matrix, (SourceEnum)(MATRIX_SOURCE_LFO1 + k), (DestinationEnum)(LFO1_FREQ + k));
    }

    // ENV
    struct LfoEnvCurveParams *lfoEnvCurves = &timbre->getParamRaw()->lfoEnvCurves;
    lfoEnv[0].init(&timbre->getParamRaw()->lfoEnv1, &lfoEnvCurves->env1AtkCurve, &lfoEnvCurves->env1RelCurve, &this->matrix, MATRIX_SOURCE_LFOENV1, (DestinationEnum)0);
    lfoEnv2[0].init(&timbre->getParamRaw()->lfoEnv2, &lfoEnvCurves->env2AtkCurve, &lfoEnvCurves->env2RelCurve, &this->matrix, MATRIX_SOURCE_LFOENV2, (DestinationEnum)LFOENV2_SILENCE);

    // Step sequencer
    for (int k = 0; k < NUMBER_OF_LFO_STEP; k++)
    {
        lfoStepSeq[k].init(stepseqparams[k], stepseqs[k], &matrix, (SourceEnum)(MATRIX_SOURCE_LFOSEQ1 + k), (DestinationEnum)(LFOSEQ1_GATE + k), (DestinationEnum)(LFOSEQ1_LENGTH + k));
    }

    setNewEffectParam(0);

    // init pan
    const int *opInfo = algoOpInformation[(int)timbre->getParamRaw()->engine1.algo];
    if (likely(opInfo[0] == 1))
    {
        float pan1 = timbre->getParamRaw()->engineMix1.panOsc1 + 1.0f;
        int pan = __USAT((int)(pan1 * 128), 8);
        pan1Left = panTable[pan];
        pan1Right = panTable[256 - pan];
    }
    if (likely(opInfo[1] == 1))
    {
        float pan1 = timbre->getParamRaw()->engineMix1.panOsc2 + 1.0f;
        int pan = __USAT((int)(pan1 * 128), 8);
        pan2Left = panTable[pan];
        pan2Right = panTable[256 - pan];
    }
    if (likely(opInfo[2] == 1))
    {
        float pan1 = timbre->getParamRaw()->engineMix2.panOsc3 + 1.0f;
        int pan = __USAT((int)(pan1 * 128), 8);
        pan3Left = panTable[pan];
        pan3Right = panTable[256 - pan];
    }
    if (likely(opInfo[3] == 1))
    {
        float pan1 = timbre->getParamRaw()->engineMix2.panOsc4 + 1.0f;
        int pan = __USAT((int)(pan1 * 128), 8);
        pan4Left = panTable[pan];
        pan4Right = panTable[256 - pan];
    }
    if (likely(opInfo[4] == 1))
    {
        float pan1 = timbre->getParamRaw()->engineMix3.panOsc5 + 1.0f;
        int pan = __USAT((int)(pan1 * 128), 8);
        pan5Left = panTable[pan];
        pan5Right = panTable[256 - pan];
    }
    if (likely(opInfo[5] == 1))
    {
        float pan1 = timbre->getParamRaw()->engineMix3.panOsc6 + 1.0f;
        int pan = __USAT((int)(pan1 * 128), 8);
        pan6Left = panTable[pan];
        pan6Right = panTable[256 - pan];
    }

    // init mixer gain
    mixerGain = currentTimbre->params_.effect1.param3;
}

void Voice::midiClockSongPositionStep(int songPosition, bool recomputeNext)
{
    if (likely(currentTimbre->isLfoUsed(0)))
    {
        lfoOsc[0].midiClock(songPosition, recomputeNext);
    }
    if (likely(currentTimbre->isLfoUsed(1)))
    {
        lfoOsc[1].midiClock(songPosition, recomputeNext);
    }
    if (likely(currentTimbre->isLfoUsed(2)))
    {
        lfoOsc[2].midiClock(songPosition, recomputeNext);
    }
    if (likely(currentTimbre->isLfoUsed(3)))
    {
        lfoEnv[0].midiClock(songPosition, recomputeNext);
    }
    if (likely(currentTimbre->isLfoUsed(4)))
    {
        lfoEnv2[0].midiClock(songPosition, recomputeNext);
    }
    if (likely(currentTimbre->isLfoUsed(5)))
    {
        lfoStepSeq[0].midiClock(songPosition, recomputeNext);
    }
    if (likely(currentTimbre->isLfoUsed(6)))
    {
        lfoStepSeq[1].midiClock(songPosition, recomputeNext);
    }
}

void Voice::midiClockContinue(int songPosition)
{
    if (likely(currentTimbre->isLfoUsed(0)))
    {
        lfoOsc[0].midiClock(songPosition, false);
    }
    if (likely(currentTimbre->isLfoUsed(1)))
    {
        lfoOsc[1].midiClock(songPosition, false);
    }
    if (likely(currentTimbre->isLfoUsed(2)))
    {
        lfoOsc[2].midiClock(songPosition, false);
    }
    if (likely(currentTimbre->isLfoUsed(3)))
    {
        lfoEnv[0].midiClock(songPosition, false);
    }
    if (likely(currentTimbre->isLfoUsed(4)))
    {
        lfoEnv2[0].midiClock(songPosition, false);
    }
    if (likely(currentTimbre->isLfoUsed(5)))
    {
        lfoStepSeq[0].midiClock(songPosition, false);
    }
    if (likely(currentTimbre->isLfoUsed(6)))
    {
        lfoStepSeq[1].midiClock(songPosition, false);
    }
}

void Voice::midiClockStart()
{
    if (likely(currentTimbre->isLfoUsed(0)))
    {
        lfoOsc[0].midiContinue();
    }
    if (likely(currentTimbre->isLfoUsed(1)))
    {
        lfoOsc[1].midiContinue();
    }
    if (likely(currentTimbre->isLfoUsed(2)))
    {
        lfoOsc[2].midiContinue();
    }
    if (likely(currentTimbre->isLfoUsed(3)))
    {
        lfoEnv[0].midiContinue();
    }
    if (likely(currentTimbre->isLfoUsed(4)))
    {
        lfoEnv2[0].midiContinue();
    }
    if (likely(currentTimbre->isLfoUsed(5)))
    {
        lfoStepSeq[0].midiContinue();
    }
    if (likely(currentTimbre->isLfoUsed(6)))
    {
        lfoStepSeq[1].midiContinue();
    }
}

#if defined(__GNUC__) && !defined(DEBUG)
#pragma GCC pop_options
#endif
