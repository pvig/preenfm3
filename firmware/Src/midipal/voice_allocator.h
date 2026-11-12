// Copyright 2011 Emilie Gillet.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// Polyphonic voice allocator, adapted from the Ambika firmware.

#ifndef MIDIPAL_VOICE_ALLOCATOR_H_
#define MIDIPAL_VOICE_ALLOCATOR_H_

#include <stdint.h>

namespace ambika {

class VoiceAllocator {
public:
    enum PriorityPolicy {
        PRIORITY_LAST_KEY = 0,
        PRIORITY_HIGH_KEY,
        PRIORITY_LOW_KEY,
        PRIORITY_BASS_PROTECT,
        PRIORITY_CYCLIC
    };

    VoiceAllocator()
        : pool_(0), lru_(0), size_(0), cyclic_allocator_(0xff), priority_mode_(PRIORITY_LAST_KEY) {
    }

    void Init(uint8_t size, uint8_t priorityMode, uint8_t* pool, uint8_t* lru) {
        size_ = size;
        pool_ = pool;
        lru_ = lru;
        Clear();
        SetPriorityMode(priorityMode);
    }

    void SetPriorityMode(uint8_t priorityMode) {
        priority_mode_ = priorityMode;
        if (priorityMode == PRIORITY_CYCLIC) {
            if (cyclic_allocator_ == 0xff) {
                cyclic_allocator_ = 0;
            }
        } else {
            cyclic_allocator_ = 0xff;
        }
    }

    uint8_t NoteOn(uint8_t note, const uint8_t* unavailable = 0, bool restored = false) {
        if (size_ == 0) {
            return 0xff;
        }

        uint8_t voice = 0xff;
        if (cyclic_allocator_ == 0xff) {
            voice = Find(note);
            if (voice != 0xff && unavailable && unavailable[voice]) {
                return 0xff;
            }

            if (voice == 0xff) {
                voice = FindOldestFree(unavailable);
            }
            if (voice == 0xff) {
                voice = SelectStealVictim(unavailable);
            }
        } else {
            voice = SelectCyclicVoice(unavailable);
        }

        if (voice == 0xff) {
            return 0xff;
        }
        pool_[voice] = 0x80 | note;
        if (restored) {
            SendToBack(voice);
        } else {
            Touch(voice);
        }
        return voice;
    }

    uint8_t Find(uint8_t note) const {
        for (uint8_t i = 0; i < size_; ++i) {
            if ((pool_[i] & 0x7f) == note) {
                return i;
            }
        }
        return 0xff;
    }

    uint8_t FindActive(uint8_t note) const {
        for (uint8_t i = 0; i < size_; ++i) {
            if ((pool_[i] & 0x7f) == note && (pool_[i] & 0x80)) {
                return i;
            }
        }
        return 0xff;
    }

    uint8_t NoteOff(uint8_t note) {
        uint8_t voice = FindActive(note);
        if (voice != 0xff) {
            pool_[voice] &= 0x7f;
            Touch(voice);
        }
        return voice;
    }

    void Clear() {
        ClearNotes();
        for (uint8_t i = 0; i < size_; ++i) {
            lru_[i] = size_ - i - 1;
        }
    }

    void ClearNotes() {
        for (uint8_t i = 0; i < size_; ++i) {
            pool_[i] = 0;
        }
    }

private:
    uint8_t FindOldestFree(const uint8_t* unavailable) const {
        uint8_t voice = 0xff;
        for (uint8_t i = 0; i < size_; ++i) {
            uint8_t v = lru_[i];
            if (v < size_ && !(pool_[v] & 0x80) && (!unavailable || !unavailable[v])) {
                voice = v;
            }
        }
        return voice;
    }

    uint8_t FindOldestActive(const uint8_t* unavailable, uint8_t skipVoice = 0xff) const {
        uint8_t voice = 0xff;
        for (uint8_t i = 0; i < size_; ++i) {
            uint8_t v = lru_[i];
            if (v < size_ && v != skipVoice && (!unavailable || !unavailable[v])) {
                voice = v;
            }
        }
        return voice;
    }

    uint8_t FindExtremePitchVoice(const uint8_t* unavailable, bool findHighest) const {
        uint8_t extreme_note = findHighest ? 0 : 0x7f;
        uint8_t extreme_voice = 0xff;
        for (uint8_t i = 0; i < size_; ++i) {
            uint8_t v = lru_[i];
            if (v < size_ && (!unavailable || !unavailable[v])) {
                uint8_t note_played = pool_[v] & 0x7f;
                // <= / >= picks the oldest voice when multiple voices play the same extreme pitch.
                bool isMoreExtreme = findHighest ? (note_played >= extreme_note) : (note_played <= extreme_note);
                if (extreme_voice == 0xff || isMoreExtreme) {
                    extreme_note = note_played;
                    extreme_voice = v;
                }
            }
        }
        return extreme_voice;
    }

    uint8_t SelectStealVictim(const uint8_t* unavailable) const {
        switch (priority_mode_) {
        case PRIORITY_HIGH_KEY:
            return FindExtremePitchVoice(unavailable, false);
        case PRIORITY_LOW_KEY:
            return FindExtremePitchVoice(unavailable, true);
        case PRIORITY_BASS_PROTECT: {
            uint8_t lowest_note = 0x7f;
            uint8_t lowest_voice = 0xff;
            uint8_t active_candidates = 0;
            for (uint8_t i = 0; i < size_; ++i) {
                uint8_t v = lru_[i];
                if (v < size_ && (!unavailable || !unavailable[v])) {
                    active_candidates++;
                    uint8_t note_played = pool_[v] & 0x7f;
                    if (note_played < lowest_note) {
                        lowest_note = note_played;
                        lowest_voice = v;
                    }
                }
            }
            if (active_candidates > 1 && lowest_voice != 0xff) {
                uint8_t victim = FindOldestActive(unavailable, lowest_voice);
                if (victim != 0xff) {
                    return victim;
                }
            }
            return FindOldestActive(unavailable);
        }
        case PRIORITY_LAST_KEY:
        default:
            return FindOldestActive(unavailable);
        }
    }

    uint8_t SelectCyclicVoice(const uint8_t* unavailable) {
        if (cyclic_allocator_ >= size_) {
            cyclic_allocator_ = 0;
        }
        // Prefer an inactive voice in round-robin order before stealing an active note.
        uint8_t candidate = cyclic_allocator_;
        for (uint8_t attempts = 0; attempts < size_; ++attempts) {
            ++candidate;
            if (candidate >= size_) {
                candidate = 0;
            }
            if (!(pool_[candidate] & 0x80) && (!unavailable || !unavailable[candidate])) {
                cyclic_allocator_ = candidate;
                return candidate;
            }
        }
        candidate = cyclic_allocator_;
        for (uint8_t attempts = 0; attempts < size_; ++attempts) {
            ++candidate;
            if (candidate >= size_) {
                candidate = 0;
            }
            if (!unavailable || !unavailable[candidate]) {
                cyclic_allocator_ = candidate;
                return candidate;
            }
        }
        return 0xff;
    }

    void Touch(uint8_t voice) {
        int8_t source = size_ - 1;
        int8_t destination = size_ - 1;
        while (source >= 0) {
            if (lru_[source] != voice) {
                lru_[destination--] = lru_[source];
            }
            --source;
        }
        lru_[0] = voice;
    }

    void SendToBack(uint8_t voice) {
        uint8_t source = 0;
        uint8_t destination = 0;
        while (source < size_) {
            if (lru_[source] != voice) {
                lru_[destination++] = lru_[source];
            }
            ++source;
        }
        lru_[size_ - 1] = voice;
    }

    uint8_t* pool_;
    uint8_t* lru_;
    uint8_t size_;
    uint8_t cyclic_allocator_;
    uint8_t priority_mode_;
};

}  // namespace ambika

#endif  // MIDIPAL_VOICE_ALLOCATOR_H_
