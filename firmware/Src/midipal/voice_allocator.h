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
    VoiceAllocator() : pool_(0), lru_(0), size_(0), cyclic_allocator_(0xff) {
    }

    void Init(uint8_t size, bool cyclicMode, uint8_t* pool, uint8_t* lru) {
        size_ = size;
        pool_ = pool;
        lru_ = lru;
        cyclic_allocator_ = cyclicMode ? 0 : 0xff;
        Clear();
    }

    uint8_t NoteOn(uint8_t note, const uint8_t* unavailable = 0, bool restored = false) {
        if (size_ == 0) {
            return 0xff;
        }

        uint8_t voice;
        if (cyclic_allocator_ == 0xff) {
            voice = Find(note);
            if (voice != 0xff && unavailable && unavailable[voice]) {
                return 0xff;
            }

            if (voice == 0xff) {
                for (uint8_t i = 0; i < size_; ++i) {
                    if (lru_[i] < size_ && !(pool_[lru_[i]] & 0x80)
                            && (!unavailable || !unavailable[lru_[i]])) {
                        voice = lru_[i];
                    }
                }
            }
            if (voice == 0xff) {
                // Bass-protect voice stealing:
                // Find the lowest active note among available candidate voices
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

                // If multiple voices are active, protect the lowest active note (bass)
                // and steal the oldest among the higher/chord notes!
                if (active_candidates > 1 && lowest_voice != 0xff) {
                    for (uint8_t i = 0; i < size_; ++i) {
                        uint8_t v = lru_[i];
                        if (v < size_ && v != lowest_voice && (!unavailable || !unavailable[v])) {
                            voice = v;
                        }
                    }
                }

                // Fallback: if only 1 voice or all other voices unavailable, steal oldest
                if (voice == 0xff) {
                    for (uint8_t i = 0; i < size_; ++i) {
                        if (lru_[i] < size_ && (!unavailable || !unavailable[lru_[i]])) {
                            voice = lru_[i];
                        }
                    }
                }
            }
        } else {
            uint8_t attempts = 0;
            do {
                ++cyclic_allocator_;
                if (cyclic_allocator_ >= size_) {
                    cyclic_allocator_ = 0;
                }
                ++attempts;
            } while (unavailable && unavailable[cyclic_allocator_] && attempts < size_);
            if (unavailable && unavailable[cyclic_allocator_]) {
                return 0xff;
            }
            voice = cyclic_allocator_;
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
        if (cyclic_allocator_ == 0xff) {
            if (voice != 0xff) {
                pool_[voice] &= 0x7f;
                Touch(voice);
            }
        } else if (voice != 0xff) {
            pool_[voice] = 0xff;
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
};

}  // namespace ambika

#endif  // MIDIPAL_VOICE_ALLOCATOR_H_
