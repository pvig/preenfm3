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

#ifndef MATRIX_H_
#define MATRIX_H_

#define MATRIX_SIZE 12

#include <math.h>
#include "Common.h"

class Matrix {
    friend class Timbre;
public:
    Matrix();
    ~Matrix();

    void init(struct MatrixRowParams *matrixRows);

    void resetSources() {
        for (int k = 0; k < MATRIX_SOURCE_MAX; k++) {
            setSource((SourceEnum) k, 0);
        }
    }

    void resetAllDestination() {
        for (int k = 0; k < DESTINATION_MAX; k++) {
            destinations[k] = 0;
        }
    }

    void resetDestination(int k) {
        destinations[k] = 0;
    }

    void computeAllDestinations() {
        // Store values before erasing destinations
        float mul1 = destinations[MTX1_MUL];
        float mul2 = destinations[MTX2_MUL];
        float mul3 = destinations[MTX3_MUL];
        float mul4 = destinations[MTX4_MUL];

        // Let's erase destinations first
        int k = 0;
        // No need to erase first one (tiny optim)
        // destinations[(int)rows[k].dest1] = 0;
        destinations[(int) rows[k++].dest2] = 0;
        destinations[(int) rows[k].dest1] = 0;
        destinations[(int) rows[k++].dest2] = 0;
        destinations[(int) rows[k].dest1] = 0;
        destinations[(int) rows[k++].dest2] = 0;
        destinations[(int) rows[k].dest1] = 0;
        destinations[(int) rows[k++].dest2] = 0;
        destinations[(int) rows[k].dest1] = 0;
        destinations[(int) rows[k++].dest2] = 0;
        destinations[(int) rows[k].dest1] = 0;
        destinations[(int) rows[k++].dest2] = 0;
        destinations[(int) rows[k].dest1] = 0;
        destinations[(int) rows[k++].dest2] = 0;
        destinations[(int) rows[k].dest1] = 0;
        destinations[(int) rows[k++].dest2] = 0;
        destinations[(int) rows[k].dest1] = 0;
        destinations[(int) rows[k++].dest2] = 0;
        destinations[(int) rows[k].dest1] = 0;
        destinations[(int) rows[k++].dest2] = 0;
        destinations[(int) rows[k].dest1] = 0;
        destinations[(int) rows[k++].dest2] = 0;
        destinations[(int) rows[k].dest1] = 0;
        destinations[(int) rows[k++].dest2] = 0;

        // For first row can be modified by destination
        if (likely(rows[0].source != MATRIX_SOURCE_NONE)) {
            applyRow<true>(0, rows[0].mul + mul1);
        }
        if (likely(rows[1].source != MATRIX_SOURCE_NONE)) {
            applyRow<false>(1, rows[1].mul + mul2);
        }
        if (likely(rows[2].source != MATRIX_SOURCE_NONE)) {
            applyRow<false>(2, rows[2].mul + mul3);
        }
        if (likely(rows[3].source != MATRIX_SOURCE_NONE)) {
            applyRow<false>(3, rows[3].mul + mul4);
        }

        // and compute other rows
        for (int r = 4; r < MATRIX_SIZE; r++) {
            if (likely(rows[r].source != MATRIX_SOURCE_NONE && rows[r].mul != 0.0f)) {
                applyRow<false>(r, rows[r].mul);
            }
        }
    }

    template <bool isFirstRow>
    inline __attribute__((always_inline)) void applyRow(int r, float mul) {
        uint8_t srcIdx = (uint8_t) rows[r].source;
        float src = sources[srcIdx];
        float val1 = src * mul;
        float val2 = val1;
        uint8_t d1 = (uint8_t) rows[r].dest1;
        uint8_t d2 = (uint8_t) rows[r].dest2;

        // StepSeq outputs linear values (step / 12.0f); convert to 2^x - 1 only for linear-Hz Osc Freq destinations.
        if (unlikely(srcIdx == MATRIX_SOURCE_LFOSEQ1 || srcIdx == MATRIX_SOURCE_LFOSEQ2)) {
            bool d1Osc = (d1 >= OSC1_FREQ && d1 <= ALL_OSC_FREQ);
            bool d2Osc = (d2 >= OSC1_FREQ && d2 <= ALL_OSC_FREQ);
            if (unlikely(d1Osc || d2Osc)) {
                float expMul = (exp2f(src) - 1.0f) * mul;
                if (d1Osc) {
                    val1 = expMul;
                }
                if (d2Osc) {
                    val2 = expMul;
                }
            }
        }

        if (isFirstRow) {
            destinations[d1] = val1;
        } else {
            destinations[d1] += val1;
        }
        destinations[d2] += val2;
    }

    void setSource(SourceEnum source, float value) __attribute__((always_inline)) {
        this->sources[source] = value;
    }

    float getDestination(DestinationEnum destination) __attribute__((always_inline)) {
        return this->destinations[destination];
    }

    // Usefull for MPE to retrieve PITCHBEND
    float getSource(SourceEnum source) {
        return this->sources[source];
    }

private:
    float sources[MATRIX_SOURCE_MAX];
    float destinations[DESTINATION_MAX];
    struct MatrixRowParams *rows;
};

#endif /* MATRIX_H_ */
