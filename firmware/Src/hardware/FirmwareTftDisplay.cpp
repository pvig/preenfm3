/*
 * Copyright 2022 Xavier Hosxe
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



#include "FirmwareTftDisplay.h"
#include "Common.h"
#include "Env.h"
#include "Osc.h"
#include <math.h>

extern DMA2D_HandleTypeDef hdma2d;
extern RNG_HandleTypeDef hrng;

extern uint16_t tftPalette565[NUMBER_OF_TFT_COLORS];

const float randomDisplay[32] = { 0.679268211, -0.559346101, 0.93210829, 0.263967825,
        0.868206376, -0.955814486, 0.767822461, -0.260222359, -0.617473204,
        -0.844648436, 0.380078337, 0.481601148, -0.807739871, 0.995462356,
        -0.854993442, 0.704896648, 0.258727388, 0.947767668, -0.12564,
        0.513101208, -0.75714895, -0.11723284, 0.026869686, -0.738605819,
        -0.262181885, 0.856453, 0.209812963, -0.976879332, -0.019152814,
        -0.213303568, 0.55, 0.379348783 };

FirmwareTftDisplay::FirmwareTftDisplay() : TftDisplay() {
   oscilloIsClean = true;
   envInQueue = 0;
   lfoInQueue = 0;
   operatorInQueue = 0;
    operatorPhaseNormalized = 0.0f;
}

FirmwareTftDisplay::~FirmwareTftDisplay() {
    // TODO Auto-generated destructor stub
}


void FirmwareTftDisplay::clearActions() {
    tftActions.clear();
    envInQueue = 0;
    lfoInQueue = 0;
    operatorInQueue = 0;
}


void FirmwareTftDisplay::initWaveFormExt(int index, float* waveform, int size) {
    waveForm[index].waveForms = waveform;
    waveForm[index].size = size;
}


void FirmwareTftDisplay::additionalActions() {

    // Just in case
    if (tftActions.getCount() == 0) {
        envInQueue = 0;
        lfoInQueue = 0;
        operatorInQueue = 0;
    }

    switch (currentAction.actionType) {
    case TFT_DRAW_OSCILLO_BACKGROUND_WAVEFORM:
        switch (currentAction.param3) {
        case 0:
            // Copy bgOscillo2 to bgOscillo (that clears current waveform in oscilo)
            if (PFM_IS_DMA2D_READY()) {

                MODIFY_REG(hdma2d.Instance->CR, DMA2D_CR_MODE, DMA2D_R2M);
                WRITE_REG(hdma2d.Instance->OCOLR, tftPalette565[COLOR_BLACK]);

                MODIFY_REG(hdma2d.Instance->OOR, DMA2D_OOR_LO, 2);
                MODIFY_REG(hdma2d.Instance->NLR, (DMA2D_NLR_NL|DMA2D_NLR_PL),
                        (98 | (158 << DMA2D_POSITION_NLR_PL)));
                WRITE_REG(hdma2d.Instance->OMAR, (uint32_t )(bgOscillo + 161));

                PFM_START_DMA2D();
                currentAction.param3 = 1;
            }
            break;
        case 1:
            if (currentAction.param1 <= TFT_NUMBER_OF_WAVEFORM_EXT) {
                oscilloBgDrawOperatorShape(waveForm[currentAction.param1].waveForms, waveForm[currentAction.param1].size);
                operatorInQueue--;
            } else if (currentAction.param1 == TFT_DRAW_ENVELOPPE) {
                oscilloBgDrawEnvelope();
                envInQueue--;
            } else if (currentAction.param1 == TFT_DRAW_LFO) {
                oscilloBgDrawLfo();
                lfoInQueue--;
            }
            currentAction.actionType = 0;
            oscilloForceNextDisplay();
            break;
        }
    }
}



void FirmwareTftDisplay::oscilloBgDrawOperatorShape(float* waveForm, int size) {

    int indexMiddle = 50 * 160;
    uint16_t oscilloColor = tftPalette565[COLOR_YELLOW];

    for (int x = 0; x < 160; x++) {
        float index = x * ((float)size) / 160.0f;
        oscilloYValue[x] = (int) (waveForm[(int)index] * 48.0f);
    }

    for (int x = 1; x < 159; x++) {
        if (oscilloYValue[x] < (oscilloYValue[x + 1] - 1)) {
            for (int oy = oscilloYValue[x]; oy < oscilloYValue[x + 1]; oy++) {
                bgOscillo[indexMiddle + x + oy * 160] = oscilloColor;
            }
        } else if (oscilloYValue[x] > (oscilloYValue[x + 1] + 1)) {
            for (int oy = oscilloYValue[x]; oy > oscilloYValue[x + 1]; oy--) {
                bgOscillo[indexMiddle + x + oy * 160] = oscilloColor;
            }
        } else {
            bgOscillo[indexMiddle + x + oscilloYValue[x] * 160] = oscilloColor;
        }
    }

    // Subtle phase marker: dotted vertical line + small dot on waveform crossing.
    int markerX = (int)(operatorPhaseNormalized * 159.0f + 0.5f);
    if (markerX < 1) {
        markerX = 1;
    } else if (markerX > 158) {
        markerX = 158;
    }

    uint16_t markerColor = tftPalette565[COLOR_LIGHT_GRAY];
    for (int y = 2; y <= 95; y += 2) {
        bgOscillo[markerX + y * 160] = markerColor;
    }

    int markerY = oscilloYValue[markerX];
    if (markerY > -48 && markerY < 48) {
        uint16_t dotColor = tftPalette565[COLOR_CYAN];
        int dotIndex = indexMiddle + markerX + markerY * 160;
        bgOscillo[dotIndex] = dotColor;
        bgOscillo[dotIndex - 1] = dotColor;
        bgOscillo[dotIndex + 1] = dotColor;
    }
}



void FirmwareTftDisplay::oscilloBgActionClear() {
    if (!oscilloIsClean) {
        TFTAction newAction;
        newAction.param1 = 255;
        newAction.param3 = 0;
        newAction.actionType = TFT_DRAW_OSCILLO_BACKGROUND_WAVEFORM;
        tftActions.insert(newAction);
        oscilloIsClean = true;
    }
}


void FirmwareTftDisplay::oscilloBgActionOperatorShape(int wfNumber) {
    if (operatorInQueue >= 2) {
        // No need to ask for new drawing
        // params have been updated and next drawing will use them
        return;
    }
    operatorInQueue++;
    TFTAction newAction;
    newAction.actionType = TFT_DRAW_OSCILLO_BACKGROUND_WAVEFORM;
    newAction.param1 = wfNumber;
    newAction.param3 = 0;
    tftActions.insert(newAction);

    oscilloIsClean = false;
}

void FirmwareTftDisplay::oscilloBgActionLfo() {
    if (lfoInQueue >= 2) {
        // No need to ask for new drawing
        // params have been updated and next drawing will use them
        return;
    }
    lfoInQueue++;

    TFTAction newAction;
    newAction.actionType = TFT_DRAW_OSCILLO_BACKGROUND_WAVEFORM;
    newAction.param1 = TFT_DRAW_LFO;
    newAction.param3 = 0;
    tftActions.insert(newAction);

    oscilloIsClean = false;
}


void FirmwareTftDisplay::oscilloBgActionEnvelope() {
    if (envInQueue >= 2) {
        // No need to ask for new drawing
        // params have been updated and next drawing will use them
        return;
    }
    envInQueue++;
    TFTAction newAction;
    newAction.actionType = TFT_DRAW_OSCILLO_BACKGROUND_WAVEFORM;
    newAction.param1 = TFT_DRAW_ENVELOPPE;
    newAction.param3 = 0;
    tftActions.insert(newAction);

    oscilloIsClean = false;
}



void FirmwareTftDisplay::oscilloBgSetEnvelope(float a, float d, float s, float r, float aL, float dL, float sL, float rL, int8_t aCurve, int8_t dCurve, int8_t sCurve, int8_t rCurve) {
    oscilParams1[0] = a;
    oscilParams1[1] = d + oscilParams1[0];
    oscilParams1[2] = s + oscilParams1[1];
    oscilParams1[3] = r + oscilParams1[2];

    oscilParams2[0] = aL;
    oscilParams2[1] = dL;
    oscilParams2[2] = sL;
    oscilParams2[3] = rL;

    envCurve[0] = aCurve;
    envCurve[1] = dCurve;
    envCurve[2] = sCurve;
    envCurve[3] = rCurve;

    oscilParams1[4] = COLOR_YELLOW;

}

void FirmwareTftDisplay::oscilloBgSetLfoEnvelope(float a, float d, float s, float r, float aL, float dL, float sL, float rL) {
    oscilParams1[0] = a;
    oscilParams1[1] = d + oscilParams1[0];
    oscilParams1[2] = s + oscilParams1[1];
    oscilParams1[3] = r + oscilParams1[2];

    oscilParams2[0] = aL;
    oscilParams2[1] = dL;
    oscilParams2[2] = sL;
    oscilParams2[3] = rL;

    envCurve[0] = 1;
    envCurve[1] = 1;
    envCurve[2] = 1;
    envCurve[3] = 1;

    oscilParams1[4] = COLOR_BLUE;
}


void FirmwareTftDisplay::oscilloBgSetLfo(float shape, float freq, float kSyn, float syncMode, float bias, float phase) {
    if (freq >= 100.0f) {
        freq = 1.0f;
    }
    oscilParams1[0] = shape;
    oscilParams1[1] = freq;
    oscilParams1[2] = kSyn;
    oscilParams1[3] = syncMode;
    oscilParams1[4] = bias;
    oscilParams1[5] = phase;

}

void FirmwareTftDisplay::oscilloBgSetOperatorPhase(float phaseDegrees) {
    if (phaseDegrees < 0.0f) {
        phaseDegrees = 0.0f;
    } else if (phaseDegrees > 360.0f) {
        phaseDegrees = 360.0f;
    }
    operatorPhaseNormalized = phaseDegrees * (1.0f / 360.0f);
}


void FirmwareTftDisplay::oscilloBgDrawEnvStep(int16_t x1, int16_t y1, int16_t x2, int16_t y2, uint16_t color, int8_t curve) {
    if (curve == CURVE_TYPE_LIN || x1 == x2) {
        oscilloBgDrawLine(x1, y1, x2, y2, color);
    } else {
        int indexLow = 98 * 160 + 1;
        float width = x2 - x1;
        float invWidth = 1.0f / width;
        float height = y2 - y1;
        int oldPy = 0;
        int py;
        for (int x = x1; x <= x2; x++) {
            if (unlikely(x == x2)) {
                py = allLfoTables[curve].table[63] * height + y1;
            } else {
                int indexInt = x - x1;
                // Let's take only 62 value because we use index+1 bellow
                float index = 62.0f * (float)indexInt * invWidth;
                // Linear interpolation
                int index2 = (int)index;
                float indexRest = index - index2;
                float yValue1 = allLfoTables[curve].table[(int)index];
                float yValue2 = allLfoTables[curve].table[(int)index + 1];
                float yInterpolated = yValue1 * (1-indexRest) + yValue2 * indexRest;
                py = yInterpolated * height + y1;
            }
            // Draw vertical line between old y value and new one
            if (unlikely(x == x1)) {
                bgOscillo[indexLow + x - py * 160] = color;
            } else {
                if (oldPy <= py) {
                    for (int y = oldPy; y <= py; y += 1) {
                        bgOscillo[indexLow + x - y * 160] = color;
                    }
                } else {
                    for (int y = py; y <= oldPy; y += 1) {
                        bgOscillo[indexLow + x - y * 160] = color;
                    }
                }
            }
            oldPy = py;
        }
    }
}


void FirmwareTftDisplay::oscilloBgDrawLine(int16_t x1, int16_t y1, int16_t x2, int16_t y2, uint16_t color) {
    int16_t deltax = 0, deltay = 0, x = 0, y = 0, xinc1 = 0, xinc2 = 0, yinc1 = 0, yinc2 = 0, den = 0, num = 0, numadd = 0, numpixels = 0,
            curpixel = 0;
    int indexLow = 98 * 160 + 1;

    deltax = ABS(x2 - x1);
    deltay = ABS(y2 - y1);
    x = x1;
    y = y1;

    if (x2 >= x1) {
        xinc1 = 1;
        xinc2 = 1;
    } else {
        xinc1 = -1;
        xinc2 = -1;
    }

    if (y2 >= y1) {
        yinc1 = 1;
        yinc2 = 1;
    } else {
        yinc1 = -1;
        yinc2 = -1;
    }

    if (deltax >= deltay) {
        xinc1 = 0;
        yinc2 = 0;
        den = deltax;
        num = deltax / 2;
        numadd = deltay;
        numpixels = deltax;
    } else {
        xinc2 = 0;
        yinc1 = 0;
        den = deltay;
        num = deltay / 2;
        numadd = deltax;
        numpixels = deltay;
    }

    for (curpixel = 0; curpixel <= numpixels; curpixel++) {
        //  SETPIXEL(x, y);
        bgOscillo[indexLow + x - y * 160] = color;

        num += numadd;
        if (num >= den) {
            num -= den;
            x += xinc1;
            y += yinc1;
        }
        x += xinc2;
        y += yinc2;
    }
}


void FirmwareTftDisplay::oscilloBgDrawVerticalLine(int16_t x1, int16_t y1, int16_t y2, uint16_t color) {
    int indexLow = 98 * 160 + 1;
    for (int y = y1; y <= y2; y += 1) {
        bgOscillo[indexLow + x1 - y * 160] = color;
    }
}

void FirmwareTftDisplay::oscilloBgDrawEnvelope() {
    float div;

    if (oscilParams1[3] < 1.0f) {
        div = 1.0f;
    } else if (oscilParams1[3] < 2.0f) {
        div = 2.0f;
    } else {
        div = ((int)(oscilParams1[3] * .2f)) * 5.0 + 5.0f ;
    }

    float scale = 158.0f / div;
    uint16_t color = tftPalette565[(int)oscilParams1[4]];
    uint16_t darkGrey = tftPalette565[COLOR_DARK_GRAY];

    // Draw the enveloppe
    oscilloBgDrawVerticalLine((int)(oscilParams1[0] * scale), 0, (int)(oscilParams2[0] * 97.0f), darkGrey);
    oscilloBgDrawVerticalLine((int)(oscilParams1[1] * scale), 0, (int)(oscilParams2[1] * 97.0f), darkGrey);
    oscilloBgDrawVerticalLine((int)(oscilParams1[2] * scale), 0, (int)(oscilParams2[2] * 97.0f), darkGrey);


    oscilloBgDrawEnvStep(0, 0, (int)(oscilParams1[0] * scale), (int)(oscilParams2[0] * 97.0f), color, envCurve[0]);
    oscilloBgDrawEnvStep((int)(oscilParams1[0] * scale), (int)(oscilParams2[0] * 97.0f), (int)(oscilParams1[1] * scale), (int)(oscilParams2[1] * 97.0f), color, envCurve[1]);
    oscilloBgDrawEnvStep((int)(oscilParams1[1] * scale), (int)(oscilParams2[1] * 97.0f), (int)(oscilParams1[2] * scale), (int)(oscilParams2[2] * 97.0f), color, envCurve[2]);
    oscilloBgDrawEnvStep((int)(oscilParams1[2] * scale), (int)(oscilParams2[2] * 97.0f), (int)(oscilParams1[3] * scale), (int)(oscilParams2[3] * 97.0f), color, envCurve[3]);

    uint8_t digits[] = { (uint8_t)(((int)div) / 10),  (uint8_t)(((int)div) % 10) };
    for (int d = 0; d < 2; d++) {
        if ((d == 1) || (digits[d] != 0)) {
            int x = 136 + d * 11;
            int y = 2;
            const uint8_t* digitBits = tftAlgo->getDigitBits(digits[d]);

            for (int j = 0; j < 5; j++) {
                const uint8_t line = digitBits[j] >> 3;
                for (int i = 0; i < 5; i++) {
                    if (((line >> (4 - i)) & 1) == 1) {
                        bgOscillo[(x + i * 2)     + (y + j * 2) * 160] = color;
                        bgOscillo[(x + i * 2 + 1) + (y + j * 2) * 160] = color;
                        bgOscillo[(x + i * 2)     + (y + j * 2 + 1) * 160] = color;
                        bgOscillo[(x + i * 2 + 1) + (y + j * 2 + 1) * 160] = color;
                    }
                }
            }
        }
    }
}




/*
 * randType == 0 : random
 * randType == 1 : Brownian
 * randType == 2 : Wandering
 * randType == 3 : Flow
 */
void FirmwareTftDisplay::oscilloFillWithRand(int randtype) {
    // Rand
    float incIndex = 1.0f / 160.0f * oscilParams1[1];
    float sample = 0.0f;
    float nextSample = 0.0f;
    float sampleLp = 0.0f;
    int8_t sampleInt = 0;

    // Adjust to enter index>1 condition when x==0
    float index = oscilParams1[4] + 1.0f - incIndex;
    int randIndex = 0;
    for (int x = 0; x < 160; x++) {
        index += incIndex;

        if (unlikely(index >= 1.0f))  {
            index -= 1.0f;
            sample =  randomDisplay[(randIndex++) % 32] * 48.0f;
            switch (randtype) {
            case 0:
                sampleInt = (int) sample;
                break;
            case 1:
                sampleLp = sample * .4f + sampleLp * 0.6f;
                sampleInt = (int) sampleLp;
                break;
            case 2:
                sampleLp = sample;
                sample = nextSample;
                nextSample = sampleLp;
                break;
            case 3:
                sampleLp = sample * .4f + sampleLp * 0.6f;
                sample = nextSample;
                nextSample = sampleLp;
                break;
            }
        }
        if (randtype <= 1) {
            oscilloYValue[x] = sampleInt;
        } else {
            oscilloYValue[x] =  (int) ((nextSample - sample) * index + sample);
        }
    }
}


void FirmwareTftDisplay::oscilloBgDrawLfo() {
    // Sclae 160 pixel = 1 second of LFO
    int indexMiddle = 50 * 160;
    uint16_t oscilloColor = tftPalette565[COLOR_BLUE];
    float phaseForPreview = oscilParams1[5];
    int previewDelayPixels = 0;

    if (phaseForPreview < 0.0f) {
        float delaySeconds = -phaseForPreview;
        // Keep delay preview bounded to the 1-second oscilloscope window.
        previewDelayPixels = (int)(delaySeconds * 160.0f + 0.999f);
        if (previewDelayPixels > 160) {
            previewDelayPixels = 160;
        }
        phaseForPreview = 0.0f;
    }

    // Shape
    switch ((int)oscilParams1[0]) {
    case 0: {
        // Sin
        float *samples = waveTables[OSC_SHAPE_SIN].table;
        int size = waveTables[OSC_SHAPE_SIN].max + 1;
        for (int x = 0; x < 160; x++) {
            float index = ((float)x) * ((float)size) / 160.0f * oscilParams1[1] + size * phaseForPreview;
            int iIndex = index;
            iIndex %= size;
            if (iIndex < 0) {
                iIndex += size;
            }
            oscilloYValue[x] = (int) (samples[iIndex] * 47.0f);
        }
        break;
    }
    case 1: {
        // Saw
        int size = 200;
        float incSample = 0.005f;
        for (int x = 0; x < 160; x++) {
            float index = ((float)x) * ((float)size) / 160.0f * oscilParams1[1] + size * phaseForPreview;
            int iIndex = index;
            iIndex %= size;
            oscilloYValue[x] = -47 + (int) ((float)iIndex) * incSample * 94.0f;
        }
        break;
    }
    case 2: {
        // Triange
        int size = 200;
        // incSample twice 1/200 for triangle
        float incSample = 0.01f;
        for (int x = 0; x < 160; x++) {
            float index = ((float)x) * ((float)size) / 160.0f * oscilParams1[1] + size * phaseForPreview;
            int iIndex = index;
            iIndex %= size;
            if (iIndex < 100) {
                oscilloYValue[x] = (int) ((float)iIndex) * incSample * 47.0f;
            } else {
                oscilloYValue[x] = 95 - (int) ((float)iIndex) * incSample * 47.0f;
            }
            oscilloYValue[x] = - 47 + oscilloYValue[x] * 2;
        }
        break;
    }
    case 3: {
        // Square
        int size = 50;
        for (int x = 0; x < 160; x++) {
            float index = ((float)x) * ((float)size) / 160.0f * oscilParams1[1] + size * phaseForPreview;
            int iIndex = index;
            iIndex %= size;
            oscilloYValue[x] = iIndex < 25 ? 47 : -47;
        }
        break;
    }
    case 4:
    case 5:
    case 6:
    case 7:
        // Rand
        oscilloFillWithRand((int)oscilParams1[0] - 4);
        break;
    case 8: {
        // Falling saw
        int size = 200;
        float incSample = 0.005f;
        for (int x = 0; x < 160; x++) {
            float index = ((float)x) * ((float)size) / 160.0f * oscilParams1[1] + size * phaseForPreview;
            int iIndex = index;
            iIndex %= size;
            oscilloYValue[x] = 47 - (int) ((float)iIndex) * incSample * 94.0f;
        }
        break;
    }
    case 9: {
        // Exponential decay from max to min over one cycle.
        const float k = 6.0f;
        const float expEnd = expf(-k);
        for (int x = 0; x < 160; x++) {
            float phase = ((float)x) / 160.0f * oscilParams1[1] + phaseForPreview;
            phase -= (int)phase;
            if (phase < 0.0f) {
                phase += 1.0f;
            }
            float norm = (expf(-k * phase) - expEnd) / (1.0f - expEnd);
            oscilloYValue[x] = (int)((norm * 2.0f - 1.0f) * 47.0f);
        }
        break;
    }
    case 10: {
        // Log-like decay: starts gently then drops faster.
        const float a = 31.0f;
        const float invLog = 1.0f / logf(1.0f + a);
        for (int x = 0; x < 160; x++) {
            float phase = ((float)x) / 160.0f * oscilParams1[1] + phaseForPreview;
            phase -= (int)phase;
            if (phase < 0.0f) {
                phase += 1.0f;
            }
            float norm = logf(1.0f + a * (1.0f - phase)) * invLog;
            oscilloYValue[x] = (int)((norm * 2.0f - 1.0f) * 47.0f);
        }
        break;
    }
    case 11: {
        // Exponential rise from min to max.
        const float k = 6.0f;
        const float expEnd = expf(-k);
        for (int x = 0; x < 160; x++) {
            float phase = ((float)x) / 160.0f * oscilParams1[1] + phaseForPreview;
            phase -= (int)phase;
            if (phase < 0.0f) {
                phase += 1.0f;
            }
            float norm = (1.0f - expf(-k * phase)) / (1.0f - expEnd);
            oscilloYValue[x] = (int)((norm * 2.0f - 1.0f) * 47.0f);
        }
        break;
    }
    case 12: {
        // Log-like rise: starts gently then rises faster.
        const float a = 31.0f;
        const float invLog = 1.0f / logf(1.0f + a);
        for (int x = 0; x < 160; x++) {
            float phase = ((float)x) / 160.0f * oscilParams1[1] + phaseForPreview;
            phase -= (int)phase;
            if (phase < 0.0f) {
                phase += 1.0f;
            }
            float norm = (1.0f - logf(1.0f + a * (1.0f - phase)) * invLog);
            oscilloYValue[x] = (int)((norm * 2.0f - 1.0f) * 47.0f);
        }
        break;
    }
    case 13: {
        // Rounded attack-decay hump.
        for (int x = 0; x < 160; x++) {
            float phase = ((float)x) / 160.0f * oscilParams1[1] + phaseForPreview;
            phase -= (int)phase;
            if (phase < 0.0f) {
                phase += 1.0f;
            }
            float y;
            if (phase < 0.5f) {
                float t = phase * 2.0f;
                float s = t * t * (3.0f - 2.0f * t);
                y = -1.0f + 2.0f * s;
            } else {
                float t = (phase - 0.5f) * 2.0f;
                float s = t * t * (3.0f - 2.0f * t);
                y = 1.0f - 2.0f * s;
            }
            oscilloYValue[x] = (int)(y * 47.0f);
        }
        break;
    }
    case 14: {
        // Rounded attack-hold-decay.
        for (int x = 0; x < 160; x++) {
            float phase = ((float)x) / 160.0f * oscilParams1[1] + phaseForPreview;
            phase -= (int)phase;
            if (phase < 0.0f) {
                phase += 1.0f;
            }
            float y;
            if (phase < 0.25f) {
                float t = phase * 4.0f;
                float s = t * t * (3.0f - 2.0f * t);
                y = -1.0f + 2.0f * s;
            } else if (phase < 0.5f) {
                y = 1.0f;
            } else {
                float t = (phase - 0.5f) * 2.0f;
                float s = t * t * (3.0f - 2.0f * t);
                y = 1.0f - 2.0f * s;
            }
            oscilloYValue[x] = (int)(y * 47.0f);
        }
        break;
    }
    case 15: {
        // S-curve decay.
        for (int x = 0; x < 160; x++) {
            float phase = ((float)x) / 160.0f * oscilParams1[1] + phaseForPreview;
            phase -= (int)phase;
            if (phase < 0.0f) {
                phase += 1.0f;
            }
            float s = phase * phase * (3.0f - 2.0f * phase);
            float y = 1.0f - 2.0f * s;
            oscilloYValue[x] = (int)(y * 47.0f);
        }
        break;
    }
    case 16: {
        // Buchla-like "plong": fast rounded attack, short hold, curved decay.
        const float attack = 0.05f;
        const float hold = 0.06f;
        const float k = 7.0f;
        const float expEnd = expf(-k);
        for (int x = 0; x < 160; x++) {
            float phase = ((float)x) / 160.0f * oscilParams1[1] + phaseForPreview;
            phase -= (int)phase;
            if (phase < 0.0f) {
                phase += 1.0f;
            }
            float y;
            if (phase < attack) {
                float t = phase / attack;
                float s = t * t * (3.0f - 2.0f * t);
                y = -1.0f + 2.0f * s;
            } else if (phase < attack + hold) {
                y = 1.0f;
            } else {
                float d = (phase - attack - hold) / (1.0f - attack - hold);
                float norm = (expf(-k * d) - expEnd) / (1.0f - expEnd);
                y = norm * 2.0f - 1.0f;
                y += 0.10f * expf(-18.0f * d) * sinf(18.849556f * d);
                if (y > 1.0f) {
                    y = 1.0f;
                } else if (y < -1.0f) {
                    y = -1.0f;
                }
            }
            oscilloYValue[x] = (int)(y * 47.0f);
        }
        break;
    }
    case 17: {
        // Buchla-like "plong" variant with doubled peak hold.
        const float attack = 0.05f;
        const float hold = 0.12f;
        const float k = 7.0f;
        const float expEnd = expf(-k);
        for (int x = 0; x < 160; x++) {
            float phase = ((float)x) / 160.0f * oscilParams1[1] + phaseForPreview;
            phase -= (int)phase;
            if (phase < 0.0f) {
                phase += 1.0f;
            }
            float y;
            if (phase < attack) {
                float t = phase / attack;
                float s = t * t * (3.0f - 2.0f * t);
                y = -1.0f + 2.0f * s;
            } else if (phase < attack + hold) {
                y = 1.0f;
            } else {
                float d = (phase - attack - hold) / (1.0f - attack - hold);
                float norm = (expf(-k * d) - expEnd) / (1.0f - expEnd);
                y = norm * 2.0f - 1.0f;
                y += 0.10f * expf(-18.0f * d) * sinf(18.849556f * d);
                if (y > 1.0f) {
                    y = 1.0f;
                } else if (y < -1.0f) {
                    y = -1.0f;
                }
            }
            oscilloYValue[x] = (int)(y * 47.0f);
        }
        break;
    }
    case 18:
    case 19:
    case 20:
    case 21:
    case 22:
    case 23:
    case 24:
    case 25:
    case 26: {
        // Extra LFO wavetable shapes from oscillator tables.
        int oscShape;
        if ((int)oscilParams1[0] <= 20) {
            oscShape = OSC_SHAPE_SIN_SQUARE + ((int)oscilParams1[0] - 18);
        } else {
            oscShape = OSC_SHAPE_USER1 + ((int)oscilParams1[0] - 21);
        }
        float* samples = waveTables[oscShape].table;
        int size = waveTables[oscShape].max + 1;
        for (int x = 0; x < 160; x++) {
            float index = ((float)x) * ((float)size) / 160.0f * oscilParams1[1] + size * phaseForPreview;
            int iIndex = (int)index;
            iIndex %= size;
            if (iIndex < 0) {
                iIndex += size;
            }
            oscilloYValue[x] = (int)(samples[iIndex] * 47.0f);
        }
        break;
    }
    }

    int syncMode = (int)(oscilParams1[3] + 0.5f);
    int shotCount = 0;
    if (syncMode >= LFO_SYNC_ONESHOT_INTERNAL_1 && syncMode <= LFO_SYNC_ONESHOT_INTERNAL_8) {
        shotCount = syncMode - LFO_SYNC_ONESHOT_INTERNAL_1 + 1;
    } else if (syncMode >= LFO_SYNC_ONESHOT_EXTERNAL_1 && syncMode <= LFO_SYNC_ONESHOT_EXTERNAL_8) {
        shotCount = syncMode - LFO_SYNC_ONESHOT_EXTERNAL_1 + 1;
    }

    float ksyncAmount = oscilParams1[2];

    // One-shot preview with hold section.
    if (shotCount > 0 && oscilParams1[1] > 0.0f) {

        float phase0 = phaseForPreview;
        while (phase0 >= 1.0f) {
            phase0 -= 1.0f;
        }
        while (phase0 < 0.0f) {
            phase0 += 1.0f;
        }

        float cyclesToHold = (float)shotCount - phase0;
        if (cyclesToHold < 0.0f) {
            cyclesToHold = 0.0f;
        }

        // Convert hold time to pixel index with ceil-like behavior so the
        // one-shot waveform is never truncated early.
        float holdPixels = cyclesToHold * 160.0f / oscilParams1[1];
        int holdX = (int)holdPixels;
        if ((float)holdX < holdPixels) {
            holdX++;
        }
        if (holdX < 0) {
            holdX = 0;
        }

        if (holdX < 160) {
            int holdY;
            switch ((int)oscilParams1[0]) {
            case 0: // Sin
                holdY = 0;
                break;
            case 1: // Saw
                holdY = 47;
                break;
            case 8: // Falling saw
                holdY = -47;
                break;
            case 2: // Triangle
                holdY = -47;
                break;
            case 9: // Exp decay
            case 10: // Log decay
            case 13: // AD
            case 14: // AHD
            case 15: // S decay
            case 16: // Buchla plong
            case 17: // Buchla plong2
                holdY = -47;
                break;
            case 11: // Exp rise
            case 12: // Log rise
                holdY = 47;
                break;
            case 3: // Square
                holdY = 47;
                break;
            case 18: // SinSquare
                holdY = (int)(waveTables[OSC_SHAPE_SIN_SQUARE].table[0] * 47.0f);
                break;
            case 19: // SinZero
                holdY = (int)(waveTables[OSC_SHAPE_SIN_ZERO].table[0] * 47.0f);
                break;
            case 20: // SinPos
                holdY = (int)(waveTables[OSC_SHAPE_SIN_POS].table[0] * 47.0f);
                break;
            case 21: // User1
                holdY = (int)(waveTables[OSC_SHAPE_USER1].table[0] * 47.0f);
                break;
            case 22: // User2
                holdY = (int)(waveTables[OSC_SHAPE_USER2].table[0] * 47.0f);
                break;
            case 23: // User3
                holdY = (int)(waveTables[OSC_SHAPE_USER3].table[0] * 47.0f);
                break;
            case 24: // User4
                holdY = (int)(waveTables[OSC_SHAPE_USER4].table[0] * 47.0f);
                break;
            case 25: // User5
                holdY = (int)(waveTables[OSC_SHAPE_USER5].table[0] * 47.0f);
                break;
            case 26: // User6
                holdY = (int)(waveTables[OSC_SHAPE_USER6].table[0] * 47.0f);
                break;
            default:
                // Random family: keep the last reached value.
                holdY = holdX > 0 ? oscilloYValue[holdX - 1] : oscilloYValue[0];
                break;
            }

            for (int x = holdX; x < 160; x++) {
                oscilloYValue[x] = holdY;
            }
        }
    }

    if (previewDelayPixels > 0) {
        int holdY = oscilloYValue[0];

        for (int x = 159; x >= previewDelayPixels; x--) {
            oscilloYValue[x] = oscilloYValue[x - previewDelayPixels];
        }
        for (int x = 0; x < previewDelayPixels; x++) {
            oscilloYValue[x] = holdY;
        }
    }

    if (ksyncAmount > 0.0f) {
        // Ksyn
        // 1/160 = 0.00627
        float kSyncInc = 1.0f / (160.0f * ksyncAmount);
        float kSync = 0;
        for (int x = previewDelayPixels; x < 160; x++) {
            if (kSync < 1) {
                oscilloYValue[x] = ((float)oscilloYValue[x] * kSync);
            } else {
                break;
            }
            kSync += kSyncInc;
        }
    }

    if (oscilParams1[4] != 0.0f) {
        for (int x = 0; x < 160; x++) {
            oscilloYValue[x] += (oscilParams1[4] * 48.0f);
        }
    }


    for (int x = 1; x < 159; x++) {
        if (oscilloYValue[x] < (oscilloYValue[x + 1] - 1)) {
            for (int oy = oscilloYValue[x]; oy < oscilloYValue[x + 1]; oy++) {
                if (likely(oy > -48 && oy < 48)) {
                    bgOscillo[indexMiddle + x - oy * 160] = oscilloColor;
                }
            }
        } else if (oscilloYValue[x] > (oscilloYValue[x + 1] + 1)) {
            for (int oy = oscilloYValue[x]; oy > oscilloYValue[x + 1]; oy--) {
                if (likely(oy > -48 && oy < 48)) {
                    bgOscillo[indexMiddle + x - oy * 160] = oscilloColor;
                }
            }
        } else {
            if (likely(oscilloYValue[x] > -48 && oscilloYValue[x] < 48)) {
                bgOscillo[indexMiddle + x - oscilloYValue[x] * 160] = oscilloColor;
            }
        }
    }

}
