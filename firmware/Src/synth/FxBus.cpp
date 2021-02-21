/*
 * FxBus.cpp
 *
 *  Created on: Feb 14, 2021
 *      Author: patvig
 */

#include "FxBus.h"


//float DlyBuffer[FxBus::BufferSize] __attribute__((section(".ram_d1")));
//float DlyBuffer[FxBus::BufferSize] __attribute__((section(".instruction_ram")));
float FxBus::DlyBuffer[FxBus::BufferSize] __attribute__((section(".ram_d1")));

FxBus::FxBus() {
}

void FxBus::init(SynthState *synthState) {
    this->synthState_ = synthState;
	delayReadPos = delayWritePos + delaySampleCount;
    for (int s = 0; s < FxBus::BufferSize; s++) {
    	DlyBuffer[ s ] = 0;
    }
}

/**
 * init before sum timbres
 */
void FxBus::mixSumInit() {

	sample = getSampleBlock();
    for (int s = 0; s < BLOCK_SIZE; s++) {
    	*(sample++) = 0;
    	*(sample++) = 0;
    }

    float prevFxTime = fxTime;
    fxType =  synthState_->fullState.masterfxConfig[MASTERFX_TYPE];
    fxTime =  synthState_->fullState.masterfxConfig[MASTERFX_TIME];
    fxSpace =  synthState_->fullState.masterfxConfig[MASTERFX_SPACE];
    fxTone =  synthState_->fullState.masterfxConfig[MASTERFX_TONE];
    fxDiffusion =  synthState_->fullState.masterfxConfig[MASTERFX_DIFFUSION];
    fxWidth =  synthState_->fullState.masterfxConfig[MASTERFX_WIDTH];

    delaySampleCount = delaySampleCount * 0.95f + fxTime * FxBus::BufferSize * 0.05f;

    if(prevFxTime != fxType) {
    	delayReadPos = delayWritePos - delaySampleCount;
    	if( delayReadPos < 0 )
    		delayReadPos += FxBus::BufferSize;
    }

}

/**
 * add timbre to bus mix
 */
void FxBus::mixSum(float *inStereo, int timbreNum) {
	float level = synthState_->mixerState.instrumentState_[timbreNum].send;

	sample = getSampleBlock();

	for (int s = 0; s < BLOCK_SIZE; s++) {
    	*(sample++) += *inStereo++ * level;
    	*(sample++) += *inStereo++ * level;
    }
}

/**
 * process fx on bus mix
 */
void FxBus::processBlock(int32_t *outBuff) {
	sample = getSampleBlock();
    float sampleMultipler = (float) 0x7fffff;

    for (int s = 0; s < BLOCK_SIZE; s++) {
    	if( delayWritePos >= FxBus::BufferSize )
    		delayWritePos -= FxBus::BufferSize;

    	if( delayReadPos >= FxBus::BufferSize )
    		delayReadPos -= FxBus::BufferSize;

    	DlyBuffer[ delayWritePos++ ] = *(sample++) - DlyBuffer[ delayReadPos ] * fxSpace;
    	DlyBuffer[ delayWritePos++ ] = *(sample++) - DlyBuffer[ delayReadPos + 1 ] * fxSpace;

    	*(outBuff++) += (int32_t) ((DlyBuffer[ delayReadPos++ ] * sampleMultipler));
    	*(outBuff++) += (int32_t) ((DlyBuffer[ delayReadPos++ ] * sampleMultipler));
    }
}
