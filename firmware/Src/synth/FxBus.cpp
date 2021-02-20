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
}

void FxBus::mixSumInit() {

	sample = getSampleBlock();
    for (int s = 0; s < BLOCK_SIZE; s++) {
    	*(sample++) = 0;
    	*(sample++) = 0;
    }
}

void FxBus::mixSum(float *inStereo, int timbreNum) {

	float level = synthState_->mixerState.instrumentState_[timbreNum].send;

	sample = getSampleBlock();

	for (int s = 0; s < BLOCK_SIZE; s++) {
    	*(sample++) += *inStereo++ * level;
    	*(sample++) += *inStereo++ * level;
    }
}

void FxBus::processBlock(int32_t *outBuff) {
	sample = getSampleBlock();
    float sampleMultipler = (float) 0x7fffff;

    for (int s = 0; s < BLOCK_SIZE; s++) {
    	if( delayWritePos >= FxBus::BufferSize )
    		delayWritePos = 0;

    	delayReadPos = delayWritePos + delaySampleCount;
    	if( delayReadPos >= FxBus::BufferSize )
    		delayReadPos -= FxBus::BufferSize;

    	DlyBuffer[ delayWritePos++ ] = *(sample++);
    	DlyBuffer[ delayWritePos++ ] = *(sample++);

    	*(outBuff++) += (int32_t) ((DlyBuffer[ delayReadPos++ ] * sampleMultipler));
    	*(outBuff++) += (int32_t) ((DlyBuffer[ delayReadPos++ ] * sampleMultipler));
    }
}
