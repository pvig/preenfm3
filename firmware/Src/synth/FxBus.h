
#ifndef FX_BUS_
#define FX_BUS_

#include "SynthStateAware.h"

class FxBus : public SynthStateAware {
public:
	FxBus();
	virtual ~FxBus() {}
    void init(SynthState *synthState);

	void mixSumInit();
    void mixSum(float *inStereo, int timbreNum);
	void processBlock(int32_t *outBuff);

	float* getSampleBlock() {
	    return sampleBlock_;
	}

	const float* getSampleBlock() const {
	    return sampleBlock_;
	}

	static const int BufferSize = 2048 * 2;

protected:

	float sampleBlock_[BLOCK_SIZE * 2];
	float *sample;
	static float DlyBuffer[BufferSize];
    int delayWritePos = 0;
    int delayReadPos;
    int delaySampleCount = 2000;
};

#endif	// end FX_BUS_
