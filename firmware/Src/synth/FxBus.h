
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

	static const int BufferSize = 10000 * 2;

protected:

	float sampleBlock_[BLOCK_SIZE * 2];
	float *sample;
	static float DlyBuffer[BufferSize];
    int delayWritePos = 0;
    float delayReadPos = 0;
    int delayReadPosInt = 0;
    int delayReadPosInt2 = 0;
    float readSpeed = 1.98f;

    int fxType = 0;
    float fxTime = 0.98;
    float fxFeedback = 0.5;
    float fxTone = 0.5;
    float fxDiffusion = 0;
    float fxWidth =  0;

    int delaySampleCount = fxTime * BufferSize;


};

#endif	// end FX_BUS_
