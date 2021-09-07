#include "hpr.h"
#include "IPlug_include_in_plug_src.h"
#include "Utils.h"
#include "styles.h"

#if IPLUG_DSP
hpr::hpr(const InstanceInfo& info)
	: Plugin(info, MakeConfig(kNumParams, kNumPresets))
	, complexBuffer(Complex(), FFT_SIZE)
{
	//GetParam(kGain)->InitDouble("Gain", 50.0, 0.0, 100.0, 0.01, "%");
	GetParam(kH)->InitDouble("H", 50.0, 0., 100.0, 0.01, "%");
	GetParam(kP)->InitDouble("P", 50.0, 0., 100.0, 0.01, "%");
	GetParam(kR)->InitDouble("R", 10.0, 0., 100.0, 0.01, "%");

	sampleBuffer.resize(FFT_SIZE, 0.0);
	workerBuffer.resize(FFT_SIZE, 0.0);
	workerBuffer2.resize(FFT_SIZE, 0.0);
	nextOutput.resize(FFT_SIZE, 0.0);
	harmonics.resize(FFT_SIZE, 0.0);
	percussive.resize(FFT_SIZE, 0.0);
	residual.resize(FFT_SIZE, 0.0);

	start.resize(FFT_SIZE / 4, 0.0);

#if IPLUG_EDITOR // http://bit.ly/2S64BDd
	mMakeGraphicsFunc = [&]() {
		return MakeGraphics(*this, PLUG_WIDTH, PLUG_HEIGHT, PLUG_FPS, GetScaleForScreen(PLUG_WIDTH, PLUG_HEIGHT));
	};

	mLayoutFunc = [&](IGraphics* pGraphics) {
		//pGraphics->AttachCornerResizer(EUIResizerMode::Scale, false);
		pGraphics->AttachPanelBackground(IColor::FromColorCode(Colors::BG));
		pGraphics->LoadFont("Roboto-Regular", ROBOTO_FN);

		const IRECT left = pGraphics->GetBounds().GetFromLeft(350).GetCentredInside(pGraphics->GetBounds());
		const int nRows = 1;
		const int nCols = 4;

		int cellIdx = -1;

		auto nextCell = [&]() {
			return left.GetGridCell(++cellIdx, nRows, nCols).GetPadded(-5.);
		};
		pGraphics->AttachControl(new IVKnobControl(nextCell(), kH, "", knob1Style()));
		pGraphics->AttachControl(new IVKnobControl(nextCell(), kP, "", knob2Style()));
		pGraphics->AttachControl(new IVKnobControl(nextCell(), kR, "", knob3Style()));
	};
#if DEBUG
	runTests();
#endif
#endif
}

#define TEST_FFT 0
void hpr::ProcessBlock(sample** inputs, sample** outputs, int nFrames)
{
	//const double gain = GetParam(kGain)->Value() / 100.;
	H = GetParam(kH)->Value() / 100.0;
	P = GetParam(kP)->Value() / 100.0;
	R = GetParam(kR)->Value() / 100.0;

	const size_t nChans = NOutChansConnected();

	constexpr size_t windowL = FFT_SIZE / 2;
	constexpr size_t hopSize = windowL / 2;

	// nFrames should be checked and aggrued with a circBuffer or similar to equate the spectrum seg length
	// now we're taking a shortcut and expecting a appropriately sized frame size
	// This is possible to do by tweaking the sound card and/or buffer settings in 
	// e.g. Reaper
	assert(FFT_SIZE == nFrames);

	{ // Tests
		// test audio and latency
		//for (auto i = 0; i < nChans; ++i)
		//	for (auto j = 0; j < nFrames; ++j)
		//		outputs[i][j] = inputs[i][j];
	
		// test mono-stereo functions
		//fromMono(outputs, sampleBuffer.data(), nChans, nFrames);
		//toMono(sampleBuffer.data(), inputs, nChans, nFrames);
		//return;
	}

	std::fill_n(sampleBuffer.begin(), nFrames, 0.0);
	std::fill_n(workerBuffer.begin(), nFrames, 0.0);

	// STFT and OLA
	// First half of a window is handled to match the block
	// overlapping frame (second half)
	toMono(sampleBuffer.data(), inputs, nChans, hopSize);

	// concatenate the samples for the block-overlapping frame
	for (auto i = 0; i < hopSize; ++i)
		workerBuffer[i] = start[i];
	for (auto i = 0; i < hopSize; ++i)
		workerBuffer[i + hopSize] = sampleBuffer[i];

	hann(sampleBuffer.data(), workerBuffer.data(), windowL);

	// Frequency-domain
	fftBuffer(sampleBuffer, windowL);

#if TEST_FFT
	ifftBuffer(sampleBuffer, windowL);
#else
	// Harmonics
	harmonicsMatrix.nextMedians(harmonics.data(), workerBuffer.data(), sampleBuffer.data(), windowL);

	// Percussive
	std::fill(workerBuffer.begin(), workerBuffer.end(), 0.0);
	medianPercussive(percussive.data(), workerBuffer.data(), sampleBuffer.data(), windowL);

	// Combine with parameters and calc residual
	calculateHPR(windowL);
#endif


	// write beginning of this window to the end of the buffer to be output now
	for (size_t i = 0; i < hopSize; ++i)
    	nextOutput[i + wi] += sampleBuffer[i];

	// Output
	if (initted)
		// With nFrames latency, we write the concluded buffer here
		fromMono(outputs, nextOutput.data(), nChans, nFrames);
	else
		initted = true;

	// write end of this window to the beginning of the buffer to be output next
	std::fill_n(nextOutput.begin(), nFrames, 0.0);
	for (size_t i = 0; i < hopSize; ++i)
		nextOutput[i] = sampleBuffer[hopSize + i];

	// Handling of the windows that fit completely inside this frame
	const int endIndexOfHops = nFrames - hopSize;
	for (wi = 0; wi < endIndexOfHops; wi += hopSize)
	{
		toMono(sampleBuffer.data(), inputs, nChans, windowL, wi);
		hann(sampleBuffer.data(), windowL);

		// Frequency-domain
		fftBuffer(sampleBuffer, windowL);

#if TEST_FFT
		ifftBuffer(sampleBuffer, windowL);
#else
		// Harmonics
		harmonicsMatrix.nextMedians(harmonics.data(), workerBuffer.data(), sampleBuffer.data(), windowL);
	
		// Percussive
		std::fill(workerBuffer.begin(), workerBuffer.end(), 0.0);
		medianPercussive(percussive.data(), workerBuffer.data(), sampleBuffer.data(), windowL);

		// Combine with parameters
		calculateHPR(windowL);
#endif

		for (auto i = 0; i < windowL; ++i)
			nextOutput[i + wi] += sampleBuffer[i];
	}

	// Handling of the block-overlapping frame (first half)
	// save last chunk for next iteration (OLA)
	toMono(start.data(), inputs, nChans, hopSize, wi);

	return;
}
#endif

