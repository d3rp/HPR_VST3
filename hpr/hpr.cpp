#include "hpr.h"
#include "IPlug_include_in_plug_src.h"

enum Colors
{
	BG = 0x343838,
	FG1 = 0x005F6B,
	FG2 = 0x008C9E,
	FG3 = 0x00B4CC,
	FG4 = 0x00DFFC
};

IColor getColor(Colors c) { return IColor::FromColorCode(c); }

IVStyle knobStyle()
{
	// Background
	// OFF/Foreground
	// ON/Pressed
	// Frame
	// Highlight
	// Shadow
	// Extra 1
	// Extra 2
	// Extra 3
	IVStyle res = DEFAULT_STYLE
		.WithColor(kBG, getColor(Colors::FG4))
		.WithColor(kFG, getColor(Colors::FG3))
		.WithColor(kPR, getColor(Colors::FG3))
		.WithColor(kHL, COLOR_RED)
		.WithColor(kSH, getColor(Colors::FG1))
		.WithColor(kFR, getColor(Colors::FG3))
		.WithColor(kX1, COLOR_WHITE)
		.WithColor(kX2, COLOR_GREEN)
		.WithColor(kX3, COLOR_BLUE)
		;
	res.labelText = IText(20, getColor(Colors::FG2));
	res.valueText = IText(16, getColor(Colors::FG1));

	return res;
}
sample median(sample * const srcDst, const int n)
{
	auto m = n / 2;
	std::nth_element(srcDst, srcDst + m, srcDst + n);
	return srcDst[m];
}


sample median(sample* const worker, sample const* const src, const int n)
{
	std::copy_n(src, n, worker);
	return median(worker, n);
}

void testMedian()
{
	{
		std::vector<sample> v{ 5, 6, 4, 3, 2, 6, 7, 9, 3 };
		std::vector<sample> w(v.size());
		auto m = median(w.data(), v.data(), v.size());
		assert(m == 5);
	}
	{
		std::vector<sample> v{ 15, 6, 4, 3, 2, 6, 7, 9, 3, -1 };
		std::vector<sample> w(v.size());
		auto m = median(w.data(), v.data(), v.size());
		assert(m == 6);
	}
}

hpr::hpr(const InstanceInfo& info)
	: Plugin(info, MakeConfig(kNumParams, kNumPresets))
	, complexBuffer(Complex(), FFT_SIZE)
	, olaBuffer(FFT_SIZE * 0.5)
{

	//GetParam(kGain)->InitDouble("Gain", 50.0, 0.0, 100.0, 0.01, "%");
	GetParam(kH)->InitDouble("H", 50.0, 0., 100.0, 0.01, "%");
	GetParam(kP)->InitDouble("P", 50.0, 0., 100.0, 0.01, "%");
	GetParam(kR)->InitDouble("R", 10.0, 0., 100.0, 0.01, "%");

	sampleBuffer.resize(FFT_SIZE, 0.0);
	workerBuffer.resize(FFT_SIZE, 0.0);
	previousBuffer.resize(FFT_SIZE, 0.0);
	nextOutputs.resize(FFT_SIZE, 0.0);
	harmonics.resize(FFT_SIZE, 0.0);
	percussive.resize(FFT_SIZE, 0.0);

	start.resize(FFT_SIZE / 4, 0.0);
	end.resize(FFT_SIZE / 4, 0.0);

#if IPLUG_EDITOR // http://bit.ly/2S64BDd
	mMakeGraphicsFunc = [&]() {
		return MakeGraphics(*this, PLUG_WIDTH, PLUG_HEIGHT, PLUG_FPS, GetScaleForScreen(PLUG_WIDTH, PLUG_HEIGHT));
	};

	mLayoutFunc = [&](IGraphics* pGraphics) {
		//pGraphics->AttachCornerResizer(EUIResizerMode::Scale, false);
		pGraphics->AttachPanelBackground(IColor::FromColorCode(Colors::FG4));
		pGraphics->LoadFont("Roboto-Regular", ROBOTO_FN);

		const IRECT left = pGraphics->GetBounds().GetFromLeft(350);
		const int nRows = 1;
		const int nCols = 4;

		int cellIdx = -1;

		auto nextCell = [&]() {
			return left.GetGridCell(++cellIdx, nRows, nCols).GetPadded(-5.);
		};
		pGraphics->AttachControl(new IVKnobControl(nextCell(), kH, "", knobStyle()));
		pGraphics->AttachControl(new IVKnobControl(nextCell(), kP, "", knobStyle()));
		pGraphics->AttachControl(new IVKnobControl(nextCell(), kR, "", knobStyle()));
	};
	//runTests();
#endif
}

void runTests()
{
	testMedian();

}

void medianPercussive(sample* const dst, sample* const worker, sample const * const src, int n)
{
	const int L = PERCUSSIVE_MEDIAN_LEN;
	const int halfP = L / 2;
	const int N = FFT_SIZE - halfP - 1;
	for (auto i = halfP; i < N; ++i)
	{
		const auto* s = &src[i - halfP];
		dst[i] = median(worker, s, L);
	}
}

void calculateResidual(sample* const harmonics, sample* const percussive, sample* const residual, sample const* const src, int nFrames)
{
	for (auto i = 0; i < nFrames; ++i)
	{
		residual[i] = src[i] - harmonics[i];
		residual[i] -= percussive[i];
	}
}

#if IPLUG_DSP
void hpr::OnIdle()
{
	//spectrumSender.TransmitData(*this);
	//textSender.TransmitData(*this);
	//textSender2.TransmitData(*this);
}

void hann(sample* const dst, sample const* const src, const int n)
{
	for (int i = 0; i < n; i++) {
		double multiplier = 0.5 * (1 - cos(2 * PI * i / (n - 1)));
		dst[i] = multiplier * src[i];
	}
}

/** in-place Hanning */
void hann(sample* const dst, const int n)
{
	for (int i = 0; i < n; i++) {
		double multiplier = 0.5 * (1 - cos(2 * PI * i / (n - 1)));
		dst[i] = multiplier * dst[i];
	}
}


// void normaliseForLogScale(sample* const srcDst, sample* const dampingBuffer, int nFrames)
// {
// 	constexpr int halfSpectrum = SPECTRUM_REAL_SIZE;
// 	constexpr double ratio = 22000.0 / halfSpectrum;
// 	static const double highF = std::log10(22000);
// 	static const double lowF = std::log10(100);
// 	static const double fRange = highF - lowF;
// 	for (int s = 0; s < halfSpectrum; ++s)
// 	{
// 		int idx = std::log2(1.0 + s / (double)halfSpectrum) * halfSpectrum;
// 		//int idx = std::log10(1.0 + ((s * ratio) / 261.0)) * 261; // wikipedia frequency log scale
// 		//int idx = 511 * (std::log10(1.0 + s * ratio) - lowF) / fRange; // stackoverflow copypaste
// 		sample newS = (srcDst[s] / 50.0) - 1.0;
// 		sample oldS = dampingBuffer[idx];
// 		srcDst[idx] = std::max(newS, 0.9 * (oldS + 1.0) - 1.0);
// 	}
// }
// void spectrogram(CircularBuffer* dst, std::valarray<Complex>& const src, const int n) noexcept
// {
// 	// takes fft'd buffer, calc energy in other buffer
// 	for (auto i = 0; i < n; ++i)
// 	{
// 		const auto x = src[i];
// 		const auto real = x.real();
// 		const auto im = x.imag();
// 		dst->push(10 * std::log(real * real + im * im));
// 	}
// }
void fromReal(std::valarray<Complex>& dst, sample const* const src, const int n)
{
	for (int s = 0; s < n; s++)
		dst[s] = Complex(src[s], 0);
}
void toReal(sample* dst, std::valarray<Complex> const& src, const int n)
{
	for (int s = 0; s < n; s++)
		dst[s] = src[s].real();
}

void fromMono(sample** dst, sample* src, const int nChans, const int n)
{
	for (int ch = 0; ch < nChans; ++ch)
		for (int s = 0; s < n; ++s)
			dst[ch][s] = 0.0;
	
	for (int ch = 0; ch < nChans; ++ch)
		for (int s = 0; s < n; ++s)
			dst[ch][s] += src[s];

}
void toMono(sample* dst, sample** src, const int nChans, const int n, const int offset = 0)
{
	const sample r = 1.0 / (sample)nChans;

	// First channel
	for (auto i = 0; i < n; ++i)
		dst[i] = src[0][offset + i] * r;

	// Consecutive channels
	int i = 0;
	for (int c = 1; c < nChans; ++c)
	{
		for (auto i = 0; i < n; ++i) 
			dst[i] = dst[i] + src[c][offset + i] * r;
	}
}

void addWithParams(sample* const dst, sample* harm, sample* perc, const double H, const double P, const double R, int n)
{
	const double scale = std::max(1.0, H + P == 0 ? 0 : 1.0 / (H + P));

	for (auto i = 0; i < n; ++i)
	 	dst[i] = harm[i] * H;

	for (auto i = 0; i < n; ++i)
		dst[i] += perc[i] * P;

	for (auto i = 0; i < n; ++i)
		dst[i] *= scale;
}

void hpr::ProcessBlock(sample** inputs, sample** outputs, int nFrames)
{
	//const double gain = GetParam(kGain)->Value() / 100.;
	const double HGain = GetParam(kH)->Value() / 100.0;
	const double PGain = GetParam(kP)->Value() / 100.0;
	const double RGain = GetParam(kR)->Value() / 100.0;
	const size_t nChans = NOutChansConnected();
	constexpr size_t windowL = FFT_SIZE / 2;
	constexpr size_t hopSize = windowL / 2;

	// nFrames should be checked and aggrued with a circBuffer or similar to equate the spectrum seg length
	// now we're taking a shortcut and expecting a appropriately sized frame size
	// This is possible to do by tweaking the sound card and/or buffer settings in 
	// e.g. Reaper
	assert(FFT_SIZE == nFrames);

	// test audio and latency
	//for (auto i = 0; i < nChans; ++i)
	//	for (auto j = 0; j < nFrames; ++j)
	//		outputs[i][j] = inputs[i][j];

	// test mono-stereo functions
	//fromMono(outputs, sampleBuffer.data(), nChans, nFrames);
	//toMono(sampleBuffer.data(), inputs, nChans, nFrames);
	//return;

	std::fill_n(sampleBuffer.begin(), nFrames, 0.0);
	std::fill_n(workerBuffer.begin(), nFrames, 0.0);
	// OLA
	toMono(sampleBuffer.data(), inputs, nChans, hopSize);

	// Last half
	// olaBuffer.push(sampleBuffer.data(), hopSize);
	// olaBuffer.getShifted(workerBuffer.data(), windowL);

	for (auto i = 0; i < hopSize; ++i)
		workerBuffer[i] = start[i];
	//std::copy_n(start.begin(), hopSize, workerBuffer.begin());
	for (auto i = 0; i < hopSize; ++i)
		workerBuffer[i + hopSize] = sampleBuffer[i];
	//std::copy_n(sampleBuffer.begin(), hopSize, workerBuffer.begin() + hopSize);

	hann(workerBuffer.data(), windowL);
	fromReal(complexBuffer, workerBuffer.data(), windowL);
	fft(complexBuffer);

	toReal(sampleBuffer.data(), complexBuffer, windowL);
	{ // Harmonics
		harmonicsMatrix.nextMedians(harmonics.data(), workerBuffer.data(), sampleBuffer.data(), windowL);
		//sampleBuffer = harmonics;
	}
	//fromReal(complexBuffer, sampleBuffer.data(), windowL);

	{ // Percussive
		std::fill(workerBuffer.begin(), workerBuffer.end(), 0.0);
		medianPercussive(percussive.data(), workerBuffer.data(), sampleBuffer.data(), windowL);
	}

	//std::copy_n(sampleBuffer.begin(), windowL, harmonics.begin());
	addWithParams(sampleBuffer.data(), harmonics.data(), percussive.data(), HGain, PGain, RGain, windowL);
	fromReal(complexBuffer, sampleBuffer.data(), windowL);

	ifft(complexBuffer);
	toReal(workerBuffer.data(), complexBuffer, windowL);

	// write beginning of this window to the end of the buffer to be output now
	for (size_t i = 0; i < hopSize; ++i)
    	previousBuffer[i + wi] += workerBuffer[i];

	// With nFrames latency, we write the concluded buffer here
	if (initted)
		fromMono(outputs, previousBuffer.data(), nChans, nFrames);
	else
		initted = true;

	// write end of this window to the beginning of the buffer to be output next
	std::fill_n(previousBuffer.begin(), nFrames, 0.0);
	for (size_t i = 0; i < hopSize; ++i)
		previousBuffer[i] = workerBuffer[hopSize + i];

	const int endIndexOfHops = 1 + nFrames - windowL;
	for (wi = 0; wi < endIndexOfHops; wi += hopSize)
	{
		toMono(sampleBuffer.data(), inputs, nChans, windowL, wi);
		hann(sampleBuffer.data(), windowL);
		fromReal(complexBuffer, sampleBuffer.data(), windowL);
		fft(complexBuffer);

		toReal(sampleBuffer.data(), complexBuffer, windowL);
		{ // Harmonics
			harmonicsMatrix.nextMedians(harmonics.data(), workerBuffer.data(), sampleBuffer.data(), windowL);
			//sampleBuffer = harmonics;
		}
		for (auto& h : harmonics)
		{
			assert(h < 4.0);
		}

		{ // Percussive
			std::fill(workerBuffer.begin(), workerBuffer.end(), 0.0);
			medianPercussive(percussive.data(), workerBuffer.data(), sampleBuffer.data(), windowL);
		}

		//std::copy_n(sampleBuffer.begin(), windowL, harmonics.begin());
		addWithParams(sampleBuffer.data(), harmonics.data(), percussive.data(), HGain, PGain, RGain, windowL);
		fromReal(complexBuffer, sampleBuffer.data(), windowL);

		ifft(complexBuffer);
		toReal(sampleBuffer.data(), complexBuffer, windowL);
		for (auto i = 0; i < windowL; ++i)
			previousBuffer[i + wi] += sampleBuffer[i];
	}

	// Push_back first half for the next round
	//wi += hopSize;
	// for (auto i = 0; i < hopSize; ++i)
	// 	olaBuffer.push(sampleBuffer.data() + hopSize, hopSize);

	toMono(start.data(), inputs, nChans, hopSize, wi);
	//std::copy_n(sampleBuffer.begin() + hopSize, hopSize, start.begin());

	return;
	//////////////////////////////////////////////////
	//  Harmonics
	//for (auto i = 0; i < 4; ++i)
		//medianHarmonics[i] = median(workerBuffer.data(), sampleBuffer.data() + (Lh * i), Lh);
	//////////////////////////////////////////////////

	//////////////////////////////////////////////////
	// Percussive
	//const int percChunks = SPECTRUM_REAL_SIZE / PERCUSSIVE_MEDIAN_LEN;
	//const int percTail = SPECTRUM_REAL_SIZE % PERCUSSIVE_MEDIAN_LEN;
	{
		


	}

	//////////////////////////////////////////////////
	// Harmonics
	// {
	// }


	//////////////////////////////////////////////////

}
#endif


void fft(CArray& x)
{
	const size_t N = x.size();
	if (N <= 1) return;

	// divide
	CArray even = x[std::slice(0, N / 2, 2)];
	CArray  odd = x[std::slice(1, N / 2, 2)];

	// conquer
	fft(even);
	fft(odd);

	// combine
	for (size_t k = 0; k < N / 2; ++k)
	{
		Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
		x[k] = even[k] + t;
		x[k + N / 2] = even[k] - t;
	}
}

void ifft(CArray& x)
{
	// conjugate the complex numbers
	x = x.apply(std::conj);

	// forward fft
	fft(x);

	// conjugate the complex numbers again
	x = x.apply(std::conj);

	// scale the numbers
	x /= x.size();
}


