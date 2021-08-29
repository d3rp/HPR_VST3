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

sample median(sample* const worker, sample const* const src, const int n)
{
	std::copy_n(src, n, worker);
	auto m = n / 2;
	std::nth_element(worker, worker + m, worker + n);
	return worker[m];
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
	, workerComplexBuffer(Complex(), SPECTRUM_REAL_SIZE)
	, olaBuffer(SPECTRUM_REAL_SIZE * 0.5)
	//, lambdaDrawBuffer(SPECTRUM_REAL_SIZE * 64)
{

	//GetParam(kGain)->InitDouble("Gain", 50.0, 0.0, 100.0, 0.01, "%");
	GetParam(kH)->InitDouble("H", 0., 0., 100.0, 0.01, "%");
	GetParam(kP)->InitDouble("P", 0., 0., 100.0, 0.01, "%");
	GetParam(kR)->InitDouble("R", 0., 0., 100.0, 0.01, "%");

	sampleBuffer.resize(SPECTRUM_SIZE, -1.0);
	workerBuffer.resize(SPECTRUM_SIZE, -1.0);
	previousBuffer.resize(SPECTRUM_SIZE, -1.0);
	nextOutputs.resize(SPECTRUM_SIZE, -1.0);

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

		//const IRECT right = pGraphics->GetBounds().GetReducedFromLeft(left.W()).GetReducedFromRight(10);
		//const IRECT rightMost = pGraphics->GetBounds().GetFromRight(10); // GetReducedFromLeft(left.W()).GetReducedFromRight(10);

		//pGraphics->AttachControl(new ITextControl(nextCell(), "HPR", IText(50, getColor(Colors::FG2))));
		//pGraphics->AttachControl(new IRTTextControl<1, float>(nextCell(), "window: %.0f", ":", "IRTTextControl", IText(20, getColor(Colors::FG2))), kCtrlTagDebug);
		//pGraphics->AttachControl(new IRTTextControl<1, float>(nextCell(), "max amp: %.0f", ":", "IRTTextControl2", IText(20, getColor(Colors::FG2))), kCtrlTagDebug2);
		//pGraphics->AttachControl(new IVKnobControl(nextCell(), kGain, "", knobStyle()));
		pGraphics->AttachControl(new IVKnobControl(nextCell(), kH, "", knobStyle()));
		pGraphics->AttachControl(new IVKnobControl(nextCell(), kP, "", knobStyle()));
		pGraphics->AttachControl(new IVKnobControl(nextCell(), kR, "", knobStyle()));
		//pGraphics->AttachControl(new IVKnobControl(b.GetCentredInside(100).GetVShifted(-100).GetHShifted(100), kDither));
		//pGraphics->AttachControl(new Spectrograph(b.GetCentredInside(100).GetVShifted(-100).GetHShifted(100)), kCtrlTagSpectrum);
	/*    pGraphics->AttachControl(new IVPlotControl(nextCell(),{
		  { COLOR_WHITE, [](double d) { return d;  } } },
		 1024, "fooPlot"), kCtrlTagSpectrum);
	*/

	//  const int realSpectrumSize = SPECTRUM_REAL_SIZE;
	//  const int spectrogramW = (lambdaDrawBuffer.size / realSpectrumSize);
	//  const float spectrogramScale = right.H() / realSpectrumSize;
	//  const float spectrogramLineW = right.W() / spectrogramW;
	//  pGraphics->AttachControl(new ILambdaControl(right,
	//      [
	//          &lambdaDrawBuffer = lambdaDrawBuffer,
	//          realSpectrumSize = realSpectrumSize,
	//          spectrogramW = spectrogramW,
	//          scale = spectrogramScale,
	//          lineW = spectrogramLineW
	//      ](ILambdaControl* pCaller, IGraphics& g, IRECT& r) {
	//          auto b = lambdaDrawBuffer.getShifted();
	//          const IColor lineColor(255, 255, 255, 255);

	//          g.PathLine(r.L, r.T, r.L, r.T);
	//          //g.PathLineTo(r.L, r.T);
	//          for (auto j = 0; j < spectrogramW; ++j)
	//          {
	//              const float x = r.L + lineW * j;
	//              for (auto i = 1; i < realSpectrumSize; ++i)
	//              {
	//                  const float y0 = r.B - scale * (i - 1);
	//                  const float y1 = r.B - scale * i;
	//                  const float ampl = b[i + j * realSpectrumSize];
	//                  g.PathRect(IRECT(x, y0, x + lineW, y1));
	//                  g.PathFill(IColor(ampl * 255, 255, 255, 255));
	//              }
	//          }
	//          
	//          // bounds from r
	//          // pCaller->GetAnimationProgress() for interpolation steps
	//          // pCaller-> mouse state from pCaller etc
	//          // g.PathTriangle(...), g.PathFill(color), g.PathFill(...), 

	//      }, 1, false, true));
		  //pGraphics->AttachControl(new IVScopeControl<1, SPECTRUM_REAL_SIZE>(rightMost, "",
		//DEFAULT_STYLE
		//  .WithColor(kBG, getColor(Colors::BG))
		// .WithColor(kFG, COLOR_WHITE)), kCtrlTagSpectrum);
	};
#endif
	testMedian();
}

void medianHorizontal(sample* dst, sample* src, int nFrames)
{
	return;
}

void medianVertical(sample* dst, sample* src, int nFrames)
{
	return;
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


void normaliseForLogScale(sample* const srcDst, sample* const dampingBuffer, int nFrames)
{
	constexpr int halfSpectrum = SPECTRUM_REAL_SIZE;
	constexpr double ratio = 22000.0 / halfSpectrum;
	static const double highF = std::log10(22000);
	static const double lowF = std::log10(100);
	static const double fRange = highF - lowF;
	for (int s = 0; s < halfSpectrum; ++s)
	{
		int idx = std::log2(1.0 + s / (double)halfSpectrum) * halfSpectrum;
		//int idx = std::log10(1.0 + ((s * ratio) / 261.0)) * 261; // wikipedia frequency log scale
		//int idx = 511 * (std::log10(1.0 + s * ratio) - lowF) / fRange; // stackoverflow copypaste
		sample newS = (srcDst[s] / 50.0) - 1.0;
		sample oldS = dampingBuffer[idx];
		srcDst[idx] = std::max(newS, 0.9 * (oldS + 1.0) - 1.0);
	}
}
void spectrogram(CircularBuffer<sample>* dst, std::valarray<Complex>& const src, const int n) noexcept
{
	// takes fft'd buffer, calc energy in other buffer
	for (auto i = 0; i < n; ++i)
	{
		const auto x = src[i];
		const auto real = x.real();
		const auto im = x.imag();
		dst->push(10 * std::log(real * real + im * im));
	}
}
void copyFromReal(std::valarray<Complex>& dst, sample const* const src, const int n)
{
	for (int s = 0; s < n; s++)
		dst[s] = Complex(src[s], 0);
}
void copyToReal(sample* dst, std::valarray<Complex> const& src, const int n)
{
	for (int s = 0; s < n; s++)
		dst[s] = src[s].real();
}

void process()
{
	// triggered when we have full spectrum worth of fresh data
}

void fromMono(sample** dst, sample* src, const int nChans, const int n)
{
	for (int c = 0; c < nChans; ++c) {
		std::copy_n(src, n, dst[c]);
	}

}
void toMono(sample* dst, sample** src, const int nChans, const int n)
{
	int r = 1.0 / (float)nChans;
	auto combine = [&r = r](const double x) { return x * r;  };

	// First channel
	std::copy_n(src[0], n, dst);
	std::for_each_n(dst, n, combine);

	// Consecutive channels
	int i = 0;
	for (int c = 1; c < nChans; ++c) {
		i = 0;
		sample* chanSrc = src[c];
		std::for_each_n(dst, n, [&i = i, &r = r, &chanSrc = chanSrc](const double x) {
			return x + chanSrc[i++] * r;
			});
	}
}
// TODO : NEXT : OLA
void hpr::ProcessBlock(sample** inputs, sample** outputs, int nFrames)
{
	const double gain = GetParam(kGain)->Value() / 100.;
	const int nChans = NOutChansConnected();

	std::fill_n(sampleBuffer.begin(), nFrames, 0.0);

	fromMono(outputs, previousBuffer.data(), nChans, nFrames);
	toMono(previousBuffer.data(), inputs, nChans, nFrames);
	return;

	// zero
	// combine to mono
	constexpr int halfSSize = SPECTRUM_REAL_SIZE * 0.5;
	toMono(workerBuffer.data(), inputs, 2, halfSSize);

	// TODO: consider mono or stereo?

	hann(sampleBuffer.data(), workerBuffer.data(), halfSSize);
	copyFromReal(workerComplexBuffer, workerBuffer.data(), halfSSize);

	fft(workerComplexBuffer);
	// magic
	ifft(workerComplexBuffer);
	copyToReal(workerBuffer.data(), workerComplexBuffer, halfSSize);

	// first half of OLA
	// on first run, it's initialised as zero
	olaBuffer.push(workerBuffer.data(), halfSSize);
	auto ob = olaBuffer.getShifted();




	//////////////////////////////////////////////////
	// pass 1/4
	hann(sampleBuffer.data(), sampleBuffer.data(), nFrames);
	for (int s = 0; s < nFrames; s++)
		workerComplexBuffer[s] = Complex(sampleBuffer[s], 0);

	fft(workerComplexBuffer);

	//spectrogram(&lambdaDrawBuffer, workerComplexBuffer, nFrames/2);

	//const auto Lh = SPECTRUM_SIZE / (2 * 4);

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
		const int L = PERCUSSIVE_MEDIAN_LEN;
		const int halfP = PERCUSSIVE_MEDIAN_LEN / 2;
		const int N = SPECTRUM_REAL_SIZE - halfP - 1;
		for (auto i = halfP; i < N; ++i)
		{
			const auto* s = &sampleBuffer[i - halfP];
			sampleBuffer[i] = median(workerBuffer.data(), s, L);
		}



	}


	//////////////////////////////////////////////////

	copyFromReal(workerComplexBuffer, sampleBuffer.data(), SPECTRUM_REAL_SIZE);
	ifft(workerComplexBuffer);

	copyToReal(outputs[0], workerComplexBuffer, nFrames);
	copyToReal(outputs[1], workerComplexBuffer, nFrames);

	// end of pass 1/4
	//////////////////////////////////////////////////

	// X = X0 + (X1 - X0)(log(V) - log(V0)) / (log(V1) - log(V0))


	//normaliseForLogScale(sampleBuffer.data(), previousBuffer.data(), nFrames);

	//medianHorizontal(workerBuffer, harmonics.data(), nFrames);
	//medianVertical(*inputs, percussive.data(), nFrames);
	//calculateResidual(*inputs, harmonics.data(), percussive.data(), residual.data(), nFrames);

	/**
	Senders for different components

	*/
	//sample* pData = sampleBuffer.data();
	//spectrumSender.ProcessBlock(&pData, nFrames/2, kCtrlTagSpectrum, 1);

	//auto currentMax = std::max_element(sampleBuffer.begin(), sampleBuffer.begin() + nFrames/2);
	//mLastOutputData.vals[0] = (float)nFrames;
	//textSender.PushData(mLastOutputData);

	//mLastOutputData2.vals[0] = 100;// (float)currentMax;
	//textSender2.PushData(mLastOutputData2);
}
#endif


BEGIN_IPLUG_NAMESPACE
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
END_IPLUG_NAMESPACE


