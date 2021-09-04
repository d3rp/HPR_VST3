#pragma once

#include "IPlug_include_in_plug_hdr.h"
#include "IControls.h"
#include "ISender.h"

#include <complex>
#include <iostream>
#include <valarray>
#include <cassert>

using namespace iplug;
using namespace igraphics;

typedef std::complex<sample> Complex;
typedef std::valarray<Complex> CArray;

void fft(CArray& x);
void ifft(CArray& x);

struct CircularBuffer
{

	CircularBuffer(const int n)
		: data(0.0, n)
		, size(n)
		, index(0)
	{
	}

	void push(sample const * const newData, const int n) noexcept
	{
		int i = 0;
		while (i < n)
		{
			while (index < size && i < n)
				data[index++] = newData[i++];

			if (index == size)
				index = 0;
		}
	}

	sample get(int i) const
	{
		// CircularBuffer[index] -> index is out of bounds even after offset!
		assert(i < (2 * size));

		const int offsetIndex = (index + i) < size ? index + i : i - size;
		return data[offsetIndex];
	}

	//sample& operator[](int i) { return get(i); }
	const sample& operator[](int i) const { return get(i); }

	void getShifted(sample* const dst, int n) const noexcept
	{
		// Oh no, should've kept books about buffer sizes..
		assert(n <= size);

		for (auto i = 0; i < size; ++i)
			dst[i] = get(i);
	}

	const unsigned int size;
	std::valarray<sample> data;
	int index;
};

const int kNumPresets = 1;
enum EParams
{
	kGain = 0,
	kH = 1,
	kP = 2,
	kR = 3,
	kNumParams
};

#define HARMONICS_MEDIAN_LEN 16
#define PERCUSSIVE_MEDIAN_LEN 16
static constexpr unsigned int FFT_SIZE = PLUG_LATENCY;
static constexpr unsigned int HARMONICS_MATRIX_SIZE = HARMONICS_MEDIAN_LEN * FFT_SIZE;


sample median(sample * const srcDst, const int n);

/**
Handles buffering for the horisontal median filter
*/
struct HarmonicsMatrix
{
	static constexpr int L = HARMONICS_MEDIAN_LEN;

	HarmonicsMatrix()
		: data()
		, tmp()
	{
		for (auto i = 0; i < HARMONICS_MATRIX_SIZE; ++i)
			data[i] = 0.0;

		for (auto i = 0; i < HARMONICS_MEDIAN_LEN; ++i)
			tmp[i] = 0.0;
	}

	void nextMedians(sample* const dst, sample* const worker, sample const* const src, int n)
	{
		// update values
		for (auto i = 0; i < n; ++i)
			data[(L * i) + idx] = std::max(src[i], 1.0);

		advanceIndices();

		if (!enoughHistory)
		{
			std::copy_n(src, n, dst);
			enoughHistory = (idx > L - 2);
		}

		for (auto i = 0; i < n; ++i)
		{
			// Copy to tmp buffer for to take median
			int j = 0;
			for (auto& e : tmp)
				e = data[(L * i) + j++];

			// save median value
			dst[i] = median(tmp.data(), HARMONICS_MEDIAN_LEN);
		}
	}

	void advanceIndices()
	{
		if (++idx < HARMONICS_MEDIAN_LEN)
			;
		else
			idx = 0;
	}

	bool enoughHistory = false;
	int idx = 0;
	std::array<sample, HARMONICS_MATRIX_SIZE> data;
	std::array<sample, HARMONICS_MEDIAN_LEN> tmp;
};

void fromReal(std::valarray<Complex>& dst, sample const* const src, const int n);
void toReal(sample* dst, std::valarray<Complex> const& src, const int n);

class hpr final : public Plugin
{
public:
	hpr(const InstanceInfo& info);

#if IPLUG_DSP // http://bit.ly/2S64BDd
	void OnIdle() override;
	void ProcessBlock(sample** inputs, sample** outputs, int nFrames) override;
#endif

	/**
	Calculate masks for H, P, R and apply to samples.
	Note: Also transforms from frequency to time-domain, because the proposal shows
	combining the signals in the time-domain (multiplication of parameter multipliers in
	frequency-domain would equate to convolution in time-domain).

	Implicit buffers are
		 harmonics	= median filtered harmonics
		 percussive = median filtered percussive values
		 H, P, R	= parameter selectors (multipliers applied to masks)
		 sampleBuffer	= output after applying masks and params
	*/
	void calculateHPR(int windowL) noexcept
	{
		// Implementation based on the paper by Driedger et al.: 
		// https://www.researchgate.net/publication/303667409_Extending_Harmonic-Percussive_Separation_of_Audio_Signals
		const sample e = std::numeric_limits<sample>::min();

		auto& harmonicsTmp = workerBuffer;
		std::copy_n(harmonics.begin(), windowL, harmonicsTmp.begin());

		auto& percussiveTmp = workerBuffer2;
		std::copy_n(percussive.begin(), windowL, percussiveTmp.begin());

		auto percIt = percussive.begin();
		std::for_each_n(harmonicsTmp.begin(), windowL, [&](sample h) {
			const auto p = *(percIt++);
			if (h / (p + e) > 1.0)
				return h;
			else
				return 0.0;
			});

		auto harmIt = harmonics.begin();
		std::for_each_n(percussiveTmp.begin(), windowL, [&](sample p) {
			const auto h = *(harmIt++);
			if (p / (h + e) >= 1.0)
				return p;
			else
				return 0.0;
			});

		const auto maxH = std::max(1.0, *std::max_element(harmonicsTmp.begin(), harmonicsTmp.begin() + windowL));
		const auto maxP = std::max(1.0, *std::max_element(percussiveTmp.begin(), percussiveTmp.begin() + windowL));
		const auto maxS = *std::max_element(sampleBuffer.begin(), sampleBuffer.begin() + windowL);
		std::for_each_n(sampleBuffer.begin(), windowL, [&](sample s) { s / (maxS + e); });

		for (int i = 0; i < windowL; ++i)
		{
			const auto h = harmonicsTmp[i] / maxH;
			const auto p = percussiveTmp[i] / maxP;

			if (h + p == 0)
				continue;

			const auto r = std::max(0.0, 1.0 - (h + p));

			assert(h > 0); assert(p > 0); assert(r > 0);
			const auto s = sampleBuffer[i];
			
			const auto xh = s * h;
			const auto xp = s * p;
			const auto xr = s * r;

			assert(xh > 0); assert(xp > 0); assert(xr > 0);
			assert(xh < 1.0); assert(xp < 1.0); assert(xr < 1.0);

			harmonics[i] = xh;
			percussive[i] = xp;
			residual[i] = xr;
		}
			// Time-domain
		ifftBuffer(harmonics, windowL);
		ifftBuffer(percussive, windowL);
		ifftBuffer(residual, windowL);

		std::for_each_n(harmonics.begin(), windowL, [mult=H](sample x) { return x * mult; });
		std::for_each_n(percussive.begin(), windowL, [mult=P](sample x) { return x * mult; });
		std::for_each_n(residual.begin(), windowL, [mult=R](sample x) { return x * mult; });
	
		assert(*std::max_element(harmonics.begin(), harmonics.begin() + windowL) < 1.0);
		assert(*std::max_element(percussive.begin(), percussive.begin() + windowL) < 1.0);
		assert(*std::max_element(residual.begin(), residual.begin() + windowL) < 1.0);

		for (int i = 0; i < windowL; ++i)
			sampleBuffer[i] = harmonics[i] + percussive[i] + residual[i];

		assert(*std::max_element(sampleBuffer.begin(), sampleBuffer.begin() + windowL) < 1.0);;
	}
	void fftBuffer(std::vector<sample>& buffer, int windowL) noexcept
	{
		fftBuffer(buffer.data(), windowL);
	}
	void fftBuffer(sample* const srcDst, int windowL) noexcept
	{
		fromReal(complexBuffer, srcDst, windowL);
		fft(complexBuffer);
		toReal(srcDst, complexBuffer, windowL);
	}

	void ifftBuffer(std::vector<sample>& buffer, int windowL) noexcept
	{
		ifftBuffer(buffer.data(), windowL);
	}
	void ifftBuffer(sample* const srcDst, int windowL) noexcept
	{
		fromReal(complexBuffer, srcDst, windowL);
		ifft(complexBuffer);
		toReal(srcDst, complexBuffer, windowL);
	}

	// Parameter multipliers to be updated from the GUI
	double H, P, R;

	bool initted = false;
	size_t wi = 0;

	std::valarray<Complex> complexBuffer;
	std::vector<sample> sampleBuffer;
	std::vector<sample> workerBuffer;
	std::vector<sample> workerBuffer2;
	std::vector<sample> nextOutput;
	std::vector<sample> start;

	std::vector<sample> harmonics;
	std::vector<sample> percussive;
	std::vector<sample> residual;

	HarmonicsMatrix harmonicsMatrix;
};
