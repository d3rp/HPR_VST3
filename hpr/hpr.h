#pragma once

#include "IPlug_include_in_plug_hdr.h"
#include "IControls.h"
#include "ISender.h"

#include <complex>
#include <iostream>
#include <valarray>

//const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

BEGIN_IPLUG_NAMESPACE
void fft(CArray& x);
void ifft(CArray& x);

#define MAX_CHUNK_SIZE 1024

struct Chunk
{
	sample const* const* const data;
	const int nFrames;
};

template <typename FloatingPointType = double>
struct CircularBuffer
{

	CircularBuffer(const int n)
		: data(0.0, n)
		, size(n)
		, index(0)
	{
	}

	void push(FloatingPointType x)
	{
		data[index++] = x;
		if (index == size)
			index = 0;
	}

	void push(FloatingPointType const * const newData, const int n) noexcept
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


	/**
	* Works mainly for non-real-time applications (allocation)
	*/
	std::valarray<FloatingPointType> getShifted() const noexcept
	{
		return data.cshift(index);
	}

	const unsigned int size;
	std::valarray<FloatingPointType> data;
	int index;
};

/**
* Twin circularbuffers that can be toggled to produce OLA
*/
struct ToggleBuffer
{
	using FType = double;
	ToggleBuffer(const int n)
		:buffers{ CircularBuffer<FType>(n), CircularBuffer<FType>(n) }
	{ }

	void push(FType&& x) { buffers[index].push(x); }

	void push(FType const * const newData, const int n) noexcept
	{
		buffers[index].push(newData, n);
	}

	std::valarray<FType> getShifted() const noexcept
	{
		return buffers[index].getShifted();

	}


	void flip() { index = (int)!bIdx; }

	CircularBuffer<FType> buffers[2];
	bool bIdx = 0;  // hack
	unsigned int index = 0;
};
END_IPLUG_NAMESPACE

const int kNumPresets = 1;
enum EParams
{
	kGain = 0,
	kH = 1,
	kP = 2,
	kR = 3,
	kNumParams
};

enum EControlTags
{
	kCtrlTagSpectrum = 0,
	kCtrlTagDebug,
	kCtrlTagDebug2
};

#define HARMONICS_MEDIAN_LEN 1
#define PERCUSSIVE_MEDIAN_LEN 200
#define SPECTRUM_SIZE 2048
#define SPECTRUM_REAL_SIZE_D SPECTRUM_SIZE * 0.5
#define Lh (SPECTRUM_SIZE / 8) // half is positive frequencies and 
static constexpr unsigned int SPECTRUM_REAL_SIZE = (int)SPECTRUM_REAL_SIZE_D;

using namespace iplug;
using namespace igraphics;

class hpr final : public Plugin
{
public:
	hpr(const InstanceInfo& info);

#if IPLUG_DSP // http://bit.ly/2S64BDd
	void OnIdle() override;
	void ProcessBlock(sample** inputs, sample** outputs, int nFrames) override;
#endif

	std::valarray<Complex> workerComplexBuffer;
	std::array<sample, Lh> medianHarmonics;
	std::vector<sample> sampleBuffer;
	std::vector<sample> workerBuffer;
	std::vector<sample> previousBuffer;
	std::vector<sample> nextOutputs;
	ToggleBuffer olaBuffer;


	std::vector<sample> harmonics;
	std::vector<sample> percussive;
	std::vector<sample> residual;

	//CircularBuffer<sample> lambdaDrawBuffer;

	//IBufferSender<1, 4, SPECTRUM_SIZE/2> spectrumSender;
	//ISender<1> textSender;
	//ISender<1> textSender2;
	//ISenderData<1> mLastOutputData = { kCtrlTagDebug, 1, 0 };
	//ISenderData<1> mLastOutputData2 = { kCtrlTagDebug2, 1, 0 };
};
