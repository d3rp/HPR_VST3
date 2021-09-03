#pragma once

#include "IPlug_include_in_plug_hdr.h"
#include "IControls.h"
#include "ISender.h"

#include <complex>
#include <iostream>
#include <valarray>
#include <cassert>

//const double PI = 3.141592653589793238460;


using namespace iplug;
using namespace igraphics;

typedef std::complex<sample> Complex;
typedef std::valarray<Complex> CArray;

void fft(CArray& x);
void ifft(CArray& x);

#define MAX_CHUNK_SIZE 1024

struct Chunk
{
	sample const* const* const data;
	const int nFrames;
};

struct CircularBuffer
{

	CircularBuffer(const int n)
		: data(0.0, n)
		, size(n)
		, index(0)
	{
	}

	// void push(sample x)
	// {
	// 	data[index++] = x;
	// 	if (index == size)
	// 		index = 0;
	// }

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

// template <typename FType = double, int Size = SPECTRUM_REAL_SIZE>
// struct CircularBufferPrealloc
// {
// 
// 	CircularBufferPrealloc()
// 	{
// 		for (auto i = 0; i < Size; ++i)
// 			data[i] = 0.0;
// 	}
// 
// 	void push(FType x)
// 	{
// 		data[index++] = x;
// 		if (index == size)
// 			index = 0;
// 	}
// 
// 	void push(FType const * const newData, const int n) noexcept
// 	{
// 		int i = 0;
// 		while (i < n)
// 		{
// 			while (index < size && i < n)
// 				data[index++] = newData[i++];
// 
// 			if (index == size)
// 				index = 0;
// 		}
// 	}
// 
// 	FType get(int i)
// 	{
// 		assertm(i < (2 * size), "CircularBuffer[index] -> index is out of bounds even after offset!");
// 		const int offsetIndex = i < size ? index + i : i - size;
// 		return data[offsetIndex];
// 	}
// 
// 	FType& operator[](int i) { return get(i); }
// 	const FType& operator[](int i) const { return get(i); }
// 
// 	std::array<FType, Size> data;
// 	int index = 0;
// };
// 

// /**
// * Twin circularbuffers that can be toggled to produce OLA
// */
// struct ToggleBuffer
// {
// 	using FType = double;
// 	ToggleBuffer(const int n)
// 		:buffers{ CircularBuffer<FType>(n), CircularBuffer<FType>(n) }
// 	{ }
// 
// 	void push(FType&& x) { buffers[index].push(x); }
// 
// 	void push(FType const * const newData, const int n) noexcept
// 	{
// 		buffers[index].push(newData, n);
// 	}
// 
// 	std::valarray<FType> getShifted() const noexcept
// 	{
// 		return buffers[index].getShifted();
// 	}
// 
// 	void getShifted(FType* const dst, int n) noexcept
// 	{
// 		buffers[index].getShifted(dst, n);
// 	}
// 
// 	void flip() { index = (int)!bIdx; }
// 
// 	CircularBuffer<FType> buffers[2];
// 	bool bIdx = 0;  // hack
// 	unsigned int index = 0;
// };

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

#define HARMONICS_MEDIAN_LEN 3
#define PERCUSSIVE_MEDIAN_LEN 16
static constexpr unsigned int SPECTRUM_SIZE = PLUG_LATENCY;
//static constexpr unsigned int SPECTRUM_REAL_SIZE_D = SPECTRUM_SIZE * 0.5;
static constexpr unsigned int Lh = (SPECTRUM_SIZE / 8); // half is positive frequencies and 
static constexpr unsigned int SPECTRUM_REAL_SIZE = SPECTRUM_SIZE; // (int)SPECTRUM_REAL_SIZE_D;
static constexpr unsigned int HARMONICS_MATRIX_SIZE = HARMONICS_MEDIAN_LEN * SPECTRUM_REAL_SIZE;


sample median(sample * const srcDst, const int n);

struct HarmonicsMatrix
{
	HarmonicsMatrix()
	{

		for (auto& e : data)
			e = 0.0;

		for (auto& e : indices)
			e = 0;

		for (auto& e : tmp)
			e = 0;

		// for (auto i = 0; i < HARMONICS_MATRIX_SIZE; ++i)
		// 	data[i] = 0.0;

		// for (auto i = 0; i < SPECTRUM_REAL_SIZE; ++i)
		// 	indices[i] = 0.0;
	}

	// void push_back(const double x);
	// void dumpMedians(sample* const dst) const;
	void nextMedians(sample* const dst, sample* const worker, sample const* const src, int n)
	{
		// update values
		for (auto i = 0; i < n; ++i)
		{
			data[indices[i]] = src[i];
			advanceIndices(i);
		}

		//const int L = HARMONICS_MEDIAN_LEN;
		
		for (auto i = 0; i < n; ++i)
		{
			// Copy to tmp buffer for to take median
			int j = 0;
			for (auto& e : tmp)
			{
				e = data[indices[i]];
				advanceIndices(i);
			}
			// save median value
			dst[i] = median(tmp.data(), n);
		}

	}

	void advanceIndices(int i)
	{
		if (indices[i] + 1 < HARMONICS_MEDIAN_LEN)
			++indices[i];
		else
			indices[i] = 0;
	}

	std::array<unsigned int, SPECTRUM_REAL_SIZE> indices;
	std::array<double, HARMONICS_MATRIX_SIZE> data;
	std::array<double, HARMONICS_MEDIAN_LEN> tmp;
};


class hpr final : public Plugin
{
public:
	hpr(const InstanceInfo& info);

#if IPLUG_DSP // http://bit.ly/2S64BDd
	void OnIdle() override;
	void ProcessBlock(sample** inputs, sample** outputs, int nFrames) override;
#endif

	bool initted = false;
	size_t wi = 0;
	std::valarray<Complex> complexBuffer;
	std::array<sample, Lh> medianHarmonics;
	std::vector<sample> sampleBuffer;
	std::vector<sample> workerBuffer;
	std::vector<sample> workerBuffer2;
	std::vector<sample> previousBuffer;
	std::vector<sample> nextOutputs;
	CircularBuffer olaBuffer;
	std::vector<sample> start;
	std::vector<sample> end;



	HarmonicsMatrix harmonicsMatrix;
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
