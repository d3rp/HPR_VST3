#pragma once
#include "hpr.h"

// Alternative from c++ cookbook (and original from a newsgroup)	
// based on public domain code that can be found on the digital signal processing newswgroup on usenet (comp.dsp). 

// TODO : Potentially much faster, as the current fft allocates
// and in profiling proves to be the major bottleneck

// #include <iostream>
// #include <complex>
// #include <cmath>
// #include <iterator>
// 
// using namespace std;
// 
// unsigned int bitReverse(unsigned int x, int log2n) {
// 	int n = 0;
// 	int mask = 0x1;
// 	for (int i = 0; i < log2n; i++) {
// 		n <<= 1;
// 		n |= (x & 1);
// 		x >>= 1;
// 	}
// 	return n;
// }
// 
// const double PI = 3.1415926536;
// 
// template<class Iter_T>
// void fft(Iter_T a, Iter_T b, int log2n)
// {
// 	typedef typename iterator_traits<Iter_T>::value_type complex;
// 	const complex J(0, 1);
// 	int n = 1 << log2n;
// 	for (unsigned int i = 0; i < n; ++i) {
// 		b[bitReverse(i, log2n)] = a[i];
// 	}
// 	for (int s = 1; s <= log2n; ++s) {
// 		int m = 1 << s;
// 		int m2 = m >> 1;
// 		complex w(1, 0);
// 		complex wm = exp(-J * (PI / m2));
// 		for (int j = 0; j < m2; ++j) {
// 			for (int k = j; k < n; k += m) {
// 				complex t = w * b[k + m2];
// 				complex u = b[k];
// 				b[k] = u + t;
// 				b[k + m2] = u - t;
// 			}
// 			w *= wm;
// 		}
// 	}
// }
// 
// int main() {
// 	typedef complex<double> cx;
// 	cx a[] = { cx(0,0), cx(1,1), cx(3,3), cx(4,4),
// 	  cx(4, 4), cx(3, 3), cx(1,1), cx(0,0) };
// 	cx b[8];
// 	fft(a, b, 3);
// 	for (int i = 0; i < 8; ++i)
// 		cout << b[i] << "\n";
// }

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



void runTests()
{
	testMedian();

}

void medianPercussive(sample* const dst, sample* const worker, sample const* const src, int n)
{
	const int L = PERCUSSIVE_MEDIAN_LEN;
	const int halfP = L / 2;
	const int N = FFT_SIZE - halfP - 1;
	for (auto i = 0; i < halfP; ++i)
		dst[i] = src[i];

	for (auto i = N - halfP; i < N; ++i)
		dst[i] = src[i];

	for (auto i = halfP; i < (N - halfP); ++i)
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

void hpr::OnIdle()
{
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

#define BETA 0.01
void addWithParams(sample* const dst, sample* harm, sample* perc, const double H, const double P, const double R, int n)
{
	const double total = H + P + R;
	const double scale = std::max(0.49, total == 0.0 ? 0.0 : 0.5 / total);


	std::fill_n(dst, n, 0.0);
	return;

	for (auto i = 0; i < n; ++i)
	{
		const sample h = harm[i];
		if (std::abs(h) > BETA)
			dst[i] = h * H;
		else
			dst[i] = h * R;
	}

	for (auto i = 0; i < n; ++i)
	{
		const sample p = perc[i];
		if (std::abs(p) > BETA)
			dst[i] += p * P;
		else
			dst[i] += p * R;
	}

	for (auto i = 0; i < n; ++i)
		dst[i] *= scale;
}

