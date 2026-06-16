//
//  fcwt.cpp
//  fCWT
//
//  Created by Lukas Arts on 21/12/2020.
//  Copyright © 2021 Lukas Arts.
/*Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

//
//  fcwt.cpp
//  fCWT-testing
//
//  Created by Lukas Arts on 21/12/2020.
//  Copyright © 2020 Lukas Arts.
/*Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include "api.h"

#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <mutex>
#include <vector>

#include <kfr/dft.hpp>

namespace {
using KfrComplex = kfr::complex<float>;

KfrComplex *as_kfr(std::complex<float> *data) {
  return reinterpret_cast<KfrComplex *>(data);
}

const KfrComplex *as_kfr(const std::complex<float> *data) {
  return reinterpret_cast<const KfrComplex *>(data);
}
} // namespace

struct fcwt::FFTPlan {
  explicit FFTPlan(int size) : inverse(static_cast<size_t>(size)) {
    temp.resize(inverse.temp_size);
  }

  void execute_inverse(std::complex<float> *out,
                       const std::complex<float> *in) const {
    inverse.execute(as_kfr(out), as_kfr(in), const_cast<kfr::u8 *>(temp.data()),
                    true);
  }

  kfr::dft_plan<float> inverse;
  std::vector<kfr::u8> temp;
};

namespace {
void execute_forward_complex(const std::complex<float> *input,
                             std::complex<float> *output, int size) {
  const kfr::dft_plan<float> plan(static_cast<size_t>(size));
  std::vector<kfr::u8> temp(plan.temp_size);
  plan.execute(as_kfr(output), as_kfr(input), temp.data(), false);
}

void execute_forward_real(const float *input, std::complex<float> *output,
                          int size) {
  const kfr::dft_plan_real<float> plan(static_cast<size_t>(size));
  std::vector<kfr::u8> temp(plan.temp_size);
  plan.execute(as_kfr(output), input, temp.data());
}
} // namespace

void fcwt::API::daughter_wavelet_multiplication(std::complex<float> *input,
                                                std::complex<float> *output,
                                                float const *mother,
                                                const float scale, int isize,
                                                bool imaginary,
                                                bool doublesided) const {
  const auto isizef = static_cast<float>(isize);
  const float endpointf = std::min(isizef / 2.0f, ((isizef * 2.0f / scale)));
  const float step = static_cast<float>(scale) / 2.0f;
  int endpoint = static_cast<int>(endpointf);
  const int athreads = std::min(threads, std::max(1, endpoint / 16));
  const float maximum = isizef - 1;
  const int s1 = isize - 1;

#ifdef AVX
  const __m256 step8 = _mm256_set1_ps(step);
  const __m256 offset = _mm256_set_ps(3, 3, 2, 2, 1, 1, 0, 0);
  const __m256 maximum8 = _mm256_set1_ps(maximum);
  const float imaginary_sign = static_cast<float>(1 - 2 * imaginary);

#ifndef SINGLE_THREAD
#pragma omp parallel for
#endif
  for (int i = 0; i < athreads; i++) {
    const int start = (endpoint * i) / athreads;
    const int end = (endpoint * (i + 1)) / athreads;

    int q1 = start;
    for (; q1 + 3 < end; q1 += 4) {
      const auto q = static_cast<float>(q1);
      const __m256 qq = _mm256_set1_ps(q);
      const U256f tmp = {_mm256_min_ps(
          maximum8, _mm256_mul_ps(step8, _mm256_add_ps(qq, offset)))};

      const float m0 = mother[static_cast<int>(tmp.a[0])];
      const float m1 = mother[static_cast<int>(tmp.a[2])];
      const float m2 = mother[static_cast<int>(tmp.a[4])];
      const float m3 = mother[static_cast<int>(tmp.a[6])];
      const __m256 wav =
          _mm256_set_ps(m3 * imaginary_sign, m3, m2 * imaginary_sign, m2,
                        m1 * imaginary_sign, m1, m0 * imaginary_sign, m0);

      const auto *in = reinterpret_cast<const float *>(input + q1);
      auto *out = reinterpret_cast<float *>(output + q1);
      _mm256_storeu_ps(out, _mm256_mul_ps(_mm256_loadu_ps(in), wav));
    }

    for (; q1 < end; q1++) {
      const float tmp = std::min(maximum, step * static_cast<float>(q1));
      output[q1].real(input[q1].real() * mother[static_cast<int>(tmp)]);
      output[q1].imag(input[q1].imag() * mother[static_cast<int>(tmp)] *
                      imaginary_sign);
    }

    if (doublesided) {
      q1 = start;
      for (; q1 + 3 < end; q1 += 4) {
        const auto q = static_cast<float>(q1);
        const __m256 qq = _mm256_set1_ps(q);
        const U256f tmp = {_mm256_min_ps(
            maximum8, _mm256_mul_ps(step8, _mm256_add_ps(qq, offset)))};

        const float m0 = mother[static_cast<int>(tmp.a[0])];
        const float m1 = mother[static_cast<int>(tmp.a[2])];
        const float m2 = mother[static_cast<int>(tmp.a[4])];
        const float m3 = mother[static_cast<int>(tmp.a[6])];
        const __m256 wav =
            _mm256_set_ps(m0, m0 * imaginary_sign, m1, m1 * imaginary_sign, m2,
                          m2 * imaginary_sign, m3, m3 * imaginary_sign);

        const int dest = s1 - q1 - 3;
        const auto *in = reinterpret_cast<const float *>(input + dest);
        auto *out = reinterpret_cast<float *>(output + dest);
        _mm256_storeu_ps(out, _mm256_mul_ps(_mm256_loadu_ps(in), wav));
      }

      for (; q1 < end; q1++) {
        const float tmp = std::min(maximum, step * static_cast<float>(q1));
        output[s1 - q1].real(input[s1 - q1].real() *
                             mother[static_cast<int>(tmp)] * imaginary_sign);
        output[s1 - q1].imag(input[s1 - q1].imag() *
                             mother[static_cast<int>(tmp)]);
      }
    }
  }
#else
#ifndef SINGLE_THREAD
#pragma omp parallel for
#endif
  for (int i = 0; i < athreads; i++) {
    int start = (endpoint * i) / athreads;
    int end = (endpoint * (i + 1)) / athreads;

    for (int q1 = start; q1 < end; q1++) {
      float q = (float)q1;
      float tmp = std::min(maximum, step * q);

      output[q1].real(input[q1].real() * mother[(int)tmp]);
      output[q1].imag(input[q1].imag() * mother[(int)tmp] *
                      (1 - 2 * imaginary));
    }

    if (doublesided) {
      for (int q1 = start; q1 < end; q1++) {
        float q = (float)q1;
        float tmp = std::min(maximum, step * q);

        output[s1 - q1].real(input[s1 - q1].real() * mother[(int)tmp] *
                             (1 - 2 * imaginary));
        output[s1 - q1].imag(input[s1 - q1].imag() * mother[(int)tmp]);
      }
    }
  }

#endif
}

void fcwt::API::create_FFT_optimization_plan(const int maxsize,
                                             const int flags) const {
  (void)flags;
  const int nt = find2power(maxsize);
  if (nt <= 10) {
    std::cerr
        << "Maxsize is too small (<=1024)... please use a larger number\n";
    return;
  }

  for (int i = 11; i <= nt; i++) {
    int n = 1 << i;

    std::vector<float> dat(n);
    std::vector<std::complex<float>> O1(n);
    std::vector<std::complex<float>> out(n);

    std::cout << "Preparing KFR forward FFT plan with N:" << n << "\n";
    execute_forward_real(dat.data(), O1.data(), n);

    std::cout << "Preparing KFR backward FFT plan with N:" << n << "\n";
    FFTPlan p_back(n);
    p_back.execute_inverse(out.data(), O1.data());

    std::cout << "KFR FFT plans for N: " << n
              << " have been prepared. KFR chooses algorithms internally and "
                 "does not use FFTW wisdom files."
              << "\n";
  }
}

void fcwt::API::create_FFT_optimization_plan(int maxsize,
                                             const std::string &flags) {
  int flag = 0;

  if (flags == "FFTW_MEASURE") {
    flag = 0;
  } else if (flags == "FFTW_PATIENT") {
    flag = 0;
  } else if (flags == "FFTW_EXHAUSTIVE") {
    flag = 0;
  } else if (flags == "FFTW_ESTIMATE") {
    flag = 0;
  } else {
    std::cerr << "Unknown flag: " << flags << std::endl;
    return;
  }
  create_FFT_optimization_plan(maxsize, flag);
}

void fcwt::API::load_FFT_optimization_plan(const int input_size) const {
  const int nt = find2power(input_size);
  const int newsize = 1 << nt;

  if (use_optimalization_schemes) {
    if (newsize <= 1024) {
      std::cout << "Inputsize is too small (N <= 1024) to use optimization."
                << std::endl;
      return;
    }

    std::cout << "KFR does not use FFTW wisdom files; using KFR's internal "
                 "planning for N="
              << newsize << "." << std::endl;
  }
}

// Convolve in time domain using a single wavelet
void fcwt::API::convolve(const FFTPlan &p, std::complex<float> *Ihat,
                         std::complex<float> *O1, std::complex<float> *out,
                         Wavelet *wav, int size, int newsize, float scale,
                         bool lastscale) {

  if (lastscale) {
    std::vector<std::complex<float>> lastscalemem(newsize);

    fftbased(p, Ihat, O1, reinterpret_cast<float *>(lastscalemem.data()),
             wav->mother.data(), newsize, scale, wav->imag_frequency,
             wav->doublesided);
    if (use_normalization)
      fft_normalize(lastscalemem.data(), newsize);
    memcpy(out, lastscalemem.data(), sizeof(std::complex<float>) * size);
  } else {
    if (!out) {
      std::cout << "OUT NOT A POINTER" << std::endl;
    }
    fftbased(p, Ihat, O1, (float *)out, wav->mother.data(), newsize, scale,
             wav->imag_frequency, wav->doublesided);
    if (use_normalization)
      fft_normalize(out, newsize);
  }
}

void fcwt::API::fftbased(const FFTPlan &p, std::complex<float> *Ihat,
                         std::complex<float> *O1, float *out, float *mother,
                         int size, float scale, bool imaginary,
                         bool doublesided) {

  void *pt = out;

  // Perform daughter wavelet generation and multiplication with the Fourier
  // transformed input signal
  daughter_wavelet_multiplication(Ihat, O1, mother, scale, size, imaginary,
                                  doublesided);

  std::size_t space = 16;
  std::align(16, sizeof(std::complex<float>), pt, space);

  p.execute_inverse(reinterpret_cast<std::complex<float> *>(pt), O1);
}

void fcwt::API::fft_normalize(std::complex<float> *out, int size) {

  int nbatch = threads;
  int batchsize = (int)ceil((float)size / ((float)threads));

  // #pragma omp parallel for
  for (int i = 0; i < nbatch; i++) {
    int start = batchsize * i;
    int end = std::min(size, batchsize * (i + 1));

    for (int i8 = start; i8 < end; i8++) {
      out[i8] = out[i8] / (float)size;
    }
  }
}

void fcwt::API::cwt(float *pinput, int psize, std::complex<float> *poutput,
                    Scales *scales, bool complexinput) {
  const std::lock_guard lock(cwt_mutex);

  const int size = psize;

  // Find nearest power of 2
  const int nt = find2power(size);
  const int newsize = 1 << nt;

  std::vector<std::complex<float>> Ihat(newsize);
  std::vector<std::complex<float>> O1(newsize);

  // //Load optimization schemes if necessary
  load_FFT_optimization_plan(size);

  // //Perform forward FFT on input signal
  std::vector<float> real_input;
  std::vector<std::complex<float>> complex_input;
  if (complexinput) {
    complex_input.resize(newsize);
    memcpy(complex_input.data(), pinput, sizeof(std::complex<float>) * size);
    execute_forward_complex(complex_input.data(), Ihat.data(), newsize);
  } else {
    real_input.resize(newsize);
    memcpy(real_input.data(), pinput, sizeof(float) * size);
    execute_forward_real(real_input.data(), Ihat.data(), newsize);
  }

  const FFTPlan pinv(newsize);

  {
    const std::lock_guard wavelet_lock(wavelet->transform_mutex);

    // Generate mother wavelet function
    wavelet->generate(newsize);

    for (int i = 1; i < (newsize >> 1); i++) {
      Ihat[newsize - i] = std::conj(Ihat[i]);
    }

    std::complex<float> *out = poutput;

    for (int i = 0; i < scales->nscales; i++) {
      // FFT-base convolution in the frequency domain
      convolve(pinv, Ihat.data(), O1.data(), out, wavelet, size, newsize,
               scales->scales[i], i == (scales->nscales - 1));
      out = out + size;
    }
  }
}

void fcwt::API::cwt(float *pinput, int psize, std::complex<float> *poutput,
                    Scales *scales) {
  cwt(pinput, psize, poutput, scales, false);
}

void fcwt::API::cwt(std::complex<float> *pinput, int psize,
                    std::complex<float> *poutput, Scales *scales) {
  cwt((float *)pinput, psize, poutput, scales, true);
}

void fcwt::API::cwt(float *pinput, int psize, Scales *scales,
                    std::complex<float> *poutput, int pn1, int pn2) {
  assert((psize * scales->nscales) == (pn1 * pn2));
  cwt(pinput, psize, poutput, scales);
}

void fcwt::API::cwt(std::complex<float> *pinput, int psize, Scales *scales,
                    std::complex<float> *poutput, int pn1, int pn2) {
  assert((psize * scales->nscales) == (pn1 * pn2));
  cwt(pinput, psize, poutput, scales);
}
