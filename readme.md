# fCWT

Highly optimized C++ implementation of the Continuous Wavelet Transform.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

fCWT uses KFR for FFT/DFT operations. CMake first tries `find_package(KFR CONFIG)`;
if KFR is not installed, it downloads KFR tag `7.0.1` by default:

```sh
cmake -S . -B build -DFCWT_FETCH_KFR=ON
```

To use an installed KFR package instead, set `CMAKE_PREFIX_PATH` and disable fetching:

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/kfr/install -DFCWT_FETCH_KFR=OFF
```

KFR is GPL/commercial licensed; check that license choice fits your distribution.

## Format

```sh
just format
```
