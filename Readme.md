# EHT- C
[![Build Status](https://travis-ci.com/AlessandroBudroni/EHT-C.svg?branch=main)](https://travis-ci.com/AlessandroBudroni/EHT-C)  [![Coverage Status](https://coveralls.io/repos/github/AlessandroBudroni/EHT-C/badge.svg?branch=main)](https://coveralls.io/github/AlessandroBudroni/EHT-C?branch=main) [![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/) 

![GitHub issues](https://img.shields.io/github/issues/AlessandroBudroni/EHT-C) ![GitHub pull requests](https://img.shields.io/github/issues-pr/AlessandroBudroni/EHT-C) ![GitHub repo size](https://img.shields.io/github/repo-size/AlessandroBudroni/EHT-C)

**Errors Hidden Trapdoor** is a public-key cryptosystem introduced with the manuscript [will be soon available](..) based on the well known LWE problem.

This repository contains a working proof-of-concept of EHT written in pure C with no external dependencies. 

### How to build
The code has been tested only on machines with Intel Core processors running Ubuntu/Linux Mint and GCC 9.3.0 as compiler. To build the repository, run the following

```
mkdir target
cd target
cmake ..
make
```
Run correctness tests
```
make test
```
To run the benchmark test for EHT-high-A, run the following

```
./test/EHT_high_A_bench_cycles

```
By default the memory is allocated both in heap and stack. One can choose the option of using only stack memory and achieve the best performances with the following
```
cmake -D FULL_STACK=ON ..
```
With the above, it may be required for some machines to increase the stack size available with `ulimit -s 131070` to run EHT-medium-\* and EHT-high-\*.

**We do not claim this code to be secure against side channel attacks and therefore it should not be used in production.**

### Licence
The code is released under the [GNU General Public License version 3](https://opensource.org/licenses/GPL-3.0).

The function `sample_error(...)` in `src/encrypt.c` was taken and modified from a function in [https://github.com/Microsoft/PQCrypto-LWEKE](https://github.com/Microsoft/PQCrypto-LWEKE)  that is relased under the [MIT Licence](https://opensource.org/licenses/MIT) .

### Authors:

- Alessandro Budroni
