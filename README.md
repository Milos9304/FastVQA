# FastVQA
Simulation framework of variational quantum algorithms focused on performance, portability and distributivity on parallel architectures. Although any quantum backends can be used, only [QuEST](https://github.com/QuEST-Kit/QuEST), a highly portable statevector simulator, has been included so far. The implementation of variational quantum algorithms found in src/vqas was inspired by [xacc](https://github.com/eclipse/xacc) library.

## Features
- implementation of a Variational Quantum Eigensolver (VQE)
	- [VQE CVaR version](https://doi.org/10.22331/q-2020-04-20-256)
	- [VQE Ascending CVaR version](https://arxiv.org/abs/2105.11766)
- implementation of [Adiabatic Quantum Computation via Parameterized Quantum Circuits](https://arxiv.org/abs/2206.04373)
- implementation of QAOA (testing)
- parallelized computation via MPI library (testing)

## Installation
```
git submodule update --init --recursive
mkdir build && cd build
cmake .. && make
make install
```

### Custom installation path
```
cmake .. -DCMAKE_INSTALL_PREFIX=destination_path
```

# Usage
See [LattiQ](https://github.com/Milos9304/LattiQ), the experimental framework for lattice problems on NISQ devices that is running on top of FastVQA and [AqcPqC] for an experimental framework for Adiabatic Quantum Computation via Parameterized Quantum Circuits built on top of FastVQA. More detailed docs and api are on plan soon.

