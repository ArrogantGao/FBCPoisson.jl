# FBCPoisson.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ArrogantGao.github.io/FBCPoisson.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ArrogantGao.github.io/FBCPoisson.jl/dev/)
[![Build Status](https://github.com/ArrogantGao/FBCPoisson.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ArrogantGao/FBCPoisson.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ArrogantGao/FBCPoisson.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ArrogantGao/FBCPoisson.jl)

This package is used to solve the Poisson equation under free boundary conditions in 3D with smoonth continuous charge distribution.

Key techniques used:
1. Truncated Green's function method to handle free boundary conditions.
2. Non-uniform FFT to efficiently compute the convolution with the Green's function.