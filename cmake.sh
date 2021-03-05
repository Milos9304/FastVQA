#!/bin/bash

CC=gcc
CXX=g++
FC=gfortran

BUILD_TYPE=Release
#INSTALL_DIR=$HOME/Quantum/XACC/installdir/xacc
BUILD_TESTS=TRUE
BUILD_EXAMPLES=TRUE
XACC_DIR="~/.xacc"

cmake \
  -DCMAKE_C_COMPILER=${CC} \
  -DCMAKE_CXX_COMPILER=${CXX} \
  -DCMAKE_Fortran_COMPILER=${FC} \
  -DXACC_DIR=${XACC_DIR} \
..
