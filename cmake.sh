#!/bin/bash

CC=gcc
CXX=g++
FC=gfortran

BUILD_TYPE=Release
XACC_DIR="~/.xacc"

cmake \
  -DCMAKE_CXX_COMPILER=${CXX} \
  -DXACC_DIR=${XACC_DIR} \
..
