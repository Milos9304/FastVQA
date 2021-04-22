#!/bin/bash

CC=icc #gcc
CXX=icpc #g++
#FC=gfortran

BUILD_TYPE=Release
XACC_DIR="~/.xacc"
GMP_H_DIR="/home/nx05/nx05/milos_pr/.local/include"
GMP_L_DIR="/home/nx05/nx05/milos_pr/.local/lib"

cmake \
  -DCMAKE_CXX_COMPILER=${CXX} \
  -DXACC_DIR=${XACC_DIR} \
  -DGMP_H_DIR=${GMP_H_DIR} \
  -DGMP_L_DIR=${GMP_L_DIR} \
..
