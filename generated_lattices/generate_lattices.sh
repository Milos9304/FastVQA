#!/bin/bash

seed=77

latticegen="${HOME}/.local/bin/latticegen -randseed ${seed}"

$latticegen q 10 1 20 b > "q_10_1_20_b"
$latticegen q 130 1 20 b > "q_130_1_20_b"

