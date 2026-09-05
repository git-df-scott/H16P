#!/usr/bin/env bash
set -euo pipefail
gcc -O3 -shared -fPIC -fopenmp audit/fable_engine/retmap.c -o audit/fable_engine/libretmap.so -lm
gcc -O3 -shared -fPIC -fopenmp audit/fable_engine/retmap_log.c -o audit/fable_engine/libretmap_log.so -lm
g++ -O2 -std=c++17 -fext-numeric-literals fastra_d1_2026_09_05/half_beta_quad.cpp -o /tmp/d1_half -lquadmath
g++ -O2 -std=c++17 -fext-numeric-literals fastra_d1_2026_09_05/full_beta_quad.cpp -o /tmp/d1_full -lquadmath
