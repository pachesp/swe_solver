#!/bin/bash
#script sourced by ../runSWE
cd SWE_solver
./compile
cd ..
./SWE_solver/build/SWE_gnu_debug_none_fwave_1 -x $size -y $size -t $type -o $output
