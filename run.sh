#!/bin/bash
#script sourced by ../<case>/runSWE<x>

cd swe_solver
. ./compile
#TODO add a way of making sure that the compiling was successful, and only then continue
cd ..
./swe_solver/build/SWE_gnu_debug_none_fwave_$solverNum -x $sizeX -y $sizeY -t $type -o $output
