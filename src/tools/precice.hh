#ifndef TOOLS_PRECICE_H
#define TOOLS_PRECICE_H

#include "precice/SolverInterface.hpp"
using namespace precice;
using namespace precice::constants;

#include "blocks/SWE_WavePropagationBlock.hh"

#include "tools/help.hh"

struct PreciceData{
  int  snd_heightId;
  int  snd_huId;
  int  snd_hvId;

  double* snd_height_db;
  double* snd_hu_db;
  double* snd_hv_db;

  int  recv_heightId;
  int  recv_huId;
  int  recv_hvId;

  double* recv_height_db;
  double* recv_hu_db;
  double* recv_hv_db;

  int * vertexIDs;
};

void recv_preCICE(SolverInterface &interface, SWE_Block &l_wavePropgationBlock, SWE_Block1D* ghoshtBlock,
                  SWE_Block1D* newBlock, PreciceData *data, int columNr, int size);

void snd_preCICE(SolverInterface &interface, SWE_Block &l_wavePropgationBlock, PreciceData *data, int columNr, int size);

#endif // TOOLS_PRECICE_H
