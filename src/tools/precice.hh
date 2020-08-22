#ifndef TOOLS_PRECICE_H
#define TOOLS_PRECICE_H

#include "precice/SolverInterface.hpp"
using namespace precice;
using namespace precice::constants;

#include "blocks/SWE_WavePropagationBlock.hh"
#include "tools/help.hh"

enum type{
    twoDthreeDdsup=2,
    twoDthreeDdsub,
    threeDtwoDdsup,
    threeDtwoDdsub
};

struct PreciceData{

  int  alphaId;
  int  prghId;
  int  velocityId;
  int  tempVelocity3dId;


  double* grid3D;
  int l_nY;
  double l_dY;
  int* vertexIDs;
  size_t simType;

  Float2D* CP_height_f2d;
  Float2D* CP_hu_f2d;
  Float2D* CP_hv_f2d;

};

void write_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  int size, int columNr);

void writeGradient_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                int size, int columNr);

void read_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, SWE_Block1D* ghoshtBlock,
                  PreciceData *data, int size, int columNr = 0);

void readGradient_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock,
                  PreciceData *data, int size, int columNr = 0);

void writeCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float time, float &time_CP, int size, int columNr);

void restoreCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float &time, float time_CP, int size, int columNr);

//-----------------------To interfoam-----------------------------------------------------

void write2Interfoam_2D3DSupercritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                int columNr, double* tempVelocity3d);

void readFromInterfoam_3D2DSupercritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
    SWE_Block1D* ghoshtBlock, int columNr = 0);

void write2Interfoam_2D3DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                int columNr, double* tempVelocity3d);

void readFromInterfoam_2D3DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                SWE_Block1D* ghoshtBlock, int columNr = 0);


void write2Interfoam_3D2DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  int columNr, double* tempVelocity3d);

void readFromInterfoam_3D2DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
    SWE_Block1D* ghoshtBlock, int columNr=0);
// void storeData_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
//                      int columNr, int size);
//
// void writeData_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
//                 int columNr, int size);
#endif // TOOLS_PRECICE_H
