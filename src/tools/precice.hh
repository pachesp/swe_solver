#ifndef TOOLS_PRECICE_H
#define TOOLS_PRECICE_H

#include "precice/SolverInterface.hpp"
using namespace precice;
using namespace precice::constants;

#include "blocks/SWE_WavePropagationBlock.hh"
#include "tools/help.hh"

enum type{
    twoDtwoDsup,
    twoDtwoDsub,
    twoDthreeDdsup,
    twoDthreeDdsub,
    threeDtwoDdsup,
    threeDtwoDdsub,
    threeDthreeDsup,
    threeDthreeDsub
};

struct PreciceData{
    int  snd_heightId;
    int  snd_huId;
    int  snd_hvId;

    double* snd_height_db;
    double* snd_hu_db;
    double* snd_hv_db;

    int recv_heightId;
    int recv_huId;
    int recv_hvId;

    double* recv_height_db;
    double* recv_hu_db;
    double* recv_hv_db;

    int* vertexIDs;




  int  alphaId;
  int  ghId;
  int  velocityId;
  int  velocityGradientId;
  double* grid;
  int l_nY;
  double l_dY;
  size_t simType;

  int heightGradId;
  int huGradId;
  int hvGradId;
  double* heightGrad_db;
  double* huGrad_db;
  double* hvGrad_db;

  int heightId;
  int huId;
  int hvId;
  double* heightS1_db;
  double* huS1_db;
  double* hvS1_db;


  Float2D* CP_height_f2d;
  Float2D* CP_hu_f2d;
  Float2D* CP_hv_f2d;
};

void write2SWE_SWE_Right_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  int size, int columNr);

void write2SWE_SWE_Left_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                int size, int columNr);

void readFromSWE_SWE_Left_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, SWE_Block1D* ghoshtBlock,
                  PreciceData *data, int size, int columNr = 0);

void readFromSWE_SWE_Right_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock,
                  PreciceData *data, int size, int columNr = 0);

void writeCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float time, float &time_CP, int size, int columNr);

void restoreCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float &time, float time_CP, int size, int columNr);

//-----------------------To interfoam-----------------------------------------------------

void write2Interfoam_2D3DSupercritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                int columNr);

void readFromInterfoam_3D2DSupercritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
    SWE_Block1D* ghoshtBlock, int columNr = 0);

void write2Interfoam_2D3DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                int columNr, std::vector<double> alpha);

std::vector<double> readFromInterfoam_2D3DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                SWE_Block1D* ghoshtBlock, int columNr = 0);


void write2Interfoam_3D2DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  int columNr);

void readFromInterfoam_3D2DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
    SWE_Block1D* ghoshtBlock, int columNr=0);

// void storeData_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
//                      int columNr, int size);
//
// void writeData_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
//                 int columNr, int size);
#endif // TOOLS_PRECICE_H
