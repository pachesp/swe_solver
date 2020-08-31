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
    threeDtwoDdsub
};

struct PreciceData{
  int  alphaId;
  int  ghId;
  int  velocityId;
  int  velocityGradientId;
  int l_nY;
  double l_dY;
  size_t simType;

  int* vertexIDs;
  double* grid;

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

class PreciceExchange {
  protected:
    SolverInterface& interface;
    SWE_Block& wavePropagationBlock;
    PreciceData* data;
    int columNr;
    SWE_Block1D* ghoshtBlock;

    std::vector<double> alphaVector;

  public:
    PreciceExchange(SolverInterface& interface,
        SWE_Block& wavePropagationBlock, PreciceData* data, int columNr, SWE_Block1D* ghoshtBlock = nullptr);
    virtual void write() = 0;
    virtual void read() = 0;

    void setAlphaVector(std::vector<double> data){
      alphaVector = data;
  }
};

class SWE_SWE_SuperCritical_Left_Exchange : public PreciceExchange{
  public:
    SWE_SWE_SuperCritical_Left_Exchange(SolverInterface& interface,  SWE_Block& wavePropagationBlock, PreciceData* data,
      int columNr, SWE_Block1D* ghoshtBlock);
      virtual void write();
      virtual void read();
};

class SWE_SWE_SuperCritical_Right_Exchange : public PreciceExchange{
  public:
    SWE_SWE_SuperCritical_Right_Exchange(SolverInterface& interface,  SWE_Block& wavePropagationBlock, PreciceData* data,
      int columNr, SWE_Block1D* ghoshtBlock);
      virtual void write();
      virtual void read();
};

class SWE_OF_SuperCritical_Exchange : public PreciceExchange{
  public:
    SWE_OF_SuperCritical_Exchange(SolverInterface& interface,  SWE_Block& wavePropagationBlock, PreciceData* data,
      int columNr, SWE_Block1D* ghoshtBlock = nullptr);
      virtual void write();
      virtual void read();
};

class SWE_OF_SubCritical_Exchange : public PreciceExchange{
  public:
    SWE_OF_SubCritical_Exchange(SolverInterface& interface,  SWE_Block& wavePropagationBlock, PreciceData* data,
      int columNr, SWE_Block1D* ghoshtBlock);
      virtual void write();
      virtual void read();
};

class OF_SWE_SuperCritical_Exchange : public PreciceExchange{
  public:
    OF_SWE_SuperCritical_Exchange(SolverInterface& interface,  SWE_Block& wavePropagationBlock, PreciceData* data,
      int columNr, SWE_Block1D* ghoshtBlock);
      virtual void write();
      virtual void read();
};

class OF_SWE_SubCritical_Exchange : public PreciceExchange{
  public:
    OF_SWE_SubCritical_Exchange(SolverInterface& interface,  SWE_Block& wavePropagationBlock, PreciceData* data,
      int columNr, SWE_Block1D* ghoshtBlock);
      virtual void write();
      virtual void read();
};


void writeCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float time, float &time_CP, int size, int columNr);

void restoreCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float &time, float time_CP, int size, int columNr);

#endif // TOOLS_PRECICE_H
