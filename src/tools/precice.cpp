#include "tools/precice.hh"
#include <iomanip>      // std::setprecision

PreciceExchange::PreciceExchange(SolverInterface& interface,
  SWE_Block& wavePropagationBlock,
  PreciceData* data,
  int columNr,
  SWE_Block1D* ghoshtBlock):
  interface(interface),
  wavePropagationBlock(wavePropagationBlock),
  data(data),
  columNr(columNr),
  ghoshtBlock(ghoshtBlock)
  {  };

//case 0
SWE_SWE_SuperCritical_Left_Exchange::SWE_SWE_SuperCritical_Left_Exchange(
  SolverInterface& interface,
  SWE_Block& wavePropagationBlock, PreciceData* data,
  int columNr, SWE_Block1D* ghoshtBlock) :
  PreciceExchange(interface, wavePropagationBlock, data, columNr, ghoshtBlock)
  {};

void SWE_SWE_SuperCritical_Left_Exchange::write(){
  std::cout << "executing 2d to 2d supercritical or subcritical Left write" << '\n';

  int nY = data->l_nY;
  double h[nY + 2]{};
  double hu[nY + 2]{};
  double hv[nY + 2]{};

  for(int i = 0; i < nY + 2 ; i++){ //Data stored column-wise
   h[i] = (double)wavePropagationBlock.getWaterHeight()[columNr][i];
   hu[i] = (double)wavePropagationBlock.getDischarge_hu()[columNr][i];
   hv[i] = (double)wavePropagationBlock.getDischarge_hv()[columNr][i];
  }

  interface.writeBlockScalarData(data->heightId, nY + 2, data->vertexIDs, h);
  interface.writeBlockScalarData(data->huId, nY + 2, data->vertexIDs, hu);
  interface.writeBlockScalarData(data->hvId, nY + 2, data->vertexIDs, hv);
}

void SWE_SWE_SuperCritical_Left_Exchange::read(){
  std::cout << "executing 2d to 2d supercritical or subcritical Left read" << '\n';

     int nY = data->l_nY;
     double hGrad[nY + 2]{};
     double huGrad[nY + 2]{};
     double hvGrad[nY + 2]{};

     interface.readBlockScalarData(data->heightGradId, nY + 2, data->vertexIDs, hGrad);
     interface.readBlockScalarData(data->huGradId, nY + 2, data->vertexIDs, huGrad);
     interface.readBlockScalarData(data->hvGradId, nY + 2, data->vertexIDs, hvGrad);

     for (int i = 0; i < nY + 2; i++) {
       wavePropagationBlock.hGrad[i] = (float)hGrad[i];
       wavePropagationBlock.huGrad[i] = (float)huGrad[i];
       wavePropagationBlock.hvGrad[i] = (float)hvGrad[i];
    }
}


//case 0
SWE_SWE_SuperCritical_Right_Exchange::SWE_SWE_SuperCritical_Right_Exchange(
  SolverInterface& interface,
  SWE_Block& wavePropagationBlock, PreciceData* data,
  int columNr, SWE_Block1D* ghoshtBlock) :
  PreciceExchange(interface, wavePropagationBlock, data, columNr, ghoshtBlock){}

void SWE_SWE_SuperCritical_Right_Exchange::write(){
  std::cout << "executing 2d to 2d supercritical Right write" << '\n';

     int nY = data->l_nY;
     double hGrad[nY + 2]{};
     double huGrad[nY + 2]{};
     double hvGrad[nY + 2]{};

     for(int i = 0; i < nY + 2; i++){ //Data stored column-wise
       hGrad[i] = (double)wavePropagationBlock.hGrad[i];
       huGrad[i] = (double)wavePropagationBlock.huGrad[i];
       hvGrad[i] = (double)wavePropagationBlock.hvGrad[i];
     }

     interface.writeBlockScalarData(data->heightGradId, nY + 2, data->vertexIDs, hGrad);
     interface.writeBlockScalarData(data->huGradId, nY + 2, data->vertexIDs, huGrad);
     interface.writeBlockScalarData(data->hvGradId, nY + 2, data->vertexIDs, hvGrad);
}

void SWE_SWE_SuperCritical_Right_Exchange::read(){
  std::cout << "executing 2d to 2d supercritical Right read" << '\n';

     int nY = data->l_nY;
     double h[nY + 2]{};
     double hu[nY + 2]{};
     double hv[nY + 2]{};

     float h_f[nY + 2]{};
     float hu_f[nY + 2]{};
     float hv_f[nY + 2]{};

     interface.readBlockScalarData(data->heightId, nY + 2, data->vertexIDs, h);
     interface.readBlockScalarData(data->huId, nY + 2, data->vertexIDs, hu);
     interface.readBlockScalarData(data->hvId, nY + 2, data->vertexIDs, hv);

     for(int i = 0; i< nY + 2; i++){
       h_f[i] = (float)h[i];
       hu_f[i] = (float)hu[i];
       hv_f[i] = (float)hv[i];
     }

     // Data sent by left neighbour from precice
     SWE_Block1D newBlock{h_f, hu_f, hv_f, 1};
     ghoshtBlock->copyFrom(&newBlock, nY + 2);
}


//case 2
SWE_OF_SuperCritical_Exchange::SWE_OF_SuperCritical_Exchange(
 SolverInterface& interface,
 SWE_Block& wavePropagationBlock, PreciceData* data,
 int columNr, SWE_Block1D* ghoshtBlock) :
 PreciceExchange(interface, wavePropagationBlock, data, columNr){}

void SWE_OF_SuperCritical_Exchange::write(){
  std::cout << "executing 2d to 3d supercritical write" << '\n';
  int nY = data->l_nY; //resolution
  double dY = data->l_dY; //cell size in y-direction
  int dim = 3;

  // Writing Alpha to OF - Begin
  double alpha[nY * nY] {};// = new double[nY * nY]; //scalar alpha holder for sending it to IF
  for (int i = 0; i < nY; i++) {
      double height = (double)wavePropagationBlock.getWaterHeight()[columNr][i];
      for (int j = 0; j < nY; j++) {
          double yCoord = data->grid[(i + j*nY) * dim + 1]; // y-coord of the IF mesh
          //Algorithm 1, Mintegen page 84
          if(height <= yCoord - 0.5 * dY ){
              alpha[i + j*nY] = 0.0;
          } else if(height >= yCoord + 0.5 * dY){
               alpha[i + j*nY] = 1.0;
          } else{ // Interpolate
              alpha[i + j*nY] = 0.5 + (height - yCoord) / dY;
          }
      }
  }
  // Writing Alpha to OF - End


  // Writing Velocity to OF - Begin

      double U[dim * nY * nY] {}; //3D velocity holder for sending it to IF

      double u, v, h;
      for(int i = 1; i <= nY ; i++){
          h = (double)(wavePropagationBlock.getWaterHeight()[columNr][i]);
          u = (double)(wavePropagationBlock.getDischarge_hu()[columNr][i]) / h;
          v = (double)(wavePropagationBlock.getDischarge_hv()[columNr][i]) / h;
          for (int j = 0; j < nY; j++) {
              U[((i-1) + j*nY) * dim + 0] = u * alpha[(i-1) + j*nY];
              U[((i-1) + j*nY) * dim + 1] = v * alpha[(i-1) + j*nY];
              U[((i-1) + j*nY) * dim + 2] = 0.0;
          }
      }

  //     // Mintgens algorithm for treating velocity with shear stress
  //     int dim2 = 2;
  //     double rho_water = 1000.0; // water density
  //     double tau[dim2 * nY * nY]{0.0}; // wall shear stress
  //     double ustar[dim2 * nY * nY]{0.0}; // modified velocity: sqrt(tau / rho)
  //     double mu = 1e-06 * rho_water; // dynamic viscosity * density (water)
  //     double du[dim * nY * nY]{0.0}; // velocity gradient from OF
  //     interface.readBlockVectorData(data->velocityGradientId, nY * nY, data->vertexIDs, du);
  //     for (int i = 0; i < nY; i++) {
  //         for (int j = 0; j < nY; j++) {
  //             tau[(i + j*nY) * dim + 0] =abs(mu * du[(i + j*nY) * dim + 0] / dY);
  //             tau[(i + j*nY) * dim + 1] =abs(mu * du[(i + j*nY) * dim + 1] / dY);
  //         }
  //     }
  //
  //     int count = 0;
  //     for (int i = 0; i < nY; i++) {
  //         for (int j = 0; j < nY; j++) {
  //             double yCoord = data->grid[(i + j*nY) * dim + 1];
  //             ustar[(i + j*nY) * dim2 + 0] = std::sqrt(tau[(i + j*nY) * dim2 + 0] / rho_water);
  //             ustar[(i + j*nY) * dim2 + 1] = std::sqrt(tau[(i + j*nY) * dim2 + 1] / rho_water);
  //             std::cout << fixed << setprecision(4) << ustar[(i + j*nY) * dim2 + 0] << ", ";
  //             // std::cout << ustar[(i + j*nY) * dim2 + 2] << ' ';
  //             U[count++] = u + ustar[(i + j*nY) * dim + 0] / 0.41 * (1+std::log(yCoord/h));
  //             U[count++] = v + ustar[(i + j*nY) * dim + 1] / 0.41 * (1+std::log(yCoord/h));
  //             U[count++] = 0.0;
  //         }
  //         std::cout  << '\n';
  //     }
  //
  //     double q3d[dim2]{};
  //     double q2d[dim2]{};
  //     double beta[dim2]{};
  //     for (int i = 0; i < nY; i++) {
  //         q2d[0] = (double)(wavePropagationBlock.getDischarge_hu()[columNr][i]);
  //         q2d[1] = (double)(wavePropagationBlock.getDischarge_hv()[columNr][i]);
  //         for (int j = 0; j < nY; j++) {
  //             double yCoord = data->grid[(i + j*nY) * dim + 1];      // y-coord of the IF mesh
  //             q3d[0] += alpha[i + j*nY] * U[(i + j*nY) * dim + 0] * yCoord;
  //             q3d[1] += alpha[i + j*nY] * U[(i + j*nY) * dim + 1] * yCoord;
  //         }
  //     }
  //
  //     std::cout << "THE BETA" << '\n';
  //     for (int n = 0; n < dim2; n++) {
  //       if (q3d[n]< 0.0001) {
  //         q3d[n] = 1.0;
  //       }
  //       beta[n] = abs(q2d[n] / q3d[n]);
  //       std::cout << beta[n] << '\n';
  //
  //       for (int i = 0; i < nY; i++) {
  //         for (int j = 0; j < nY; j++) {
  //           q3d[n] = 0.0;
  //           if (beta[n] > 1) {
  //             std::cout << "INTO THE WILE" << '\n';
  //             while (beta[n] > 1 ) {
  //               double yCoord = data->grid[(i + j*nY) * dim + 1];
  //               U[(i + j*nY) * dim + n] = U[(i + j*nY) * dim + n] * beta[n];
  //               std::cout<< fixed << setprecision(3) << U[(i + j*nY) * dim + n]  << ", ";
  //               q3d[n] += alpha[i + j*nY] * U[(i + j*nY) * dim + n] * yCoord;
  //               std::cout << "q is: " << '\n';
  //               std::cout << q3d[n] << ',';
  //               if (q3d[n]< 0.0001) {
  //                 q3d[n] = 1.0;
  //               }
  //               beta[n] = abs(q2d[n] / q3d[n]);
  //             }
  //           }
  //         // else if (beta[n] < 1) { //TODO missing to take care of this case, so far this breaks it all
  //          //     std::cout << "beta[" << n << "]" << beta[n] << '\n';
  //          //     while (beta[n] < 1 ) {
  //          //         for (int i = 0; i < nY; i++) {
  //          //             q3d[n] = 0.0;
  //          //             for (int j = 0; j < nY; j++) {
  //          //                 double yCoord = data->grid[(i + j*nY) * dim + 1];
  //          //                 U[(i + j*nY) * dim + n] = U[(i + j*nY) * dim + n] / beta[n];
  //          //                 q3d[n] += alpha[i + j*nY] * U[(i + j*nY) * dim + n] * yCoord;
  //          //             }
  //          //             beta[n] = q2d[n] / q3d[n];
  //          //         }
  //          //     }
  //          // }
  //         }
  //       }
  //     }
  //
  //     double umax[dim2]{0.0};
  //     for (int n = 0; n < dim2; n++) {
  //         for (int i = 0; i < nY; i++) {
  //             for (int j = 0; j < nY; j++) {
  //                 if (umax[n] < alpha[i*nY + j] * U[i*nY + j]) {
  //                     umax[n] = alpha[i*nY + j] * U[i*nY + j];
  //                 }
  //             }
  //             for (int j = 0; j < nY; j++) {
  //                 if (alpha[i*nY + j] < 0.001) {
  //                     U[i*nY + j] = umax[n];
  //                 }
  //             }
  //         }
  //     }
  // // Writing Velocity to OF - End
  // std::cout << "VELOCITY" << '\n';
  // for (int i = 0; i < nY; i++) {
  //   for (int j = 0; j < nY; j++) {
  //     std::cout <<fixed << setprecision(2) << U[(i*nY + j)*dim] <<", ";
  //   }
  //   std::cout  << '\n';
  // }

  // Exchange data
  interface.writeBlockScalarData(data->alphaId, nY * nY, data->vertexIDs, alpha);
  interface.writeBlockVectorData(data->velocityId, nY * nY, data->vertexIDs, U);
}

void SWE_OF_SuperCritical_Exchange::read(){
  std::cout << "executing 2d to 3d supercritical read. (empty)" << '\n';
 // intentionally empty
}


//case 3
SWE_OF_SubCritical_Exchange::SWE_OF_SubCritical_Exchange(
 SolverInterface& interface,
 SWE_Block& wavePropagationBlock, PreciceData* data,
 int columNr, SWE_Block1D* ghoshtBlock) :
 PreciceExchange(interface, wavePropagationBlock, data, columNr, ghoshtBlock){}

void SWE_OF_SubCritical_Exchange::write(){
    std::cout << "executing 2d to 3d subcritical write" << '\n';

    // Writing Velocity to OF - Begin
    int nY = data->l_nY; //resolution
    int dim = 3;

    // Algorithm 1
    double U[dim * nY * nY] {}; //3D velocity holder for sending it to IF
    double u, v, h;
    for(int i = 1; i <= nY ; i++){
        h = (double)(wavePropagationBlock.getWaterHeight()[columNr][i]);
        u = (double)(wavePropagationBlock.getDischarge_hu()[columNr][i]) / h;
        v = (double)(wavePropagationBlock.getDischarge_hv()[columNr][i]) / h;

        for (int j = 0; j < nY; j++) {
            double alphatmp = alphaVector[(i-1) + j*nY];
            U[((i-1) + j*nY) * dim + 0] = alphatmp * u;
            U[((i-1) + j*nY) * dim + 1] = alphatmp * v;
            U[((i-1) + j*nY) * dim + 2] = 0.0;
        }
    }

    interface.writeBlockVectorData(data->velocityId, nY * nY, data->vertexIDs, U);
  }

void SWE_OF_SubCritical_Exchange::read(){
    std::cout << "executing 2d to 3d subcritical read" << '\n';

    int nY = data->l_nY; //resolution
    double dY = data->l_dY; //cell size in y-direction

    // read alpha
    double alpha[nY * nY]{};
    interface.readBlockScalarData(data->alphaId, nY * nY, data->vertexIDs, alpha);

    float h[nY+2] {};
    for (int i = 0; i < nY; i++) {
        for (int j = 1; j <= nY; j++) {
            h[j] += (float)alpha[i*nY + (j-1)] * dY;
        }
    }

    std::vector<double> alphav;
    for (int i = 0; i < nY * nY; i++) {
        alphav.push_back(alpha[i]);
    }

    //gradient 0 for discharges hu and hv
    float hu[nY+2] {};
    float hv[nY+2] {};
    for (int i = 0; i < nY+2; i++) {
            hu[i] = wavePropagationBlock.getDischarge_hu()[columNr][i]; // TODO hardcode the column
            hv[i] = wavePropagationBlock.getDischarge_hv()[columNr][i];
    }

    SWE_Block1D newBlock{h, hu, hv, 1};
    ghoshtBlock->copyFrom(&newBlock, nY+2);

    setAlphaVector(alphav);
  }


//case 4
OF_SWE_SuperCritical_Exchange::OF_SWE_SuperCritical_Exchange(
 SolverInterface& interface,
 SWE_Block& wavePropagationBlock, PreciceData* data,
 int columNr, SWE_Block1D* ghoshtBlock) :
 PreciceExchange(interface, wavePropagationBlock, data, columNr, ghoshtBlock){}

void OF_SWE_SuperCritical_Exchange::write(){
  std::cout << "executing 3d to 2d supercritical write. (empty)" << '\n';
  //intentionally empty
  }

void OF_SWE_SuperCritical_Exchange::read(){
    std::cout << "executing 3d to 2d supercritical read" << '\n';

    int nY = data->l_nY; //resolution
    float dY = (float)data->l_dY; //cell size in y-direction
    int dim = 3; //3-dimension

    // read alpha
    double alpha[nY * nY]{};
    interface.readBlockScalarData(data->alphaId, nY * nY, data->vertexIDs, alpha);
    float h[nY+2] {0.0f};
    for (int i = 0; i < nY; i++) {
        for (int j = 1; j <= nY; j++) {
            h[j] += (float)alpha[i*nY + (j-1)] * dY;
        }
    }

    // read velocity
    double U[dim * nY * nY]{};
    interface.readBlockVectorData(data->velocityId, nY * nY, data->vertexIDs, U);
    float hu[nY+2] {};
    float hv[nY+2] {};
    for (int i = 0; i < nY; i++) {
        for (int j = 1; j <= nY; j++) {
                hu[j] += (float)U[(i*nY + (j-1)) * dim + 0] * (float)alpha[i*nY + (j-1)] * dY;
                hv[j] += (float)U[(i*nY + (j-1)) * dim + 1] * (float)alpha[i*nY + (j-1)] * dY;
        }
    }

    SWE_Block1D newBlock{h, hu, hv, 1};
    ghoshtBlock->copyFrom(&newBlock, nY+2);
  }


//case 5
OF_SWE_SubCritical_Exchange::OF_SWE_SubCritical_Exchange(
 SolverInterface& interface,
 SWE_Block& wavePropagationBlock, PreciceData* data,
 int columNr, SWE_Block1D* ghoshtBlock) :
 PreciceExchange(interface, wavePropagationBlock, data, columNr, ghoshtBlock){}

void OF_SWE_SubCritical_Exchange::write(){
    std::cout << "executing 3d to 2d subcritical write" << '\n';
    int nY = data->l_nY; //resolution
    double g = 9.81;
    double h;
    double gh[nY * nY]{};
    for (int i = 0; i < nY; i++) {
        h = wavePropagationBlock.getWaterHeight()[columNr][i];
        for (int j = 0; j < nY; j++) {
            gh[i + j*nY] = g * h;           // Set p_rgh  eq 5.16 mintgen
        }
    }

    interface.writeBlockScalarData(data->ghId, nY * nY, data->vertexIDs, gh);
  }

void OF_SWE_SubCritical_Exchange::read(){
    std::cout << "executing 3d to 2d subcritical read" << '\n';
    int nY = data->l_nY; //resolution
    double dY = data->l_dY; //cell size in y-direction
    int dim = 3;

    // read alpha
    double alpha[nY * nY]{};
    interface.readBlockScalarData(data->alphaId, nY * nY, data->vertexIDs, alpha);
    // read velocity
    double U[dim * nY * nY]{};
    interface.readBlockVectorData(data->velocityId, nY * nY, data->vertexIDs, U);
    float hu[nY+2] {};
    float hv[nY+2] {};
    for (int i = 0; i < nY; i++) {
        for (int j = 1; j <= nY; j++) {
                hu[j] += (float)U[(i*nY + (j-1)) * dim + 0] * (float)alpha[i*nY + (j-1)] * dY;
                hv[j] += (float)U[(i*nY + (j-1)) * dim + 1] * (float)alpha[i*nY + (j-1)] * dY;
        }
    }

    //gradient 0 for h
    float h[nY+2] {};
    for (int i = 0; i < nY+2; i++) {
            h[i] = wavePropagationBlock.getWaterHeight()[1][i];
    }

    SWE_Block1D newBlock{h, hu, hv, 1};
    ghoshtBlock->copyFrom(&newBlock, nY+2);
  }



void writeCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float time, float &time_CP, int size, int columNr)
{
    if(data->CP_height_f2d) delete data->CP_height_f2d;
   if(data->CP_hu_f2d) delete data->CP_hu_f2d;
   if(data->CP_hv_f2d) delete data->CP_hv_f2d;


   time_CP = time;

   data->CP_height_f2d = new Float2D(wavePropagationBlock.getWaterHeight(), false);
   data->CP_hu_f2d = new Float2D(wavePropagationBlock.getDischarge_hu(), false);
   data->CP_hv_f2d = new Float2D(wavePropagationBlock.getDischarge_hv(), false);
}

void restoreCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float &time, float time_CP, int size, int columNr)
{
    time = time_CP;

    for(int j = 0; j < size ; j++){
        for(int i = 0; i < size ; i++){

            wavePropagationBlock.getWaterHeight()[i][j]  =  (*(data->CP_height_f2d))[i][j];
            wavePropagationBlock.getDischarge_hu()[i][j] =  (*(data->CP_hu_f2d))[i][j];
            wavePropagationBlock.getDischarge_hv()[i][j] =  (*(data->CP_hv_f2d))[i][j];
        }
    }
}
