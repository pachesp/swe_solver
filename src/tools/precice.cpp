#include "tools/precice.hh"
#include <iomanip>      // std::setprecision

void writeCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float time, float &time_CP, int size, int columNr){}
void restoreCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float &time, float time_CP, int size, int columNr){}

//(2)
void write2Interfoam_2D3DSupercritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  int columNr)
{
    std::cout << "executing 2d to 3d supercritical write" << '\n';
    int nY = data->l_nY; //resolution
    double dY = data->l_dY; //cell size in y-direction
    int dim = 3;
    int count = 0;

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
        count = 0;
        double u, v, h;
        for(int i = 1; i <= nY ; i++){
            h = (double)(wavePropagationBlock.getWaterHeight()[columNr][i]);
            u = (double)(wavePropagationBlock.getDischarge_hu()[columNr][i]) / h;
            v = (double)(wavePropagationBlock.getDischarge_hv()[columNr][i]) / h;
            for (int j = 0; j < nY; j++) {
                U[((i-1) + j*nY) * dim + 0] = u;
                U[((i-1) + j*nY) * dim + 1] = v;
                U[((i-1) + j*nY) * dim + 2] = 0.0;
            }
        }

        // Mintgens algorithm
        int dim2 = 2;
        double rho_water = 1000.0; // water density
        double tau[dim2 * nY * nY]{0.0}; // wall shear stress
        double ustar[dim2 * nY * nY]{0.0}; // modified velocity: sqrt(tau / rho)
        double mu = 1e-06 * rho_water; // dynamic viscosity * density (water)
        double du[dim * nY * nY]{0.0}; // velocity gradient from OF
        interface.readBlockVectorData(data->velocityGradientId, nY * nY, data->vertexIDs, du);
        for (int i = 0; i < nY; i++) {
            for (int j = 0; j < nY; j++) {
                tau[(i + j*nY) * dim + 0] =abs(mu * du[(i + j*nY) * dim + 0] / dY);
                tau[(i + j*nY) * dim + 1] =abs(mu * du[(i + j*nY) * dim + 1] / dY);
            }
        }

        count = 0;
        for (int i = 0; i < nY; i++) {
            for (int j = 0; j < nY; j++) {
                double yCoord = data->grid[(i + j*nY) * dim + 1];
                ustar[(i + j*nY) * dim2 + 0] = std::sqrt(tau[(i + j*nY) * dim2 + 0] / rho_water);
                ustar[(i + j*nY) * dim2 + 1] = std::sqrt(tau[(i + j*nY) * dim2 + 1] / rho_water);
                U[count++] = u + ustar[(i + j*nY) * dim + 0] / 0.41 * (1+std::log(yCoord/h));
                U[count++] = v + ustar[(i + j*nY) * dim + 1] / 0.41 * (1+std::log(yCoord/h));
                U[count++] = 0.0;
            }
        }

        double q3d[dim2]{};
        double q2d[dim2]{};
        double beta[dim2]{};
        for (int i = 0; i < nY; i++) {
            q2d[0] = (double)(wavePropagationBlock.getDischarge_hu()[columNr][i]);
            q2d[1] = (double)(wavePropagationBlock.getDischarge_hv()[columNr][i]);
            for (int j = 0; j < nY; j++) {
                double yCoord = data->grid[(i + j*nY) * dim + 1];      // y-coord of the IF mesh
                q3d[0] += alpha[i + j*nY] * U[(i + j*nY) * dim + 0] * yCoord;
                q3d[1] += alpha[i + j*nY] * U[(i + j*nY) * dim + 1] * yCoord;
            }
        }

        for (int n = 0; n < dim2; n++) {
            beta[n] = abs(q2d[n] / q3d[n]);
            if (beta[n] > 1) {
                while (beta[n] > 1 ) {
                    for (int i = 0; i < nY; i++) {
                        q3d[n] = 0.0;
                        for (int j = 0; j < nY; j++) {
                            double yCoord = data->grid[(i + j*nY) * dim + 1];
                            U[(i + j*nY) * dim + n] = U[(i + j*nY) * dim + n] * beta[n];
                            q3d[n] += alpha[i + j*nY] * U[(i + j*nY) * dim + n] * yCoord;
                        }
                        beta[n] = abs(q2d[n] / q3d[n]);
                    }
                }
            }
            // else if (beta[n] < 1) { //TODO missing to take care of this case, so far this breaks it all
            //     std::cout << "beta[" << n << "]" << beta[n] << '\n';
            //     while (beta[n] < 1 ) {
            //         for (int i = 0; i < nY; i++) {
            //             q3d[n] = 0.0;
            //             for (int j = 0; j < nY; j++) {
            //                 double yCoord = data->grid[(i + j*nY) * dim + 1];
            //                 U[(i + j*nY) * dim + n] = U[(i + j*nY) * dim + n] / beta[n];
            //                 q3d[n] += alpha[i + j*nY] * U[(i + j*nY) * dim + n] * yCoord;
            //             }
            //             beta[n] = q2d[n] / q3d[n];
            //         }
            //     }
            // }
        }

        double umax[dim2]{};
        for (int n = 0; n < dim2; n++) {
            for (int i = 0; i < nY; i++) {
                for (int j = 0; j < nY; j++) {
                    if (umax[n] < alpha[i*nY + j] * U[i*nY + j]) {
                        umax[n] = alpha[i*nY + j] * U[i*nY + j];
                    }
                }
                for (int j = 0; j < nY; j++) {
                    if (alpha[i*nY + j] < 0.001) {
                        U[i*nY + j] = umax[n];
                    }
                }
            }
        }
    // Writing Velocity to OF - End

    // Exchange data
    if (data->simType == twoDthreeDdsup) { // 2d to 3d supercritical TODO maybe remove this if
        interface.writeBlockScalarData(data->alphaId, nY * nY, data->vertexIDs, alpha);
    }
    interface.writeBlockVectorData(data->velocityId, nY * nY, data->vertexIDs, U);
}

//(4)
void readFromInterfoam_3D2DSupercritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
    SWE_Block1D* ghoshtBlock, int columNr)
{
    std::cout << "executing 3d to 2d supercritical read" << '\n';

    int nY = data->l_nY; //resolution
    double dY = data->l_dY; //cell size in y-direction
    int dim = 3; //3-dimension

    // read alpha
    double alpha[nY * nY]{};
    interface.readBlockScalarData(data->alphaId, nY * nY, data->vertexIDs, alpha);
    float h[nY+2] {};
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
            if (alpha[i*nY + (j-1)] > 0.01) {
                hu[j] += (float)U[(i*nY + (j-1)) * dim + 0] * (float)alpha[i*nY + (j-1)] * dY;
                hv[j] += (float)U[(i*nY + (j-1)) * dim + 1] * (float)alpha[i*nY + (j-1)] * dY;
            }
        }
    }

    SWE_Block1D newBlock{h, hu, hv, 1};
    ghoshtBlock->copyFrom(&newBlock, nY+2);
}

//(3)
void write2Interfoam_2D3DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  int columNr, std::vector<double> alpha)
{
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
            double alphatmp = alpha[(i-1) + j*nY];
            U[((i-1) + j*nY) * dim + 0] = alphatmp * u;
            U[((i-1) + j*nY) * dim + 1] = alphatmp * v;
            U[((i-1) + j*nY) * dim + 2] = 0.0;
        }
    }

    // // Mintgens algorithm
    // double U[dim * nY * nY] {}; //3D velocity holder for sending it to IF
    // double dY = data->l_dY; //cell size in y-direction
    // int count = 0;
    // double u, v, h;
    // for(int i = 1; i <= nY ; i++){
    //     h = (double)(wavePropagationBlock.getWaterHeight()[columNr][i]);
    //     u = (double)(wavePropagationBlock.getDischarge_hu()[columNr][i]) / h;
    //     v = (double)(wavePropagationBlock.getDischarge_hv()[columNr][i]) / h;
    //     for (int j = 0; j < nY; j++) {
    //         U[((i-1) + j*nY) * dim + 0] = u;
    //         U[((i-1) + j*nY) * dim + 1] = v;
    //         U[((i-1) + j*nY) * dim + 2] = 0.0;
    //     }
    // }
    //
    // // Mintegens algorithm
    // int dim2 = 2;
    // double rho_water = 1000.0; // water density
    // double tau[dim2 * nY * nY]{0.0}; // wall shear stress
    // double ustar[dim2 * nY * nY]{0.0}; // modified velocity: sqrt(tau / rho)
    // double mu = 1e-06 * rho_water; // dynamic viscosity * density (water)
    // double du[dim * nY * nY]{0.0}; // velocity gradient from OF
    // interface.readBlockVectorData(data->velocityGradientId, nY * nY, data->vertexIDs, du);
    // for (int i = 0; i < nY; i++) {
    //     for (int j = 0; j < nY; j++) {
    //         tau[(i + j*nY) * dim + 0] =abs(mu * du[(i + j*nY) * dim + 0] / dY);
    //         tau[(i + j*nY) * dim + 1] =abs(mu * du[(i + j*nY) * dim + 1] / dY);
    //     }
    // }
    //
    // count = 0;
    // for (int i = 0; i < nY; i++) {
    //     for (int j = 0; j < nY; j++) {
    //         double yCoord = data->grid[(i + j*nY) * dim + 1];
    //         ustar[(i + j*nY) * dim2 + 0] = std::sqrt(tau[(i + j*nY) * dim2 + 0] / rho_water);
    //         ustar[(i + j*nY) * dim2 + 1] = std::sqrt(tau[(i + j*nY) * dim2 + 1] / rho_water);
    //         U[count++] = u + ustar[(i + j*nY) * dim + 0] / 0.41 * (1+std::log(yCoord/h));
    //         U[count++] = v + ustar[(i + j*nY) * dim + 1] / 0.41 * (1+std::log(yCoord/h));
    //         U[count++] = 0.0;
    //     }
    // }
    //
    // double q3d[dim2]{};
    // double q2d[dim2]{};
    // double beta[dim2]{};
    // for (int i = 0; i < nY; i++) {
    //     q2d[0] = (double)(wavePropagationBlock.getDischarge_hu()[columNr][i]);
    //     q2d[1] = (double)(wavePropagationBlock.getDischarge_hv()[columNr][i]);
    //     for (int j = 0; j < nY; j++) {
    //         double yCoord = data->grid[(i + j*nY) * dim + 1];      // y-coord of the IF mesh
    //         q3d[0] += alpha[i + j*nY] * U[(i + j*nY) * dim + 0] * yCoord;
    //         q3d[1] += alpha[i + j*nY] * U[(i + j*nY) * dim + 1] * yCoord;
    //     }
    // }
    //
    // for (int n = 0; n < dim2; n++) {
    //     beta[n] = abs(q2d[n] / q3d[n]);
    //     if (beta[n] > 1) {
    //         while (beta[n] > 1 ) {
    //             for (int i = 0; i < nY; i++) {
    //                 q3d[n] = 0.0;
    //                 for (int j = 0; j < nY; j++) {
    //                     double yCoord = data->grid[(i + j*nY) * dim + 1];
    //                     U[(i + j*nY) * dim + n] = U[(i + j*nY) * dim + n] * beta[n];
    //                     q3d[n] += alpha[i + j*nY] * U[(i + j*nY) * dim + n] * yCoord;
    //                 }
    //                 beta[n] = abs(q2d[n] / q3d[n]);
    //             }
    //         }
    //     }
    //     // else if (beta[n] < 1) { //TODO missing to take care of this case, so far this breaks it all
    //     //     std::cout << "beta[" << n << "]" << beta[n] << '\n';
    //     //     while (beta[n] < 1 ) {
    //     //         for (int i = 0; i < nY; i++) {
    //     //             q3d[n] = 0.0;
    //     //             for (int j = 0; j < nY; j++) {
    //     //                 double yCoord = data->grid[(i + j*nY) * dim + 1];
    //     //                 U[(i + j*nY) * dim + n] = U[(i + j*nY) * dim + n] / beta[n];
    //     //                 q3d[n] += alpha[i + j*nY] * U[(i + j*nY) * dim + n] * yCoord;
    //     //             }
    //     //             beta[n] = q2d[n] / q3d[n];
    //     //         }
    //     //     }
    //     // }
    // }
    //
    // double umax[dim2]{};
    // for (int n = 0; n < dim2; n++) {
    //     for (int i = 0; i < nY; i++) {
    //         for (int j = 0; j < nY; j++) {
    //             if (umax[n] < alpha[i*nY + j] * U[i*nY + j]) { // get max
    //                 umax[n] = alpha[i*nY + j] * U[i*nY + j];
    //             }
    //         }
    //         for (int j = 0; j < nY; j++) {
    //             if (alpha[i*nY + j] < 0.001) {
    //                 U[i*nY + j] = umax[n];
    //             }
    //         }
    //     }
    // }

    interface.writeBlockVectorData(data->velocityId, nY * nY, data->vertexIDs, U);
}

//(3)
std::vector<double> readFromInterfoam_2D3DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  SWE_Block1D* ghoshtBlock, int columNr)
{
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

    return alphav;
}

//(5)
void write2Interfoam_3D2DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  int columNr)
{
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

//(5)
void readFromInterfoam_3D2DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
    SWE_Block1D* ghoshtBlock, int columNr)
{
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
