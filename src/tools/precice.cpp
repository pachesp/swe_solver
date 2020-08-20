#include "tools/precice.hh"
#include <iomanip>      // std::setprecision


void write2Interfoam_supercritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  int columNr, double* tempVelocity3d)
{
    int nY = data->l_nY; //resolution
    double dY = data->l_dY; //cell size in y-direction
    int dim = 3;
    int dim2 = 2;
    int count = 0;

    // Writing Alpha to OF - Begin
    double alpha[nY * nY] {};// = new double[nY * nY]; //scalar alpha holder for sending it to IF
    for (int i = 0; i < nY; i++) {
        double height = (double)wavePropagationBlock.getWaterHeight()[columNr][i];
        for (int j = 0; j < nY; j++) {
            double yCoord = data->grid3D[(i + j*nY) * dim + 1]; // y-coord of the IF mesh
            //Algorithm 1, Mintegen page 84
            if(height <= yCoord - 0.5 * dY ){
                alpha[i + j*nY] = 0.0;
            } else if(height >= yCoord + 0.5 * dY){
                 alpha[i + j*nY] = 1.0;
            } else{ // Interpolate
                alpha[i + j*nY] = 0.5 + (height - yCoord) / dY;
            }

            //Algorithm 2
            // if(height <= yCoord){
            //     alpha[i + j*nY] = 0.0;
            // } else{
            //     alpha[i + j*nY] = 1.0;
            // }

            // std::cout << "comparing yCoord: " << yCoord << " and Waterheight: " << height;
            // std::cout << " alpha is: " << alpha[i + j*nY] << '\n';
        }
    }
    // Writing Alpha to OF - End


    // Writing Velocity to OF - Begin

        //normal Algorithm
        // double U[dim * nY * nY] {}; //3D velocity holder for sending it to IF
        // count = 0;
        // int countAlpha = 0;
        // double u, v, h;
        // for(int i = 1; i <= nY ; i++){
        //     // divide hu and hv by h for sending the velocities to interfoam if alpha > 0.001
        //     if (alpha[countAlpha++] > 0.001){
        //         h = (double)(wavePropagationBlock.getWaterHeight()[columNr][i]);
        //         u = (double)(wavePropagationBlock.getDischarge_hu()[columNr][i]) / h;
        //         v = (double)(wavePropagationBlock.getDischarge_hv()[columNr][i]) / h;
        //
        //     } else {
        //         h = (double)(wavePropagationBlock.getWaterHeight()[columNr][i]);
        //         u = 0.0;
        //         v = 0.0;
        //     }
        //
        //     for (int j = 0; j < nY; j++) {
        //         U[count++] = u;
        //         U[count++] = v;
        //         U[count++] = 0.0;
        //     }
        // }

        // Mintegens algorithm
        double U[dim * nY * nY] {}; //3D velocity holder for sending it to IF
        double u, v, u_, v_, h, h_; // auxiliar variables
        double rho_water = 1000.0; // water density
        double tau[dim2 * nY * nY]{0.0}; // wall shear stress
        double ustar[dim2 * nY * nY]{0.0}; // modified velocity: sqrt(tau / rho)
        double mu = 1e-06 * rho_water; // dynamic viscosity * density (water)

        // //handle du from 3d domain
        // double du[dim2 * nY * nY]{0.0}; // velocity gradient
        // interface.readBlockVectorData(data->TempVelocity3dId, nY * nY, data->vertexIDs, tempVelocity3d);
        // for (int i = 0; i < nY; i++) {
        //     for ( int j = 0; j < nY; j++) {
        //         for (int n = 0; n < 3; n++) {
        //             std::cout << tempVelocity3d[i*nY + j + n] << ' ';
        //         }
        //         std::cout  << ",\t";
        //     }
        //     std::cout << '\n';
        // }
        //
        // std::cin.get();
        // for (int i = 0; i < nY; i++) {
        //     for (int j = 0; j < nY; j++) {
        //         if (i == 0 || i == nY || j == 0 || j == nY){
        //             du[(i*nY + j) * dim + 0] = tempVelocity3d[(i*nY + j) * dim + 0];
        //             du[(i*nY + j) * dim + 1] = tempVelocity3d[(i*nY + j) * dim + 1];
        //         }else{
        //             u = tempVelocity3d[(i*nY + j) * dim + 0];
        //             v = tempVelocity3d[(i*nY + j) * dim + 1];
        //             u_ = tempVelocity3d[((i-1)*nY + j) * dim + 0];
        //             v_ = tempVelocity3d[(i*nY + (j-1)) * dim + 1];
        //             du[(i*nY + j) * dim + 0] = u - u_;
        //             du[(i*nY + j) * dim + 1] = v - v_;
        //         }
        //         tau[(i + j*nY) * dim + 0] =abs(mu * du[(i + j*nY) * dim + 0] / dY);
        //         tau[(i + j*nY) * dim + 1] =abs(mu * du[(i + j*nY) * dim + 1] / dY);
        //     }
        //     std::cout << '\n';
        // }


        // //handle du  and tau from 2d
        double du[dim2 * nY]{0.0}; // velocity gradient
        for (int i = 0; i < nY; i++) {
            h = (double)(wavePropagationBlock.getWaterHeight()[columNr][i]);
            u = (double)(wavePropagationBlock.getDischarge_hu()[columNr][i]) / h;
            v = (double)(wavePropagationBlock.getDischarge_hv()[columNr][i]) / h;
            h_ = (double)(wavePropagationBlock.getWaterHeight()[columNr-1][i]);
            u_ = (double)(wavePropagationBlock.getDischarge_hu()[columNr-1][i]) / h_;
            v_ = (double)(wavePropagationBlock.getDischarge_hv()[columNr-1][i]) / h_;
            du[i*dim2 + 0] = u - u_;
            du[i*dim2 + 1] = v - v_;
            for (int j = 0; j < nY; j++) {
                tau[(i + j*nY) * dim2 + 0] =abs(mu * du[i*dim2 + 0] / dY);
                tau[(i + j*nY) * dim2 + 1] =abs(mu * du[i*dim2 + 1] / dY);
            }
        }

        count = 0;
        for (int i = 0; i < nY; i++) {
            for (int j = 0; j < nY; j++) {
                double yCoord = data->grid3D[(i + j*nY) * dim + 1];
                ustar[(i + j*nY) * dim2 + 0] = std::sqrt(tau[(i + j*nY) * dim2 + 0] / rho_water);
                ustar[(i + j*nY) * dim2 + 1] = std::sqrt(tau[(i + j*nY) * dim2 + 1] / rho_water);
                U[count++] = u + ustar[(i + j*nY) * dim + 0] / 0.41 * (1+std::log(yCoord/h));
                U[count++] = v + ustar[(i + j*nY) * dim + 1] / 0.41 * (1+std::log(yCoord/h));
                U[count++] = 0.0;
                // std::cout << " Ux[" << i << "][" << j << "]: " << U[(i + j*nY) * dim + 0] << ",";
                // std::cout << " Uy[" << i << "][" << j << "]: " << U[(i + j*nY) * dim + 1] << "\n";
            }
        }

        double q3d[dim2]{};
        double q2d[dim2]{};
        double beta[dim2]{};
        for (int i = 0; i < nY; i++) {
            q2d[0] = (double)(wavePropagationBlock.getDischarge_hu()[columNr][i]);
            q2d[1] = (double)(wavePropagationBlock.getDischarge_hv()[columNr][i]);
            for (int j = 0; j < nY; j++) {
                double yCoord = data->grid3D[(i + j*nY) * dim + 1];      // y-coord of the IF mesh
                q3d[0] += alpha[i + j*nY] * U[(i + j*nY) * dim + 0] * yCoord;
                q3d[1] += alpha[i + j*nY] * U[(i + j*nY) * dim + 1] * yCoord;
            }
        }

        for (int n = 0; n < dim2; n++) {
            beta[n] = q2d[n] / q3d[n];
            if (beta[n] > 1) {
                while (beta[n] > 1 ) {
                    for (int i = 0; i < nY; i++) {
                        q3d[n] = 0.0;
                        for (int j = 0; j < nY; j++) {
                            double yCoord = data->grid3D[(i + j*nY) * dim + 1];
                            U[(i + j*nY) * dim + n] = U[(i + j*nY) * dim + n] * beta[n];
                            q3d[n] += alpha[i + j*nY] * U[(i + j*nY) * dim + n] * yCoord;
                        }
                        beta[n] = q2d[n] / q3d[n];
                    }
                }
            }
            // else if (beta[n] < 1) { //TODO missing to take care of this case, so far this breaks it all
            //     std::cout << "beta[" << n << "]" << beta[n] << '\n';
            //     while (beta[n] < 1 ) {
            //         for (int i = 0; i < nY; i++) {
            //             q3d[n] = 0.0;
            //             for (int j = 0; j < nY; j++) {
            //                 double yCoord = data->grid3D[(i + j*nY) * dim + 1];
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
    interface.writeBlockScalarData(data->alphaId, nY * nY, data->vertexIDs, alpha);
    interface.writeBlockVectorData(data->velocityId, nY * nY, data->vertexIDs, U);

}

void readFromInterfoam_supercritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
    SWE_Block1D* ghoshtBlock, int columNr)
{

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
    for (int i = 1; i <= nY; i++) {
        for (int j = 0; j < nY; j++) {
            hu[j] += (float)U[(i*nY + (j-1)) * dim + 0] * (float)alpha[i*nY + (j-1)] * dY;
            hv[j] += (float)U[(i*nY + (j-1)) * dim + 1] * (float)alpha[i*nY + (j-1)] * dY;
        }
    }

    SWE_Block1D newBlock{h, hu, hv, 1};
    ghoshtBlock->copyFrom(&newBlock, nY+2);

}

void write2Interfoam_subcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  int columNr, double* tempVelocity3d)
{
    // // Writing prgh to OF - Begin
    // double rho_water = 1000.0;
    // double rho_air = 1.0;
    // double g = 9.81;
    //
    // double p_rgh[nY * nY]{};
    // count = 0;
    // for (int i = 1; i <= nY; i++) {
    //     double rho_mixed = rho_water * alpha[i] + rho_air * (1 - alpha[i]);
    //     for (int j = 1; j <= nY; j++) {
    //         double yCoord = data->grid3D[count * dim + 1];   // y-coord of the IF mesh
    //         p_rgh[count] = rho_mixed * g * yCoord;           // Set p_rgh  eq 5.16 mintgen
    //         count++;
    //     }
    // }
    // // Writing prgh to OF - End

    // interface.writeBlockScalarData(data->prghId, nY * nY, data->vertexIDs, p_rgh);
}

void readFromInterfoam_subcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  SWE_Block1D* ghoshtBlock, int columNr)
{
    int nY = data->l_nY; //resolution
    // double dY = data->l_dY; //cell size in y-direction
    int dim = 3; //3-dimension

    double tempVelocity3d[dim * nY * nY];
    interface.readBlockVectorData(data->tempVelocity3dId, nY * nY, data->vertexIDs, tempVelocity3d);

    for (int i = 0; i < nY; i++) {
        for (int j = 0; j < nY; j++) {
            std::cout << "v:" << tempVelocity3d[i*nY + j]  << ", ";
        }
        std::cout  << '\n';
    }

    // double alpha[dim * nY * nY];
    // interface.readBlockVectorData(data->alphaId, nY * nY, data->vertexIDs, alpha);
    //
    // for (int i = 0; i < nY; i++) {
    //     for (int j = 0; j < nY; j++) {
    //         std::cout << "alpha:" << alpha[i*nY + j]  << ", ";
    //     }
    //     std::cout  << '\n';
    // }

}

void writeCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float time, float &time_CP, int size, int columNr){}
void restoreCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float &time, float time_CP, int size, int columNr){}
