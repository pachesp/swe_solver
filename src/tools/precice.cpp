#include "tools/precice.hh"
#include <iomanip>      // std::setprecision

void writeCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float time, float &time_CP, int size, int columNr){}
void restoreCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float &time, float time_CP, int size, int columNr){}

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

        // Mintegens algorithm
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
                double yCoord = data->grid3D[(i + j*nY) * dim + 1];
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
    if (data->simType == twoDthreeDdsup) { // 2d to 3d supercritical
        interface.writeBlockScalarData(data->alphaId, nY * nY, data->vertexIDs, alpha);
    }
    interface.writeBlockVectorData(data->velocityId, nY * nY, data->vertexIDs, U);
}

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

void write2Interfoam_2D3DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  int columNr)
{
    std::cout << "executing 2d to 3d subcritical write" << '\n';
    write2Interfoam_2D3DSupercritical_preCICE(interface, wavePropagationBlock, data, columNr, tempVelocity3d);
}

void readFromInterfoam_2D3DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
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


    //gradient 0 for discharges hu and hv
    float hu[nY+2] {};
    float hv[nY+2] {};
    for (int i = 0; i < nY+2; i++) {
            hu[i] = wavePropagationBlock.getDischarge_hu()[columNr][i];
            hv[i] = wavePropagationBlock.getDischarge_hv()[columNr][i];
    }

    SWE_Block1D newBlock{h, hu, hv, 1};
    ghoshtBlock->copyFrom(&newBlock, nY+2);

}

void write2Interfoam_3D2DSubcritical_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  int columNr)
{
    // double dY = data->l_dY;
    // double rho_water = 1000.0;
    // double rho_air = 1.0;
    // double rho_mixed[nY * nY];
    // double hh;
    // int dim = 3;
    // double p_rgh[nY * nY]{};
    std::cout << "executing 3d to 2d subcritical write" << '\n';
    int nY = data->l_nY; //resolution
    double g = 9.81;
    double h;
    double gh[nY * nY]{};
    for (int i = 0; i < nY; i++) {
        h = wavePropagationBlock.getWaterHeight()[columNr][i];
        // std::cout  <<fixed << setprecision(2) << "h[: " << i << "]: "<<  h  << '\n';
        for (int j = 0; j < nY; j++) {
            // double yCoord = data->grid3D[(i + j*nY) * dim + 1]; // y-coord of the IF mesh
            // rho_mixed[i + j*nY] = rho_water * alphav[i + j*nY] + rho_air * (1 - alphav[i + j*nY]);
            // std::cout <<fixed << setprecision(2) << "alhpa[" << i << "]["<< j << "]: " << alphav[i + j*nY] << '\n';
            // p_rgh[i + j*nY] = rho_mixed[i + j*nY] * g * h;           // Set p_rgh  eq 5.16 mintgen
            gh[i + j*nY] = g * h;           // Set p_rgh  eq 5.16 mintgen
        }
        // std::cout  << '\n';
    }
    // std::cout << "\nh" << '\n';
    // for (int i = 0; i < nY; i++) {
    //     std::cout <<fixed << setprecision(1) << h << ", ";
    // }
    // std::cout  << '\n';
    //
    // std::cout << "\nAlpha Vector" << '\n';
    // for (int i = 0; i < nY; i++) {
    //     for (int j = 0; j < nY; j++) {
    //         std::cout <<fixed << setprecision(2)<< abs(alphav[i*nY + j])  <<", ";
    //     }
    //     std::cout  << '\n';
    // }
    //
    // std::cout << "\nrhoMixed Vector" << '\n';
    // for (int i = 0; i < nY; i++) {
    //     for (int j = 0; j < nY; j++) {
    //         std::cout <<fixed << setprecision(2)<< rho_mixed[i*nY + j]  <<", ";
    //     }
    //     std::cout  << '\n';
    // }
    //
    // std::cout << "\nPressure" << '\n';
    // for (int i = 0; i < nY; i++) {
    //     for (int j = 0; j < nY; j++) {
    //     std::cout <<fixed << setprecision(1) << p_rgh[i*nY + j] << ", ";
    //     }
    //     std::cout  << '\n';
    // }
    // std::cout  << '\n';
    // // std::cin.get();
    interface.writeBlockScalarData(data->ghId, nY * nY, data->vertexIDs, gh);
}

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
    // std::vector<double> alphav;//(alpha, alpha+nY); //copy array to vector
    // for (int i = 0; i < nY * nY; i++) {
    //     alphav.push_back(alpha[i]);
    // }
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

    // std::cout << "SENDING THIS alpha" << '\n';
    // for (int i = 0; i < nY; i++) {
    //     for (int j = 0; j < nY; j++) {
    //         std::cout << alpha[i +j*nY]  <<", ";
    //     }
    //     std::cout  << '\n';
    // }
    //
    // std::cout << "SENDING THIS alphaV" << '\n';
    // for (int i = 0; i < nY; i++) {
    //     for (int j = 0; j < nY; j++) {
    //         std::cout << alphav[i +j*nY]  <<", ";
    //     }
    //     std::cout  << '\n';
    // }

    //gradient 0 for h
    float h[nY+2] {};
    for (int i = 0; i < nY+2; i++) {
            h[i] = wavePropagationBlock.getWaterHeight()[1][i];
    }

    SWE_Block1D newBlock{h, hu, hv, 1};
    ghoshtBlock->copyFrom(&newBlock, nY+2);
    // return alphav;
}
