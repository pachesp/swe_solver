#include "tools/precice.hh"
#include <iomanip>      // std::setprecision

// void write_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data, int size, int columNr){
//
//     for(int i = 0; i < size ; i++){ //Data stored column-wise
//       data->snd_height_db[i] = (double)wavePropagationBlock.getWaterHeight()[columNr][i];
//       data->snd_hu_db[i] = (double)wavePropagationBlock.getDischarge_hu()[columNr][i];
//       data->snd_hv_db[i] = (double)wavePropagationBlock.getDischarge_hv()[columNr][i];
//     }
//
//     interface.writeBlockScalarData(data->snd_heightId, size, data->vertexIDs, data->snd_height_db);
//     interface.writeBlockScalarData(data->snd_huId, size, data->vertexIDs, data->snd_hu_db);
//     interface.writeBlockScalarData(data->snd_hvId, size, data->vertexIDs, data->snd_hv_db);
// }

// void writeGradient_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data, int size, int columNr){
//
//     for(int i = 0; i < size ; i++){ //Data stored column-wise
//
//       data->snd_height_db[i] = (double)wavePropagationBlock.hGrad[i];
//       data->snd_hu_db[i] = (double)wavePropagationBlock.huGrad_SWE[i];
//       data->snd_hv_db[i] = (double)wavePropagationBlock.hvGrad_SWE[i];
//     }
//
//     interface.writeBlockScalarData(data->snd_heightId, size, data->vertexIDs, data->snd_height_db);
//     interface.writeBlockScalarData(data->snd_huId, size, data->vertexIDs, data->snd_hu_db);
//     interface.writeBlockScalarData(data->snd_hvId, size, data->vertexIDs, data->snd_hv_db);
// }

// void read_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock,
//                   SWE_Block1D* ghoshtBlock, PreciceData *data, int size, int columNr){
//
//     interface.readBlockScalarData(data->recv_prghId, size, data->vertexIDs, data->recv_prgh_db);
//     // interface.readBlockScalarData(data->recv_huId, size, data->vertexIDs, data->recv_hu_db);
//     // interface.readBlockScalarData(data->recv_hvId, size, data->vertexIDs, data->recv_hv_db);
//
//     float* h = doublePointer2floatPointer(data->recv_prgh_db, size);
//     float* hu = doublePointer2floatPointer(data->recv_hu_db, size);
//     float* hv = doublePointer2floatPointer(data->recv_hv_db, size);
//
//     // Data sent by left neighbour from precice
//     SWE_Block1D* newBlock = new SWE_Block1D{h, hu, hv, 1};
//     ghoshtBlock->copyFrom(newBlock, size );
//
//     if(newBlock) delete newBlock;
//     if(h) delete h;
//     if(hu) delete hu;
//     if(hv) delete hv;
// }
//
// void readGradient_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock,
//                    PreciceData *data, int size, int columNr){
//
//     interface.readBlockScalarData(data->recv_prghId, size, data->vertexIDs, data->recv_prgh_db);
//     interface.readBlockScalarData(data->recv_huId, size, data->vertexIDs, data->recv_hu_db);
//     interface.readBlockScalarData(data->recv_hvId, size, data->vertexIDs, data->recv_hv_db);
//
//     for (int i = 0; i < size; i++) {
//         wavePropagationBlock.hGrad[i] = data->recv_prgh_db[i];
//         wavePropagationBlock.huGrad_SWE[i] = data->recv_hu_db[i];
//         wavePropagationBlock.hvGrad_SWE[i] = data->recv_hv_db[i];
//     }
//
// }

void writeCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float time, float &time_CP, int size, int columNr){

    if(data->CP_height_f2d) delete data->CP_height_f2d;
    if(data->CP_hu_f2d) delete data->CP_hu_f2d;
    if(data->CP_hv_f2d) delete data->CP_hv_f2d;

    time_CP = time;

    data->CP_height_f2d = new Float2D(wavePropagationBlock.getWaterHeight(), false);
    data->CP_hu_f2d = new Float2D(wavePropagationBlock.getDischarge_hu(), false);
    data->CP_hv_f2d = new Float2D(wavePropagationBlock.getDischarge_hv(), false);

}

void restoreCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float &time, float time_CP, int size, int columNr){
    time = time_CP;

    for(int j = 0; j < size ; j++){
        for(int i = 0; i < size ; i++){

            wavePropagationBlock.getWaterHeight()[i][j]  =  (*(data->CP_height_f2d))[i][j];
            wavePropagationBlock.getDischarge_hu()[i][j] =  (*(data->CP_hu_f2d))[i][j];
            wavePropagationBlock.getDischarge_hv()[i][j] =  (*(data->CP_hv_f2d))[i][j];
        }
    }
}

//---------------------------To interfoam--------------------------

void write2Interfoam_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  int columNr)
{
    int nY = data->l_nY; //resolution
    double dY = data->l_dY; //cell size in y-direction
    int dim = 3; //3-dimension
    int count = 0;

    // Writing Alpha to OF - Begin
    double alpha[nY * nY] {};// = new double[nY * nY]; //scalar alpha holder for sending it to IF
    for (int i = 0; i < nY; i++) {
        double height = (double)wavePropagationBlock.getWaterHeight()[columNr][i];
        for (int j = 0; j < nY; j++) {
            double yCoord = data->grid3D[(i + j*nY) * dim + 1];      // y-coord of the IF mesh
            //Algorithm 1, Mintegen page 84
            if(height <= yCoord - 0.5 * dY ){
                alpha[i + j*nY] = 0.0;
            } else if(height >= yCoord + 0.5 * dY){
                 alpha[i + j*nY] = 1.0;
            } else{
                alpha[i + j*nY] = 0.5 + (height - yCoord) / dY;	// Interpolate
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

    // Writing Velocity to OF - Begin
    double U[dim * nY * nY] {}; //3D velocity holder for sending it to IF
    count = 0;
    int countAlpha = 0;
    double u, v;
    for(int i = 1; i <= nY ; i++){
        // divide hu and hv by h for sending the velocities to interfoam if alpha > 0.001
        if (alpha[countAlpha++] > 0.001){ //TODO change to the Algorithm from mintgen
            u = (double)(wavePropagationBlock.getDischarge_hu()[columNr][i] / wavePropagationBlock.getWaterHeight()[columNr][i]);
            v = (double)(wavePropagationBlock.getDischarge_hv()[columNr][i] / wavePropagationBlock.getWaterHeight()[columNr][i]);
        } else {
            u = 0.0;
            v = 0.0;
        }

        for (int j = 0; j < nY; j++) {
            U[count++] = u;
            // std::cout << "Ux[" << i << "][" << j << "] " << U[count] << ',';
            U[count++] = v;
            // std::cout << "Uy[" << i << "][" << j << "] " << U[count] << ',';
            U[count++] = 0.0;
            // std::cout << "Uz[" << i << "][" << j << "] " << U[count];
        }
        // std::cout  << '\n';
    }
    // Writing Velocity to OF - End

    // Exchange data
    interface.writeBlockScalarData(data->alphaID, nY * nY, data->vertexIDs, alpha);
    // interface.writeBlockScalarData(data->prghID, nY * nY, data->vertexIDs, p_rgh);
    interface.writeBlockVectorData(data->velocityID, nY * nY, data->vertexIDs, U);

}


void readFromInterfoam_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  SWE_Block1D* ghoshtBlock, int columNr)
{

    int nY = data->l_nY; //resolution
    double dY = data->l_dY; //cell size in y-direction
    int dim = 3; //3-dimension

    double alpha[nY * nY]{};
    interface.readBlockScalarData(data->alphaID, nY * nY, data->vertexIDs, alpha);
    float h[nY+2] {};
    for (int i = 0; i < nY; i++) {
        for (int j = 1; j <= nY; j++) {
            h[j] += (float)alpha[i*nY + (j-1)] * dY;
            // std::cout << "alpha:" << alpha[i*nY + (j-1)] << " dY: " << dY << '\n';
            // std::cout << "h[" << i << "][" << j << "]= " << h[j] << '\n';
        }
        // std::cout  << '\n';
    }

    // //Debug
    // for (int j = 0; j < nY + 2; j++) {
    //     std::cout << "H: " << h[j] << '\n';
    // }

    double U[dim * nY * nY]{};
    interface.readBlockVectorData(data->velocityID, nY * nY, data->vertexIDs, U);
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


// void storeData_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data, int columNr, int size){
//     for(int i = 0; i < size ; i++){
//       data->snd_height_db[i] = wavePropagationBlock.getWaterHeight().float2D2doublePointer()[columNr * size + i];
//       data->snd_hu_db[i] = wavePropagationBlock.getDischarge_hu().float2D2doublePointer()[columNr * size + i];
//       data->snd_hv_db[i] = wavePropagationBlock.getDischarge_hv().float2D2doublePointer()[columNr * size + i];
//     }
// }
//
// void writeData_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data, int columNr, int size){
//     interface.writeBlockScalarData(data->snd_heightId, size, data->vertexIDs, data->snd_height_db);
//     interface.writeBlockScalarData(data->snd_huId, size, data->vertexIDs, data->snd_hu_db);
//     interface.writeBlockScalarData(data->snd_hvId, size, data->vertexIDs, data->snd_hv_db);
// }


// std::cout << "\n\nCOORDINATES " << "\n";
// count = 0;
//     for (int i = 0; i <= nY-1 ; i++) {
//         for (int j = 0; j <= nY-1; j++) {
//             std::cout << std::fixed << std::setprecision(2)<< data->grid3D[count * dim + 1] << ", ";
//             count++;
//         }
//         std::cout << '\n';
//     }
//
// std::cout << "\n\nWaterHeight" << "\n";
//     for (int i = 1; i <= nY ; i++) {
//             std::cout << std::fixed << std::setprecision(2)
//             << (double)wavePropagationBlock.getWaterHeight()[columNr][i] << ", ";
//             count++;
//         }
// std::cout <<  '\n';
//
// std::cout << "\n\nALPHA" << "\n";
// for (int i = 0; i <= nY-1; i++) {
//     for (int j = 0; j <= nY-1; j++) {
//         std::cout << std::fixed << std::setprecision(2)<< alpha[i*nY + j] << ", ";
//     }
//     std::cout << '\n';
// }
