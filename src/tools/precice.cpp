#include "tools/precice.hh"

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
    double* alpha = new double[nY * nY]; //scalar alpha holder for sending it to IF
    for (int i = 0; i < nY; i++) {
        double height = (double)wavePropagationBlock.getWaterHeight()[columNr][i];
        for (int j = 0; j < nY; j++) {
            double yCoord = data->grid3D[count * dim + 1];      // y-coor of the IF mesh
            if(height < yCoord - 0.5 * dY ){
                alpha[count] = 0.0;
            } else if(height > yCoord + 0.5 * dY){
                alpha[count] = 1.0;
            } else{
                alpha[count] = 0.5 + (height - yCoord) / dY;	// Interpolate
            }
            count++;
        }
    }
    // Writing Alpha to OF - End


    // Writing prgh to OF - Begin
    double rho_water = 1000.0;
    double rho_air = 1.0;
    double g = 9.81;

    double* p_rgh = new double[nY * nY];
    count = 0;
    for (int i = 0; i < nY; i++) {
        double rho_mixed = rho_water * alpha[i] + rho_air * (1 - alpha[i]);
        for (int j = 0; j < nY; j++) {
            double yCoord = data->grid3D[count * dim + 1];   // y-coord of the IF mesh
            p_rgh[count] = rho_mixed * g * yCoord;           // Set p_rgh  eq 5.16 mintgen
            count++;
        }
    }
    // Writing prgh to OF - End

    // Writing Velocity to OF - Begin
    double* U = new double[dim * nY * nY]; //3D velocity holder for sending it to IF
    count = 0;
    for(int i = 1; i < nY ; i++){
        double hu = (double)(wavePropagationBlock.getDischarge_hu()[columNr][i]);
        double hv = (double)(wavePropagationBlock.getDischarge_hv()[columNr][i]);
        for (int j = 0; j < nY; j++) {
            double yCoord = data->grid3D[count * dim + 1]; // y-coord of the IF mesh
            // divide hu and hv by h for sending the velocities to interfoam
            U[count * dim + 0] = hu / yCoord;
            U[count * dim + 1] = hv / yCoord;
            U[count * dim + 2] = 0.0;
            count++;
        }
    }
    // Writing Velocity to OF - End

    // Exchange data
    interface.writeBlockScalarData(data->snd_alphaId, nY * nY, data->vertexIDs, alpha);
    interface.writeBlockScalarData(data->snd_prghId, nY * nY, data->vertexIDs, p_rgh);
    interface.writeBlockVectorData(data->snd_VelocityId, nY * nY, data->vertexIDs, U);

    delete [] alpha;
    delete [] p_rgh;
    delete [] U;

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
