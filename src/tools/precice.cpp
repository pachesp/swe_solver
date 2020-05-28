#include "tools/precice.hh"

void write_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data, int size, int columNr){

    for(int i = 0; i < size ; i++){ //Data stored column-wise
      data->snd_height_db[i] = wavePropagationBlock.getWaterHeight().float2D2doublePointer()[columNr * size + i];
      data->snd_hu_db[i] = wavePropagationBlock.getDischarge_hu().float2D2doublePointer()[columNr * size + i];
      data->snd_hv_db[i] = wavePropagationBlock.getDischarge_hv().float2D2doublePointer()[columNr * size + i];
    }

    interface.writeBlockScalarData(data->snd_heightId, size, data->vertexIDs, data->snd_height_db);
    interface.writeBlockScalarData(data->snd_huId, size, data->vertexIDs, data->snd_hu_db);
    interface.writeBlockScalarData(data->snd_hvId, size, data->vertexIDs, data->snd_hv_db);
}

void read_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock,
                  SWE_Block1D* ghoshtBlock, PreciceData *data, int size, int columNr){

    interface.readBlockScalarData(data->recv_heightId, size, data->vertexIDs, data->recv_height_db);
    interface.readBlockScalarData(data->recv_huId, size, data->vertexIDs, data->recv_hu_db);
    interface.readBlockScalarData(data->recv_hvId, size, data->vertexIDs, data->recv_hv_db);

    // Data sent by left neighbour from precice
    SWE_Block1D* newBlock = new SWE_Block1D{ doublePointer2floatPointer(data->recv_height_db, size),
                                doublePointer2floatPointer(data->recv_hu_db, size),
                                doublePointer2floatPointer(data->recv_hv_db, size), 1
                              };
    ghoshtBlock->copyFrom(newBlock, size );
    delete newBlock;
}

void writeCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float time, float &time_CP, int size, int columNr){
  time_CP = time;
  for(int i = 0; i < size ; i++){
    data->CP_height_db[i] = wavePropagationBlock.getWaterHeight().float2D2doublePointer()[columNr * size + i];
    data->CP_hu_db[i] = wavePropagationBlock.getDischarge_hu().float2D2doublePointer()[columNr * size + i];
    data->CP_hv_db[i] = wavePropagationBlock.getDischarge_hv().float2D2doublePointer()[columNr * size + i];
  }
}

void restoreCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float &time, float time_CP, int size, int columNr){
  time = time_CP;
  for(int i = 0; i < size ; i++){
    wavePropagationBlock.getWaterHeight().float2D2doublePointer()[columNr * size + i] = data->CP_height_db[i];
    wavePropagationBlock.getDischarge_hu().float2D2doublePointer()[columNr * size + i] = data->CP_hu_db[i];
    wavePropagationBlock.getDischarge_hv().float2D2doublePointer()[columNr * size + i] = data->CP_hv_db[i];
  }
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
