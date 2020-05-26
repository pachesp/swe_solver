#include "tools/precice.hh"

void write_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data, int columNr, int size){


    for(int i = 0; i < size ; i++){
      data->snd_height_db[i] = wavePropagationBlock.getWaterHeight().float2D2doublePointer()[columNr * size + i];
      data->snd_hu_db[i] = wavePropagationBlock.getDischarge_hu().float2D2doublePointer()[columNr * size + i];
      data->snd_hv_db[i] = wavePropagationBlock.getDischarge_hv().float2D2doublePointer()[columNr * size + i];
    }

     // std::cout << "sending column "<< columNr << " to neighbour" << '\n';
     // for(int i = 0; i <size; i++){
     //   std::cout << data->snd_height_db[i] << '\n';
     // }

    interface.writeBlockScalarData(data->snd_heightId, size, data->vertexIDs, data->snd_height_db);
    interface.writeBlockScalarData(data->snd_huId, size, data->vertexIDs, data->snd_hu_db);
    interface.writeBlockScalarData(data->snd_hvId, size, data->vertexIDs, data->snd_hv_db);

  //   interface.markActionFulfilled(actionWriteInitialData());
  // }

}


void storeData_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data, int columNr, int size){

    for(int i = 0; i < size ; i++){
      data->snd_height_db[i] = wavePropagationBlock.getWaterHeight().float2D2doublePointer()[columNr * size + i];
      data->snd_hu_db[i] = wavePropagationBlock.getDischarge_hu().float2D2doublePointer()[columNr * size + i];
      data->snd_hv_db[i] = wavePropagationBlock.getDischarge_hv().float2D2doublePointer()[columNr * size + i];
    }

    std::cout << "sending column "<< columNr << " to neighbour" << '\n';
    for(int i = 0; i <size; i++){
      std::cout << data->snd_height_db[i] << '\n';
    }

}

void writeData_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data, int columNr, int size){

    std::cout << "sending column Corrobore "<< columNr << " to neighbour" << '\n';
    for(int i = 0; i <size; i++){
      std::cout << data->snd_height_db[i] << '\n';
    }

    interface.writeBlockScalarData(data->snd_heightId, size, data->vertexIDs, data->snd_height_db);
    interface.writeBlockScalarData(data->snd_huId, size, data->vertexIDs, data->snd_hu_db);
    interface.writeBlockScalarData(data->snd_hvId, size, data->vertexIDs, data->snd_hv_db);


}

void read_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock,
                  SWE_Block1D* ghoshtBlock, SWE_Block1D* newBlock, PreciceData *data, int columNr, int size){

  // if (interface.isReadDataAvailable()) {

    interface.readBlockScalarData(data->recv_heightId, size, data->vertexIDs, data->recv_height_db);
    interface.readBlockScalarData(data->recv_huId, size, data->vertexIDs, data->recv_hu_db);
    interface.readBlockScalarData(data->recv_hvId, size, data->vertexIDs, data->recv_hv_db);

    // std::cout << "receiving column "<< columNr <<" from neighbour" <<'\n';
    // for(int i = 0; i <size; i++){
    //   std::cout << data->recv_height_db[i] << '\n';
    // }

    // Data sent by left neighbour from precice
    newBlock = new SWE_Block1D{ doublePointer2floatPointer(data->recv_height_db, size),
                                doublePointer2floatPointer(data->recv_hu_db, size),
                                doublePointer2floatPointer(data->recv_hv_db, size), 1 };

    // deepSWE_Block1DCopy(newBlock, ghoshtBlock, l_nY+2);
    ghoshtBlock->copyFrom(newBlock, size );
    delete newBlock;

  // }
}
