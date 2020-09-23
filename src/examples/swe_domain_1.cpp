/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *         Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._Michael_Bader)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * Basic setting of SWE, which uses a wave propagation solver and an artificial or ASAGI scenario on a single block.
 */
#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>

#include "blocks/SWE_WavePropagationBlock.hh"
#include "writer/VtkWriter.hh"
#include "scenarios/SWE_simple_scenarios.hh"

#include "tools/args.hh"
#include "tools/help.hh"
#include "tools/Logger.hh"
#include "tools/ProgressBar.hh"

// #define deb

//precice
#include "precice/SolverInterface.hpp"
using namespace precice;
using namespace precice::constants;

#include "tools/precice.hh"

//Add second drop fall. only for 2d 3d supercritical case
float secondDrop(float x, float y){
    //TODO avoid parameters hardcoding
    float side = 10.f;
    float offsetX = 0.f;
    float offsetY = 0.f;
    float a = x - ((side * 0.5) + offsetX);
    float b = y - ((side * 0.5) + offsetY);
    bool circ = sqrt(a * a + b * b) < (side * 0.1);
    return circ ? 20.f: 5.0f;
}

/**
 * Main program for the simulation on a single SWE_WavePropagationBlock.
 */
int main( int argc, char** argv ) {
  /**
   * Initialization.
   */
  // Parse command line parameters
  tools::Args args;

  args.addOption("grid-size-x", 'x', "Number of cells in x direction");
  args.addOption("grid-size-y", 'y', "Number of cells in y direction");
  args.addOption("simType", 't', "twoDthreeDdsup = 2 \n twoDthreeDdsub = 3 \n threeDtwoDdsup = 4 \n threeDtwoDdsub =5");
  args.addOption("output-basepath", 'o', "Output base file name");

  tools::Args::Result ret = args.parse(argc, argv);

  switch (ret)
  {
  case tools::Args::Error:
	  return 1;
  case tools::Args::Help:
	  return 0;
  default:
    break;
  }

  //! number of grid cells in x- and y-direction.
  int l_nX, l_nY;

  //whether simulation is swe-of of-swe swe-swe sub or supercritical
  size_t simType;

  //! l_baseName of the plots.
  std::string l_baseName;

  // read command line parameters
  l_nX = args.getArgument<int>("grid-size-x");
  l_nY = args.getArgument<int>("grid-size-y");
  simType = args.getArgument<size_t>("simType");
  l_baseName = args.getArgument<std::string>("output-basepath");

  // create a simple artificial scenario
  SWE_Scenario* l_scenario = nullptr;

    if (simType == twoDtwoDsup || simType == twoDtwoDsub) { // if 2d to 2d supercritical or supercritical (0 and 1)
      std::cout << "Executing scenario 2d-2d supercritical or subcritical" << '\n';
        l_scenario = new SWE_SWE_Supercritical_Left_Scenario();
    }
    else if(simType == twoDthreeDdsup){ // if 2d to 3d supercritical (2)
      std::cout << "Executing scenario 2d-3d supercritical" << '\n';
        l_scenario = new SWE_OF_Supercritical_Scenario();
    }
    else if(simType == twoDthreeDdsub){ // if 2d to 3d subcritical (3)
        std::cout << "Executing scenario 2d-3d subcritical" << '\n';
        l_scenario = new SWE_OF_Subcritical_Scenario();
    }
    else if(simType == threeDtwoDdsup){ // if 3d to 2d supercritical (4)
      std::cout << "Executing scenario 3d-2d supercritical" << '\n';
        l_scenario = new OF_SWE_Supercritical_Scenario();
    }
    else if(simType == threeDtwoDdsub){ // if 3d to 2d subcritical (5)
      std::cout << "Executing scenario 3d-2d subcritical" << '\n';
        l_scenario = new OF_SWE_Subcritical_Scenario();
    }
     else{
        assert(false && "unidentified scenario");
    }

  //! number of checkpoints for visualization (at each checkpoint in time, an output file is written).
  int l_numberOfCheckPoints = l_scenario->numberOfCheckpoints();

  //! time when the simulation ends.
  float l_endSimulation = l_scenario->endSimulation();

  //! size of a single cell in x- and y-direction
  float l_dX, l_dY;

  // compute the size of a single cell
  l_dX = (l_scenario->getBoundaryPos(BND_RIGHT) - l_scenario->getBoundaryPos(BND_LEFT) )/l_nX;
  l_dY = (l_scenario->getBoundaryPos(BND_TOP) - l_scenario->getBoundaryPos(BND_BOTTOM) )/l_nY;

  // create a single wave propagation block
  SWE_WavePropagationBlock l_wavePropgationBlock(l_nX,l_nY,l_dX,l_dY);

  //! origin of the simulation domain in x- and y-direction
  float l_originX, l_originY;

  // get the origin from the scenario
  l_originX = l_scenario->getBoundaryPos(BND_LEFT);
  l_originY = l_scenario->getBoundaryPos(BND_BOTTOM);

  //! maximum allowed time step width.
  float l_maxTimeStepWidth  = l_scenario->maxTimeStepWidth();

  // holds time at checkpoints
  float time_CP;

  // initialize the wave propagation block
  l_wavePropgationBlock.initScenario(l_originX, l_originY, l_scenario);

  //! checkpoints when output files are written.
  float* l_checkPoints = new float[l_numberOfCheckPoints+1];

  // compute the checkpoints in time
  for(int cp = 0; cp <= l_numberOfCheckPoints; cp++) {
     l_checkPoints[cp] = cp*(l_endSimulation/l_numberOfCheckPoints);
  }

  SWE_Block1D* l_leftGhostCells  = l_wavePropgationBlock.grabGhostLayer(BND_LEFT);
  SWE_Block1D* l_rightGhostCells  = l_wavePropgationBlock.grabGhostLayer(BND_RIGHT);

  // Init fancy progressbar
  tools::ProgressBar progressBar(l_endSimulation);

  // write the output at time zero
  tools::Logger::logger.printOutputTime((float) 0.);
  progressBar.update(0.);

  std::string l_fileName = generateBaseFileName(l_baseName,0,0);
  //boundary size of the ghost layers
  io::BoundarySize l_boundarySize = {{1, 1, 1, 1}};

  // consturct a VtkWriter
  io::VtkWriter l_writer( l_fileName,
		  l_wavePropgationBlock.getBathymetry(),
		  l_boundarySize,
		  l_nX, l_nY,
		  l_dX, l_dY,
          simType,
          l_originX, l_originY
      );


  // Write zero time step
  l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                          l_wavePropgationBlock.getDischarge_hu(),
                          l_wavePropgationBlock.getDischarge_hv(),
                          (float) 0.);

  /**
   * Simulation.
   */
  // print the start message and reset the wall clock time
  progressBar.clear();
  tools::Logger::logger.printStartMessage();
  tools::Logger::logger.initWallClockTime(time(NULL));

  // +++++++++++++++++ preCICE config - BEGIN +++++++++++++++++
  std::string configFileName("precice-config.xml");
  std::string solverName = "SWE_solver";
  SolverInterface interface(solverName, configFileName, 0, 1);
  int dimensions = interface.getDimensions();
  int meshID = interface.getMeshID("SWE_solver-Mesh");
  int count=0;
  int meshSize;
  int* vertexIDs = nullptr;
  double* grid = nullptr;
  PreciceData* preciceData = new PreciceData{};

if(simType == twoDtwoDsup || simType == twoDtwoDsub)
{
      vertexIDs = new int[l_nY +2];
      grid = new double[dimensions * (l_nY+2)];
      for (int j=0; j <= l_nY + 1 ; j++){
        grid[count++] = 0;
        grid[count++] = j;
      }

      meshSize = l_nY+2;

      preciceData->heightId = interface.getDataID("heightS1", meshID);
      preciceData->huId = interface.getDataID("huS1", meshID);
      preciceData->hvId = interface.getDataID("hvS1", meshID);
      preciceData->heightGradId = interface.getDataID("heightGrad", meshID);
      preciceData->huGradId = interface.getDataID("huGrad", meshID);
      preciceData->hvGradId = interface.getDataID("hvGrad", meshID);

}
else{

  int alphaId = interface.getDataID("Alpha", meshID);
  int ghId = interface.getDataID("Gh", meshID);
  int velocityId = interface.getDataID("Velocity", meshID);
  int velocityGradientId = interface.getDataID("VelocityGradient", meshID);

  //set coordinates on face centres of the interfoam left boundary
  vertexIDs = new int[l_nY * l_nY];
  grid = new double[dimensions * l_nY * l_nY];

  double xBoundary;
  if (simType == twoDthreeDdsup || simType == twoDthreeDdsub) { //if 2d to 3d
      xBoundary = 10 + l_dX;
  } else if (simType == threeDtwoDdsup || simType == threeDtwoDdsub) { //if 3d to 2d
      xBoundary = 10 - l_dX;
  } else {
      assert(false);
  }

  for (int i = 1; i <= l_nY; i++){
      for (int j = 1; j <= l_nY; j++){
        grid[count++] = xBoundary;                  // x
        grid[count++] = 0 + (i - 0.5) * l_dY;       // y
        grid[count++] = 0 + (j - 0.5) * l_dY;       // z
      }
  }

  meshSize = l_nY * l_nY;
  preciceData->alphaId = alphaId;
  preciceData->ghId = ghId;
  preciceData->velocityId = velocityId;
  preciceData->velocityGradientId = velocityGradientId;

}

preciceData->l_nY = l_nY;
preciceData->l_dY = (double)l_dY;
preciceData->simType = simType;
preciceData->grid = grid;
preciceData->vertexIDs = vertexIDs;

interface.setMeshVertices(meshID, meshSize , grid, vertexIDs);

float precice_dt = interface.initialize();

 // +++++++++++++++++ preCICE config - END +++++++++++++++++

//overload write and read functions according to the flow regime
PreciceExchange* preciceExchange = nullptr;
if (simType == twoDtwoDsup || simType == twoDtwoDsub){
  std::cout << "Data exchange 2d to 2d supercritical or subcritical" << '\n';
  preciceExchange = new SWE_SWE_SuperCritical_Left_Exchange{interface, l_wavePropgationBlock,
     preciceData, l_nY, l_leftGhostCells};
}
else if (simType == twoDthreeDdsup) {
  std::cout << "Data exchange 2d to 3d supercritical" << '\n';
  preciceExchange = new SWE_OF_SuperCritical_Exchange{interface, l_wavePropgationBlock,
     preciceData, l_nY};
}
else if (simType == twoDthreeDdsub) {
  std::cout << "executing 2d to 3d subcritical" << '\n';
  preciceExchange = new SWE_OF_SubCritical_Exchange{interface, l_wavePropgationBlock,
     preciceData, l_nY, l_rightGhostCells};
}
else if (simType == threeDtwoDdsup) {
  std::cout << "executing 3d to 2d supercritical" << '\n';
  preciceExchange = new OF_SWE_SuperCritical_Exchange{interface, l_wavePropgationBlock,
     preciceData, 0, l_leftGhostCells};
}
else if (simType == threeDtwoDdsub) {
  std::cout << "executing 3d to 2d subcritical" << '\n';
  preciceExchange = new OF_SWE_SubCritical_Exchange{interface, l_wavePropgationBlock,
     preciceData, 1, l_leftGhostCells};
}
else{
  assert(false && "wrong simulation type in preciceExchange");
}

 if(simType == twoDtwoDsup || simType == twoDtwoDsub){
       if (interface.isActionRequired(actionWriteInitialData())) {
         std::cout << "solver1 action write initial data" << '\n';
         preciceExchange->write();
         interface.markActionFulfilled(actionWriteInitialData());
       }
  }

  interface.initializeData();

  //! simulation time.
  float l_t = 0.0;
  progressBar.update(l_t);
  unsigned int l_iterations = 0;
  int chkpt=1;
  bool setSecondDrop = true;
  std::vector<double> alphav;
  while(interface.isCouplingOngoing()){

      //add second droplet after 2 seconds. Only for 2D-3D supercritical
      if (l_t > 2.0 && setSecondDrop && simType == twoDthreeDdsup) {
          l_wavePropgationBlock.setWaterHeight(secondDrop);
          setSecondDrop = false;
     }

      if(interface.isActionRequired(actionWriteIterationCheckpoint())) {
        writeCheckpoint(preciceData, l_wavePropgationBlock, l_t, time_CP, l_nX+2, l_nY);
        interface.markActionFulfilled(actionWriteIterationCheckpoint());
      }

      // read precice values
      preciceExchange->read();

      // set values in ghost cells:
      l_wavePropgationBlock.setGhostLayer();

      // reset the cpu clock
      tools::Logger::logger.resetClockToCurrentTime("Cpu");

      // approximate the maximum time step
      // TODO: This calculation should be replaced by the usage of the wave speeds occuring during the flux computation
      // Remark: The code is executed on the CPU, therefore a "valid result" depends on the CPU-GPU-synchronization.
      l_wavePropgationBlock.computeMaxTimestep();

      // compute numerical flux on each edge
      l_wavePropgationBlock.computeNumericalFluxes();

       // l_maxTimeStepWidth = l_wavePropgationBlock.getMaxTimestep();

      // update the cell values
      l_wavePropgationBlock.updateUnknowns(l_maxTimeStepWidth);

      // write precice values
      preciceExchange->write();

      l_maxTimeStepWidth = std::min(l_maxTimeStepWidth, precice_dt);

      precice_dt = interface.advance(l_maxTimeStepWidth);

      // update the cpu time in the logger
      tools::Logger::logger.updateTime("Cpu");

      // update simulation time with time step width.
      l_t += l_maxTimeStepWidth;
      l_iterations++;

      // print the current simulation time
      progressBar.clear();
      tools::Logger::logger.printSimulationTime(l_t);
      progressBar.update(l_t);

      if (interface.isActionRequired(actionReadIterationCheckpoint())) {
        restoreCheckpoint(preciceData, l_wavePropgationBlock, l_t, time_CP, l_nX+2, l_nY);
        interface.markActionFulfilled(actionReadIterationCheckpoint());
      }else{
        if(l_t >= l_checkPoints[chkpt] - l_maxTimeStepWidth && l_t < l_checkPoints[chkpt] + l_maxTimeStepWidth){
            progressBar.clear();
            tools::Logger::logger.printOutputTime(l_t);
            progressBar.update(l_t);
            l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                                    l_wavePropgationBlock.getDischarge_hu(),
                                    l_wavePropgationBlock.getDischarge_hv(),
                                    l_t);
            chkpt++;
        }
    }
  }

  interface.finalize();

  /**
   * Finalize.
   */
  // write the statistics message
  progressBar.clear();
  tools::Logger::logger.printStatisticsMessage();

  // print the cpu time
  tools::Logger::logger.printTime("Cpu", "CPU time");

  // print the wall clock time (includes plotting)
  tools::Logger::logger.printWallClockTime(time(NULL));

  // printer iteration counter
  tools::Logger::logger.printIterationsDone(l_iterations);

  return 0;
}
