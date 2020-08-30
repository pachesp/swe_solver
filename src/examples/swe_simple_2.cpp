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
  SWE_SWE_Supercritical_Right_Scenario* l_scenario = new SWE_SWE_Supercritical_Right_Scenario();

  if(simType != twoDtwoDsup){
      assert(false && "wrong type in SWE_SWE_Supercritical_Right_Scenario");
  }

  //! number of checkpoints for visualization (at each checkpoint in time, an output file is written).
  // int l_numberOfCheckPoints = 50;
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

  //***************preCICE**************************
  std::string configFileName("precice-config.xml");
  std::string solverName = "SWE_solver_right";
  SolverInterface interface(solverName, configFileName, 0, 1);
  int dimensions = interface.getDimensions();
  int meshID = interface.getMeshID("SWE_solver_right-Mesh");

  int* vertexIDs = new int[l_nY + 2];
  double* grid = new double[dimensions * (l_nY + 2)];
  int count=0;
  for (int j=0; j <= l_nY + 1 ; j++){
      grid[count++] = 0;
      grid[count++] = j;
  }

  interface.setMeshVertices(meshID, (l_nY + 2) , grid, vertexIDs);

  float precice_dt = interface.initialize();

  PreciceData* preciceData = new PreciceData();
  preciceData->heightId = interface.getDataID("heightS1", meshID);
  preciceData->huId = interface.getDataID("huS1", meshID);
  preciceData->hvId = interface.getDataID("hvS1", meshID);
  preciceData->heightGradId = interface.getDataID("heightGrad", meshID);
  preciceData->huGradId = interface.getDataID("huGrad", meshID);
  preciceData->hvGradId = interface.getDataID("hvGrad", meshID);
  preciceData->l_nY = l_nY;
  preciceData->l_dY = (double)l_dY;
  preciceData->simType = simType;
  preciceData->vertexIDs = vertexIDs;

  float time_CP;  //TODO

  //***************preCICE**************************

  // initialize the wave propagation block
  l_wavePropgationBlock.initScenario(l_originX, l_originY, l_scenario);


  //! checkpoints when output files are written.
  float* l_checkPoints = new float[l_numberOfCheckPoints+1];

  // compute the checkpoints in time
  for(int cp = 0; cp <= l_numberOfCheckPoints; cp++) {
     l_checkPoints[cp] = cp*(l_endSimulation/l_numberOfCheckPoints);
  }

  SWE_Block1D* l_leftGhostCells  = l_wavePropgationBlock.grabGhostLayer(BND_LEFT);
  // SWE_Block1D* l_leftGhostCells  = l_wavePropgationBlock.grabEdge(BND_LEFT);

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
          l_originX, l_originY );

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

  l_wavePropgationBlock.setBoundaryType(BND_LEFT, INFLOW_COUPLE);

//
// // //this apparently comes from the config file
// //   if (interface.isActionRequired(actionWriteInitialData())) {
// //     std::cout << "SWE_solver_right action write initial data" << '\n';
// //     write2SWE_SWE_Left_preCICE(interface, l_wavePropgationBlock, preciceData, l_nX+2, 1);
// //     interface.markActionFulfilled(actionWriteInitialData());
// //   }
//
  interface.initializeData();


  if (interface.isReadDataAvailable()) {
    std::cout << "precicedt = " << precice_dt <<'\n';
    std::cout << " SWE_solver_right Read data Available outside loop" << '\n';
    readFromSWE_SWE_Left_preCICE(interface, l_wavePropgationBlock, l_leftGhostCells, preciceData, l_nX +2);
  }

  //! simulation time.
  float l_t = 0.0;
  progressBar.update(l_t);
  unsigned int l_iterations = 0;
  int chkpt=1;
  while(interface.isCouplingOngoing()){

        if(interface.isActionRequired(actionWriteIterationCheckpoint())) {
            writeCheckpoint(preciceData, l_wavePropgationBlock, l_t, time_CP, l_nX+2, 1);
            interface.markActionFulfilled(actionWriteIterationCheckpoint());
          }

        // set values in ghost cells:
        l_wavePropgationBlock.setGhostLayer();

        l_wavePropgationBlock.calculateGradient(l_leftGhostCells);

        write2SWE_SWE_Left_preCICE(interface, l_wavePropgationBlock, preciceData, l_nX+2, 1);

        // reset the cpu clock
        tools::Logger::logger.resetClockToCurrentTime("Cpu");

        // approximate the maximum time step
        // TODO: This calculation should be replaced by the usage of the wave speeds occuring during the flux computation
        // Remark: The code is executed on the CPU, therefore a "valid result" depends on the CPU-GPU-synchronization.
        l_wavePropgationBlock.computeMaxTimestep();

        // compute numerical flux on each edge
        l_wavePropgationBlock.computeNumericalFluxes();

        //! maximum allowed time step width.
        // float l_maxTimeStepWidth = l_wavePropgationBlock.getMaxTimestep();
        float l_maxTimeStepWidth = 0.125;

        // update the cell values
        l_wavePropgationBlock.updateUnknowns(l_maxTimeStepWidth);


        l_maxTimeStepWidth = std::min(l_maxTimeStepWidth, precice_dt );

        precice_dt = interface.advance(l_maxTimeStepWidth);

        readFromSWE_SWE_Left_preCICE(interface, l_wavePropgationBlock, l_leftGhostCells,  preciceData,  l_nX +2);

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
          restoreCheckpoint(preciceData, l_wavePropgationBlock, l_t, time_CP, l_nX+2, 1);
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
