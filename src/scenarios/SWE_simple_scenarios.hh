/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema, Tobias Schnabel
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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
 * TODO
 */

#ifndef __SWE_SIMPLE_SCENARIOS_H
#define __SWE_SIMPLE_SCENARIOS_H

#include <cmath>

#include "SWE_Scenario.hh"

/**
 * Scenario "Radial Dam Break":
 * elevated water in the center of the domain
 */
class SWE_RadialDamBreakScenario : public SWE_Scenario {
  public:
    virtual float getBathymetry(float x, float y) {
       return 0.f;
    };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) {return OUTFLOW;};

    virtual float getWaterHeight(float x, float y, float offsetX, float offsetY) {
       return ( sqrt( (x-((side * 0.5) + offsetX)) * (x-((side * 0.5) + offsetX)) +
                      (y-((side * 0.5) + offsetY)) * (y-((side * 0.5) + offsetY)) ) < side * 0.1 ) ? 15.f: 10.0f;
    };

    virtual float getBoundaryPos(BoundaryEdge edge) {
        if (edge==BND_LEFT || edge==BND_BOTTOM)
           return 0.0f;
        else
           return side;
     };
};

// 0 and 1 (supercritial and subcritical)
class SWE_SWE_Supercritical_Left_Scenario : public SWE_RadialDamBreakScenario{
  public:
    float side = 1000.f;

    virtual int numberOfCheckpoints(){return 50; };

    virtual float maxTimeStepWidth() {return 0.125;};

    virtual float endSimulation() { return 60.f; }

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) {
        if (edge == BND_RIGHT) {
            return OUT_INFLOW_COUPLE;
        }else{
            return OUTFLOW;
        }
    };

    virtual float getWaterHeight(float x, float y, float offsetX = 0, float offsetY = 0) {
        float b = (y-((side * 0.5)+offsetY));
        bool circ_1  = false;
        float a_0 = (x-((side * 0.875) +offsetX));
        bool circ_2 = sqrt(a_0*a_0 + b*b) < (side * 0.1);
       return ( circ_1 || circ_2 ) ? 15.f: 10.0f;
    };

    virtual float getBoundaryPos(BoundaryEdge edge) {
        if (edge==BND_LEFT || edge==BND_BOTTOM)
           return 0.0f;
        else
           return side;
     };
};

// 0
class SWE_SWE_Supercritical_Right_Scenario : public SWE_SWE_Supercritical_Left_Scenario{
  public:

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) {
        if (edge == BND_LEFT) {
            return INFLOW_COUPLE;
        }else{
            return OUTFLOW;
        }
      };

    virtual float maxTimeStepWidth() {return 0.125;};

    virtual float getWaterHeight(float x, float y, float offsetX, float offsetY) {
      return 10.0f;
    };

    virtual float getBoundaryPos(BoundaryEdge edge) {
       if (edge==BND_LEFT || edge==BND_TOP)
          return side;
       else if(edge==BND_BOTTOM)
          return 0;
       else
       return side * 2;
    };
};

// 1
class SWE_SWE_Subcritical_Right_Scenario : public SWE_SWE_Supercritical_Right_Scenario{
    public:
      virtual float getWaterHeight(float x, float y, float offsetX = 0, float offsetY = 0) {
        float b = (y-((side * 0.5)+offsetY));
        bool circ_1  = false;
        float a_0 = (x-((side * 0.375) +offsetX));
        bool circ_2 = sqrt(a_0*a_0 + b*b) < (side * 0.1);
        return ( circ_1 || circ_2 ) ? 20.f: 10.0f;
      }

      virtual float maxTimeStepWidth() {return 0.125;};

};

//2
class SWE_OF_Supercritical_Scenario : public SWE_RadialDamBreakScenario{//(2)
  public:
    virtual float getWaterHeight(float x, float y, float offsetX, float offsetY) {
        float a = x - (side * 0.5) + offsetX;
        float b = y - (side * 0.5) + offsetY;
        bool circ = sqrt(a * a + b * b) < (side * 0.1);
       return circ ? 20.f: 5.0f;
    };

    virtual float maxTimeStepWidth() {return 0.001;}; //TODO

    virtual float endSimulation() { return 5.f; };

    virtual int numberOfCheckpoints(){return 100; }; // to see double splash

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) {return OUTFLOW;};

};

//3
class SWE_OF_Subcritical_Scenario : public SWE_RadialDamBreakScenario {//(3)
  public: //TODO

    virtual float getWaterHeight(float x, float y, float offsetX, float offsetY) {
         // return ( x < 0.5  ? 15.f: 10.0f);
         return 5.f;
    };

    virtual float getVeloc_u(float x, float y){
          // return ( x < 0.5  ? 1.f: 0.0f);
          return 0.f;
    };

    virtual float getVeloc_v(float x, float y){
          // return ( x < 0.5  ? 0.1f: 0.0f);
          return 0.f;
    };

    virtual float endSimulation() { return 10.f; };

    virtual int numberOfCheckpoints(){return 100; }; // to see second wave (TODO second drop setFieldsValue in OF)

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) {
        if (edge == BND_RIGHT) {
            return PASSIVE;
        }
        else if(edge == BND_LEFT){
            return INFLOW;
        }
        else{
            return WALL;
        }
        assert(false && "wrong boundary type in SWE_OF_Subcritical_Scenario");
    };
};

//4
class OF_SWE_Supercritical_Scenario : public SWE_RadialDamBreakScenario { //(4)
  public:
    virtual float getWaterHeight(float x, float y, float offsetX, float offsetY) {
     return 5.0f;
    };

    virtual float getBoundaryPos(BoundaryEdge edge) {
       if (edge==BND_LEFT || edge==BND_TOP)
          return side;
       else if(edge==BND_BOTTOM)
          return 0;
       else
       return side * 2;
    };

    virtual float endSimulation() { return 8.f; };

    virtual float maxTimeStepWidth() {return 0.001;};

    virtual int numberOfCheckpoints(){return 80; };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) {
        if (edge == BND_LEFT) {
            return PASSIVE;
        }else{
            return OUTFLOW;
        }
        assert(false && "wrong boundary type in OF_SWE_Supercritical_Scenario");
    };

};

//5
class OF_SWE_Subcritical_Scenario : public SWE_RadialDamBreakScenario { //(5)
  public:
    float sideX = 7.f;
    float sideY = 5.f; //z actually

    virtual float getWaterHeight(float x, float y, float offsetX, float offsetY) {
        // float a = x - ((side * 0.5) + offsetX);
        // float b = y - ((side * 0.5) + offsetY);
        // bool circ = sqrt(a * a + b * b) < (2.f);
        // return (x > side+4 && x < side + 6) ? 15.f: 5.0f;
        return 2.0f;
    };

    virtual float getVeloc_u(float x, float y, float offsetX, float offsetY){
        // float side = 10.0f;
        // float a = x - ((side * 0.5) + offsetX);
        // float b = y - ((side * 0.5) + offsetY);
        // bool circ = sqrt(a * a + b * b) < (2.f);
        // if (circ) {
        //     return 1*(a/sqrt(a*a + b*b));
        // } else{
        //     return 0;
        // }
        return 0.f;
    };

    virtual float getVeloc_v(float x, float y, float offsetX, float offsetY){
        // float side = 10.0f;
        // float a = x - ((side * 0.5) + offsetX);
        // float b = y - ((side * 0.5) + offsetY);
        // bool circ = sqrt(a * a + b * b) < (2.f);
        // if (circ) {
        //     return 1*(b/sqrt(a*a + b*b));
        // } else{
        //     return 0;
        // }
        return 0.f;
    };

    virtual float getBoundaryPos(BoundaryEdge edge) {
      if (edge==BND_LEFT){
          return sideX;
      }
      else if (edge == BND_RIGHT) {
          return sideX*2;
      }
      else if(edge==BND_BOTTOM){
          return 0;
      }
       else{
       return sideY;
        }
    };

    virtual float endSimulation() { return 5.f; };

    virtual int numberOfCheckpoints(){return 100; };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) {
        if (edge == BND_LEFT) {
            return PASSIVE;
        }else{
            return WALL;
        }
        assert(false && "wrong boundary type in OF_SWE_Subcritical_Scenario");
    };

};


/**
 * Scenario "Bathymetry Dam Break":
 * uniform water depth, but elevated bathymetry in the centre of the domain
 */
class SWE_BathymetryDamBreakScenario : public SWE_Scenario {

  public:

    float getBathymetry(float x, float y) {
       return ( std::sqrt( (x-500.f)*(x-500.f) + (y-500.f)*(y-500.f) ) < 50.f ) ? -255.f: -260.f;
    };

	virtual float endSimulation() { return (float) 15; };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return OUTFLOW; };

    /** Get the boundary positions
     *
     * @param i_edge which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(BoundaryEdge i_edge) {
       if ( i_edge == BND_LEFT )
         return (float)0;
       else if ( i_edge == BND_RIGHT)
         return (float)1000;
       else if ( i_edge == BND_BOTTOM )
         return (float)0;
       else
         return (float)1000;
    };

    /**
     * Get the water height at a specific location.
     *
     * @param i_positionX position relative to the origin of the bathymetry grid in x-direction
     * @param i_positionY position relative to the origin of the bathymetry grid in y-direction
     * @return water height (before the initial displacement)
     */
    float getWaterHeight( float i_positionX,
                          float i_positionY, float offsetX =0, float offsetY = 0 ) {
      return (float) 260;
    }
};

/**
 * Scenario "Sea at Rest":
 * flat water surface ("sea at rest"),
 * but non-uniform bathymetry (id. to "Bathymetry Dam Break")
 * test scenario for "sea at rest"-solution
 */
class SWE_SeaAtRestScenario : public SWE_Scenario {

  public:

    float getWaterHeight(float x, float y, float offsetX = 0, float offsetY = 0) {
       return ( sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ) < 0.1f ) ? 9.9f: 10.0f;
    };
    float getBathymetry(float x, float y) {
       return ( sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ) < 0.1f ) ? 0.1f: 0.0f;
    };

};

/**
 * Scenario "Splashing Pool":
 * intial water surface has a fixed slope (diagonal to x,y)
 */
class SWE_SplashingPoolScenario : public SWE_Scenario {

  public:

    float getBathymetry(float x, float y) {
       return -250.f;
    };

    float getWaterHeight(float x, float y, float offsetX = 0, float offsetY = 0) {
    	return 250.0f+(5.0f-(x+y)/200);
    };

	virtual float endSimulation() { return (float) 15; };

    /** Get the boundary positions
     *
     * @param i_edge which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(BoundaryEdge i_edge) {
       if ( i_edge == BND_LEFT )
         return (float)0;
       else if ( i_edge == BND_RIGHT)
         return (float)1000;
       else if ( i_edge == BND_BOTTOM )
         return (float)0;
       else
         return (float)1000;
    };

};

/**
 * Scenario "Splashing Cone":
 * bathymetry forms a circular cone
 * intial water surface designed to form "sea at rest"
 * but: elevated water region in the centre (similar to radial dam break)
 */
class SWE_SplashingConeScenario : public SWE_Scenario {

  public:

    float getWaterHeight(float x, float y, float offsetX = 0, float offsetY = 0) {
       float r = sqrt( (x-0.5f)*(x-0.5f) + (y-0.5f)*(y-0.5f) );
       float h = 4.0f-4.5f*(r/0.5f);

       if (r<0.1f) h = h+1.0f;

       return (h>0.0f) ? h : 0.0f;
    };

    float getBathymetry(float x, float y) {
       float r = sqrt( (x-0.5f)*(x-0.5f) + (y-0.5f)*(y-0.5f) );
       return 1.0f+9.0f*( (r < 0.5f) ? r : 0.5f);
    };

    float waterHeightAtRest() { return 4.0f; };
    float endSimulation() { return 0.5f; };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return OUTFLOW; };
};

#endif
