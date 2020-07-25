/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

struct spring {
  /* index pair is used in two ways: 
     1. for elastic springs, stores *INDECIES* of two end mass points
     2. for collision springs, stores *collision point on cube*, 
        and *point of contact on collision object*   */
  pair<point, point> indexpair;
  double rest_length;
};

void createStructuralSprings(vector<spring> & structuralSprings);
void createShearSprings(vector<spring> & shearSprings);
void createBendSprings(vector<spring> & bendSprings);

point computeHookForce(spring s, double kHook,struct world * jello);
point computeDampingForce(spring s, double kDamping, point vA, point vB, struct world * jello);

void computeAcceleration(struct world * jello, struct point a[8][8][8], vector<spring> * structuralSprings, 
                                                                        vector<spring> * shearSprings, 
                                                                        vector<spring> * bendSprings);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello, vector<spring> * structuralSprings, 
                                 vector<spring> * shearSprings, 
                                 vector<spring> * bendSprings);
void RK4(struct world * jello, vector<spring> * structuralSprings, 
                               vector<spring> * shearSprings, 
                               vector<spring> * bendSprings);

#endif

