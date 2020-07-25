/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1

*/

#include "jello.h"
#include "physics.h"
#include <utility>
#include <vector>
#include <iostream>

using namespace std;

void print_point(point p) {
  cout<<"Values: "<<"("<<p.x<<", "<<p.y<<", "<<p.z<<")"<<endl;
}

void print_spring_ind_pair(spring s) {
  cout<<"Values of spring index pair: "<<"("<<s.indexpair.first.x<<", "<<s.indexpair.first.y<<", "<<s.indexpair.first.z<<"), ";
  cout<<"("<<s.indexpair.second.x<<", "<<s.indexpair.second.y<<", "<<s.indexpair.second.z<<"). ";
  cout<<"Rest length is: "<<s.rest_length<<endl;
}

// modified based on PROCESS_NEIGHBOR used to record a mass point's neighbors' undeformed position.
// points inside the cude are also considered to fully model the mass-spring system.
int i,j,k,ip,jp,kp;
point pN;
#define GET_NEIGHBOUR_PT(di,dj,dk) \
  ip=i+(di);\
  jp=j+(dj);\
  kp=k+(dk);\
  if( ! ( (ip>7) || (ip<0) ||(jp>7) || (jp<0) ||(kp>7) || (kp<0) )) \
  {\
    pN.x = 1.0 * ip / 7;\
    pN.y = 1.0 * jp / 7;\
    pN.z = 1.0 * kp / 7;\
  }else{\
    pN.x = -404; /*neighbor not found*/\
  }\

// computes 1D index in force filed of a 3D index (fi,fj,fk) given resolution of force field.
// Returns integer result in ind
#define GET_FORCEFIELD_IND(fi, fj, fk, resolution, ind) \
  ind = fi*resolution*resolution + fj*resolution + fk;

// overload == operator for spring index pair equality comparison used by std::find.
// Purpose of this is to avoid duplicate inserting of a spring (i.e. spring (a,b)
// and (b,a) are a same spring)
bool operator==(const spring &lhs, const spring rhs) {
    if (  (lhs.indexpair.first.x == rhs.indexpair.second.x
         &&lhs.indexpair.first.y == rhs.indexpair.second.y
         &&lhs.indexpair.first.z == rhs.indexpair.second.z)
        &&(lhs.indexpair.second.x == rhs.indexpair.first.x
         &&lhs.indexpair.second.y == rhs.indexpair.first.y
         &&lhs.indexpair.second.z == rhs.indexpair.first.z) ) {
        return true;
      }
    return false;
}

// check if a spring is new and unique by finding it in the vector of springs
bool isNewSpring (vector<spring> *sp_vector, spring s) {
  return !(find(sp_vector->begin(), sp_vector->end(), s) != sp_vector->end());
}

// create structural springs, called ONCE in myinit()
void createStructuralSprings(vector<spring> & structuralSprings) {
  spring currSpring;
  currSpring.rest_length = 1.0/7.0;
  point springIndA, springIndB;
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        GET_NEIGHBOUR_PT(1,0,0);
        if(pN.x!=-404) {
          pMAKE(i,j,k,springIndA);
          pMAKE(i+1,j,k, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&structuralSprings,currSpring)) structuralSprings.push_back(currSpring);
        }
        GET_NEIGHBOUR_PT(0,1,0);
        if(pN.x!=-404) {
          pMAKE(i,j,k,springIndA);
          pMAKE(i,j+1,k, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&structuralSprings,currSpring)) structuralSprings.push_back(currSpring);
        }
        GET_NEIGHBOUR_PT(0,0,1);
        if(pN.x!=-404) {
          pMAKE(i,j,k,springIndA);
          pMAKE(i,j,k+1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&structuralSprings,currSpring)) structuralSprings.push_back(currSpring);
        }
        GET_NEIGHBOUR_PT(-1,0,0);
        if(pN.x!=-404) {
          pMAKE(i,j,k,springIndA);
          pMAKE(i-1,j,k, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&structuralSprings,currSpring)) structuralSprings.push_back(currSpring);
        }
        GET_NEIGHBOUR_PT(0,-1,0);
        if(pN.x!=-404) {
          pMAKE(i,j,k,springIndA);
          pMAKE(i,j-1,k, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&structuralSprings,currSpring)) structuralSprings.push_back(currSpring);
        }
        GET_NEIGHBOUR_PT(0,0,-1);
        if(pN.x!=-404) {
          pMAKE(i,j,k,springIndA);
          pMAKE(i,j,k-1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&structuralSprings,currSpring)) structuralSprings.push_back(currSpring);
        }
      }
}

// create shear springs, called ONCE in myinit()
void createShearSprings(vector<spring> & shearSprings) {
  spring currSpring;
  point springIndA, springIndB;
  // there are only two different rest length
  double short_rest_length, long_rest_length;
  short_rest_length = 1.0/7.0*sqrt(2);
  long_rest_length = 1.0/7.0*sqrt(3);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        // make curr (i,j,k) as first component (point) of a spring
        pMAKE(i,j,k,springIndA);

        // for below 16 neighbors, they have rest length of 1/7*sqrt(2)
        GET_NEIGHBOUR_PT(1,1,0);
        if(pN.x!=-404) {
          pMAKE(i+1,j+1,k, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = short_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(-1,1,0);
        if(pN.x!=-404) {
          pMAKE(i-1,j+1,k, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = short_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(-1,-1,0);
        if(pN.x!=-404) {
          pMAKE(i-1,j-1,k, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = short_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(1,-1,0);
        if(pN.x!=-404) {
          pMAKE(i+1,j-1,k, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = short_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        /***** divider *******/
        GET_NEIGHBOUR_PT(0,1,1);
        if(pN.x!=-404) {
          pMAKE(i,j+1,k+1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = short_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(0,-1,1);
        if(pN.x!=-404) {
          pMAKE(i,j-1,k+1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = short_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(0,-1,-1);
        if(pN.x!=-404) {
          pMAKE(i,j-1,k-1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = short_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(0,1,-1);
        if(pN.x!=-404) {
          pMAKE(i,j+1,k-1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = short_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        /***** divider *******/
        GET_NEIGHBOUR_PT(1,0,1);
        if(pN.x!=-404) {
          pMAKE(i+1,j,k+1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = short_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(-1,0,1);
        if(pN.x!=-404) {
          pMAKE(i-1,j,k+1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = short_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(-1,0,-1);
        if(pN.x!=-404) {
          pMAKE(i-1,j,k-1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = short_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(1,0,-1);
        if(pN.x!=-404) {
          pMAKE(i+1,j,k-1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = short_rest_length;
            shearSprings.push_back(currSpring);
          }
        }

        // for below 8 neighbors, they have length of 1/7*sqrt(3)
        /***** divider *******/
        GET_NEIGHBOUR_PT(1,1,1);
        if(pN.x!=-404) {
          pMAKE(i+1,j+1,k+1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = long_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(-1,1,1);
        if(pN.x!=-404) {
          pMAKE(i-1,j+1,k+1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = long_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(-1,-1,1);
        if(pN.x!=-404) {
          pMAKE(i-1,j-1,k+1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = long_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(1,-1,1);
        if(pN.x!=-404) {
          pMAKE(i+1,j-1,k+1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = long_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        /***** divider *******/
        GET_NEIGHBOUR_PT(1,1,-1);
        if(pN.x!=-404) {
          pMAKE(i+1,j+1,k-1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = long_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(-1,1,-1);
        if(pN.x!=-404) {
          pMAKE(i-1,j+1,k-1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = long_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(-1,-1,-1);
        if(pN.x!=-404) {
          pMAKE(i-1,j-1,k-1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = long_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        GET_NEIGHBOUR_PT(1,-1,-1);
        if(pN.x!=-404) {
          pMAKE(i+1,j-1,k-1, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&shearSprings,currSpring)) {
            currSpring.rest_length = long_rest_length;
            shearSprings.push_back(currSpring);
          }
        }
        /***** end *******/
      }
}

// create bend springs, called ONCE in myinit()
void createBendSprings(vector<spring> & bendSprings) {
  spring currSpring;
  currSpring.rest_length = 1.0/7.0 * 2;
  point springIndA, springIndB;
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        GET_NEIGHBOUR_PT(2,0,0);
        if(pN.x!=-404) {
          pMAKE(i,j,k,springIndA);
          pMAKE(i+2,j,k, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&bendSprings,currSpring)) bendSprings.push_back(currSpring);
        }
        GET_NEIGHBOUR_PT(0,2,0);
        if(pN.x!=-404) {
          pMAKE(i,j,k,springIndA);
          pMAKE(i,j+2,k, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&bendSprings,currSpring)) bendSprings.push_back(currSpring);
        }
        GET_NEIGHBOUR_PT(0,0,2);
        if(pN.x!=-404) {
          pMAKE(i,j,k,springIndA);
          pMAKE(i,j,k+2, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&bendSprings,currSpring)) bendSprings.push_back(currSpring);
        }
        GET_NEIGHBOUR_PT(-2,0,0);
        if(pN.x!=-404) {
          pMAKE(i,j,k,springIndA);
          pMAKE(i-2,j,k, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&bendSprings,currSpring)) bendSprings.push_back(currSpring);
        }
        GET_NEIGHBOUR_PT(0,-2,0);
        if(pN.x!=-404) {
          pMAKE(i,j,k,springIndA);
          pMAKE(i,j-2,k, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&bendSprings,currSpring)) bendSprings.push_back(currSpring);
        }
        GET_NEIGHBOUR_PT(0,0,-2);
        if(pN.x!=-404) {
          pMAKE(i,j,k,springIndA);
          pMAKE(i,j,k-2, springIndB);
          currSpring.indexpair = make_pair(springIndA, springIndB);
          if(isNewSpring(&bendSprings,currSpring)) bendSprings.push_back(currSpring);
        }
      }
}

// computes Hook's force exerted on the first component of a spring
point computeHookElastic(spring s, double kHook, struct world * jello) {
  point normL, resultForce;
  double length, scaler, test_length;
  int ptAi, ptAj, ptAk,
      ptBi, ptBj, ptBk;
  // although indexpairs are pair<point, point>, where each point have
  // double x,y,z, they are initialized to be integer for all non-collision
  // springs. Hence using (int) here does not having any rounding.
  ptAi = (int) s.indexpair.first.x;
  ptAj = (int) s.indexpair.first.y;
  ptAk = (int) s.indexpair.first.z;
  ptBi = (int) s.indexpair.second.x;
  ptBj = (int) s.indexpair.second.y;
  ptBk = (int) s.indexpair.second.z; 

  // calculate L based on the actual position values of the two mass 
  // points connected by spring s and normalize L
  pDIFFERENCE(jello->p[ptAi][ptAj][ptAk], jello->p[ptBi][ptBj][ptBk], normL);
  pNORMALIZE(normL);
  scaler = - kHook * (length - s.rest_length);
  pMULTIPLY(normL, scaler, resultForce);
  // Fhook = -kHook * (|L|-s.rest_length) * normalized L
  return resultForce;
}

// computes Hook's force for a collision spring
point computeHookCollision(spring s, double kHook) {
  point normL, resultForce;
  double length, scaler;
  // calculate L based on the actual position values of the two mass 
  // points connected by spring s and normalize L
  pDIFFERENCE(s.indexpair.first, s.indexpair.second, normL);
  pNORMALIZE(normL);
  scaler = - kHook * (length-s.rest_length);
  pMULTIPLY(normL, scaler, resultForce);
  // Fhook = -kHook * (|L|-s.rest_length) * normalized L
  return resultForce;
}

// computes damping force exerted on the first component of a spring
point computeDampElastic(spring s, double kDamping, point vA, point vB, struct world * jello) {
  point normL, vAminusvB, resultForce;
  double length, scaler, dotProd;
  // calculate L
  pDIFFERENCE(jello->p[(int) s.indexpair.first.x][(int) s.indexpair.first.y][(int) s.indexpair.first.z],
              jello->p[(int) s.indexpair.second.x][(int) s.indexpair.second.y][(int) s.indexpair.second.z],
              normL);
  // calculate vA - vB
  pDIFFERENCE(vA, vB, vAminusvB);
  // calculate (vA-vB) dot L
  dotProd = vAminusvB.x * normL.x + vAminusvB.y * normL.y + vAminusvB.z * normL.z;
  // normalize L
  pNORMALIZE(normL);
  scaler = - kDamping * dotProd / length;
  // calculate resultForce
  pMULTIPLY(normL, scaler, resultForce);
  // Fdamping = -kDamping * ((vA-vB) dot L) / |L| * normalized L
  return resultForce;
}

// computes damping force for a collision spring
point computeDampCollision(spring s, double kDamping, point vA) {
  point vB, normL, vAminusvB, resultForce;
  double length, scaler, dotProd;
  // vB is collision pt's velocity, which is 0
  pMAKE(0,0,0, vB);
  // calculate L
  pDIFFERENCE(s.indexpair.first, s.indexpair.second, normL);
  // calculate vA - vB
  pDIFFERENCE(vA, vB, vAminusvB);
  // calculate (vA-vB) dot L
  dotProd = vAminusvB.x * normL.x + vAminusvB.y * normL.y + vAminusvB.z * normL.z;
  // normalize L
  pNORMALIZE(normL);
  scaler = - kDamping * dotProd / length;
  // calculate resultForce
  pMULTIPLY(normL, scaler, resultForce);
  // Fdamping = -kDamping * ((vA-vB) dot L) / |L| * normalized L
  return resultForce;
}

// calculate forces on mass points given a vector of strings
// updates force array 'F'
void iterateSprings(vector<spring> * vec_of_springs, struct point F[8][8][8], struct world * jello) {
  point FhookOnA, FdampingOnA, FfinalOnA, FhookOnB, FdampingOnB, FfinalOnB;
  int pt_A_ind_x, pt_A_ind_y, pt_A_ind_z,
      pt_B_ind_x, pt_B_ind_y, pt_B_ind_z;
  // for each spring in a vector of springs
  for (int i=0; i<vec_of_springs->size(); i++) {
    // get indecies of two mass point of current spring
    pt_A_ind_x = (int) vec_of_springs->at(i).indexpair.first.x;
    pt_A_ind_y = (int) vec_of_springs->at(i).indexpair.first.y;
    pt_A_ind_z = (int) vec_of_springs->at(i).indexpair.first.z;
    pt_B_ind_x = (int) vec_of_springs->at(i).indexpair.second.x;
    pt_B_ind_y = (int) vec_of_springs->at(i).indexpair.second.y;
    pt_B_ind_z = (int) vec_of_springs->at(i).indexpair.second.z;

    // compute Hook and damping forces exerted on the first point of the spring
    FhookOnA = computeHookElastic(vec_of_springs->at(i), jello->kElastic, jello);

    FdampingOnA = computeDampElastic(vec_of_springs->at(i), jello->dElastic, 
                                      jello->v[pt_A_ind_x][pt_A_ind_y][pt_A_ind_z],
                                      jello->v[pt_B_ind_x][pt_B_ind_y][pt_B_ind_z], jello);

    // sum up Hook and damping forces exerted on the first point of the spring
    pSUM(FhookOnA, FdampingOnA, FfinalOnA);
    // negate Hook and damping forces, they are now forces exerted on the 
    // second point of the spring by Newton's 3rd law
    pNEG(FhookOnA, FhookOnB);
    pNEG(FdampingOnA, FdampingOnB);
    pSUM(FhookOnB, FdampingOnB, FfinalOnB);

    // update the two points' force in F
    F[pt_A_ind_x][pt_A_ind_y][pt_A_ind_z].x += FfinalOnA.x;
    F[pt_A_ind_x][pt_A_ind_y][pt_A_ind_z].y += FfinalOnA.y;
    F[pt_A_ind_x][pt_A_ind_y][pt_A_ind_z].z += FfinalOnA.z;

    F[pt_B_ind_x][pt_B_ind_y][pt_B_ind_z].x += FfinalOnB.x;
    F[pt_B_ind_x][pt_B_ind_y][pt_B_ind_z].y += FfinalOnB.y;
    F[pt_B_ind_x][pt_B_ind_y][pt_B_ind_z].z += FfinalOnB.z;
  }
}

// calculate penalty force when collision happens
// updates force array 'F'
void collisionResponse (spring collisionSpring, struct point F[8][8][8], point ptOfContact, struct world * jello) {
  point collisionFinalForce;
  // calculate collision Hook and damping force and sum them up
  pSUM(computeHookCollision(collisionSpring, jello->kCollision), 
       computeDampCollision(collisionSpring, jello->dCollision, jello->v[i][j][k]), 
       collisionFinalForce);
  // update F
  F[i][j][k].x += collisionFinalForce.x;
  F[i][j][k].y += collisionFinalForce.y;
  F[i][j][k].z += collisionFinalForce.z;
}

// detects collision and handles them by calling collisionResponse()
// Note: the indexpair for a collision spring actually stores positions of
// a collision point and point of collision on a wall
void detectWallCollision (struct point F[8][8][8], struct world * jello) {
  spring collisionSpring;
  collisionSpring.rest_length = 0;
  point ptOfContact, collisionHookForce, collisionDampForce, collisionFinalForce;
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        // collision with +x wall
        if (jello->p[i][j][k].x > 2) {
          pMAKE(2, jello->p[i][j][k].y, jello->p[i][j][k].z, ptOfContact);
          collisionSpring.indexpair = make_pair(jello->p[i][j][k], ptOfContact);
          collisionResponse(collisionSpring, F, ptOfContact, jello);
        }
        // collision with -x wall
        if (jello->p[i][j][k].x < -2) {
          pMAKE(-2, jello->p[i][j][k].y, jello->p[i][j][k].z, ptOfContact);
          collisionSpring.indexpair = make_pair(jello->p[i][j][k], ptOfContact);
          collisionResponse(collisionSpring, F, ptOfContact, jello);
        }
        // collision with +y wall
        if (jello->p[i][j][k].y > 2) {
          pMAKE(jello->p[i][j][k].x, 2, jello->p[i][j][k].z, ptOfContact);
          collisionSpring.indexpair = make_pair(jello->p[i][j][k], ptOfContact);
          collisionResponse(collisionSpring, F, ptOfContact, jello);
        }
        // collision with -y wall
        if (jello->p[i][j][k].y < -2) {
          pMAKE(jello->p[i][j][k].x, -2, jello->p[i][j][k].z, ptOfContact);
          collisionSpring.indexpair = make_pair(jello->p[i][j][k], ptOfContact);
          collisionResponse(collisionSpring, F, ptOfContact, jello);
        }
        // collision with +z wall
        if (jello->p[i][j][k].z > 2) {
          pMAKE(jello->p[i][j][k].x, jello->p[i][j][k].y, 2, ptOfContact);
          collisionSpring.indexpair = make_pair(jello->p[i][j][k], ptOfContact);
          collisionResponse(collisionSpring, F, ptOfContact, jello);
        }
        // collision with -x wall
        if (jello->p[i][j][k].z < -2) {
          pMAKE(jello->p[i][j][k].x, jello->p[i][j][k].y, -2, ptOfContact);
          collisionSpring.indexpair = make_pair(jello->p[i][j][k], ptOfContact);
          collisionResponse(collisionSpring, F, ptOfContact, jello);
        }
      }
}

void detectPlaneCollision(struct point F[8][8][8], struct world *jello) {
  spring collisionSpring;
  collisionSpring.rest_length = 0;
  // plane normal with correct sign has been calculated in jello.cpp's myinit()
  point n = inclinedPlaneNormal.first;

  point ptOfContact;
  double currValueOfEq;

  // calculate ptOfContact is to find orthogonal projection
  // of a point p on jello cube onto plane defined by a unit normal [a,b,c]
  // hence: orthog_proj(p) = p - ((p-v) dot n)*n
  // reference: https://math.stackexchange.com/questions/444968/project-a-point-in-3d-on-a-given-plane
  // v is just a point on plane used for calculating orthogonal projection of a point on cube onto 
  // collision plane
  point v, PminusV;
  double PminusVdotN;
  pMAKE(0,0,-jello->d/jello->c, v);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        currValueOfEq = jello->a * jello->p[i][j][k].x + 
                        jello->b * jello->p[i][j][k].y + 
                        jello->c * jello->p[i][j][k].z + jello->d;
        // if negSide flag is true, cube is on the negative side of the plane
        // so negate value of equation to make <0 comparision hold
        if (inclinedPlaneNormal.second == true) currValueOfEq=-currValueOfEq;
        if(currValueOfEq < 0)
        {
          // below is the implementation of the reference link
          pDIFFERENCE(jello->p[i][j][k], v, PminusV);
          PminusVdotN = PminusV.x * n.x + PminusV.y * n.y + PminusV.z * n.z;
          pMULTIPLY(n, PminusVdotN, ptOfContact);
          pDIFFERENCE(jello->p[i][j][k], ptOfContact, ptOfContact);

          // make the collision spring and response
          collisionSpring.indexpair = make_pair(jello->p[i][j][k], ptOfContact);
          collisionResponse(collisionSpring, F, ptOfContact, jello);
        }
      }
}

int checkCubePtIndexAndGetFFArrInd (double pos_value, double integralGridLength) {
  // check if a point's component is out of negative world boundry.
  // If yes, return correspond force field array index as 0.
  if (pos_value < -2) return 0;
  // check if a point's component is out of positive world boundry.
  // If yes, return correspond force field array index as furtherest possible positive boundry
  if (pos_value >  2) return (int) 4 / integralGridLength; // i.e. res - 1
  // now the point's component is within the world boundary.
  // Proof: Consider p.x. Known p.x >= (-2 + i * 4 / (res - 1), consider equality.
  // -> (p.x+2) / (4/res-1) = i -> (p.x+2)/integralGridLength.
  // Take floor() to only consider the "minimum" point on all x,y,z direction in a cube
  return (int) floor((pos_value + 2) / integralGridLength);
}

void interpolateForceFiled(struct point F[8][8][8], struct world * jello) {
  double res = jello->resolution;
  // length of a integral grid, which is h
  double integraGridLength = 4 / (res - 1);
  // indecies
  int fi,fj,fk,           // for force field array (1D array)
      alpha, beta, gamma; // for the point of interpolation (3D array)
  // for storing value get from GET_FORCEFIELD_IND macro
  int ind;
  // eight corners of the integral grid;
  point F000, F001, F010, F011, F100, F101, F110, F111;
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        // compute nearest force field array index for curr point
        fi = checkCubePtIndexAndGetFFArrInd(jello->p[i][j][k].x, integraGridLength);
        fj = checkCubePtIndexAndGetFFArrInd(jello->p[i][j][k].y, integraGridLength);
        fk = checkCubePtIndexAndGetFFArrInd(jello->p[i][j][k].z, integraGridLength);

        // compute alpha, beta, gamma coefficients for interpolating curr point's 
        // force F(alpha, beta, gamma)
        alpha = ( jello->p[i][j][k].x - (-2 + fi * integraGridLength) ) / integraGridLength;
        beta  = ( jello->p[i][j][k].y - (-2 + fj * integraGridLength) ) / integraGridLength;
        gamma = ( jello->p[i][j][k].z - (-2 + fk * integraGridLength) ) / integraGridLength;

        // get eight surrounding force vectors around curr point
        GET_FORCEFIELD_IND(fi  ,fj  ,fk  ,res, ind);
        F000 = jello->forceField[ind];

        GET_FORCEFIELD_IND(fi  ,fj  ,fk+1,res, ind);
        F001 = jello->forceField[ind];

        GET_FORCEFIELD_IND(fi  ,fj+1,fk  ,res, ind);
        F010 = jello->forceField[ind];

        GET_FORCEFIELD_IND(fi  ,fj+1,fk+1,res, ind);
        F011 = jello->forceField[ind];

        GET_FORCEFIELD_IND(fi+1,fj  ,fk  ,res, ind);
        F100 = jello->forceField[ind];

        GET_FORCEFIELD_IND(fi+1,fj  ,fk+1,res, ind);
        F101 = jello->forceField[ind];

        GET_FORCEFIELD_IND(fi+1,fj+1,fk  ,res, ind);
        F110 = jello->forceField[ind];

        GET_FORCEFIELD_IND(fi+1,fj+1,fk+1,res, ind);
        F111 = jello->forceField[ind];

        // multiply corresponding coefficients to each of the eight surrounding force vectors
        pMULTIPLY(F000, (1 - alpha) * (1 - beta) * (1 - gamma), F000);
        pMULTIPLY(F001, (1 - alpha) * (1 - beta) *      gamma , F001);
        pMULTIPLY(F010, (1 - alpha) *      beta  * (1 - gamma), F010);
        pMULTIPLY(F011, (1 - alpha) *      beta  *      gamma , F011);
        pMULTIPLY(F100,      alpha  * (1 - beta) * (1 - gamma), F100);
        pMULTIPLY(F101,      alpha  * (1 - beta) *      gamma , F101);
        pMULTIPLY(F110,      alpha  *      beta  * (1 - gamma), F110);
        pMULTIPLY(F111,      alpha  *      beta  *      gamma , F111);

        // do the trilinear interpolation to update force on curr point
        F[i][j][k].x += (F000.x + F001.x + F010.x + F011.x + F100.x + F101.x + F110.x + F111.x);
        F[i][j][k].y += (F000.y + F001.y + F010.y + F011.y + F100.y + F101.y + F110.y + F111.y);
        F[i][j][k].z += (F000.z + F001.z + F010.z + F011.z + F100.z + F101.z + F110.z + F111.z);
      }
}

// apply equal forces on all points of the jello given user drag
// direction and magnitude of the user drag has been determined in input.cpp
// here we only need to normalize the defined direction and multiply calculated
// magnitude to the normalize direction. 
void handleUserDrag(struct point F[8][8][8], struct world * jello) {
  double length;
  point currDragForce = dragForceDirection;
  pNORMALIZE(currDragForce);
  // magnify drag force for worlds with <0.001 timestep
  dragForceMagnitude = dragForceMagnitude * (1.0/jello->dt)/1000;
  pMULTIPLY(currDragForce, dragForceMagnitude, currDragForce);
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
	    for (k=0; k<=7; k++)
      {
        F[i][j][k].x += currDragForce.x;
        F[i][j][k].y += currDragForce.y;
        F[i][j][k].z += currDragForce.z;
      }

  // clear magnitude for future time steps
  dragForceMagnitude = 0;
}

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8], vector<spring> * structuralSprings, 
                                                                        vector<spring> * shearSprings, 
                                                                        vector<spring> * bendSprings)
{
  point F[8][8][8];
  // initialize F to 0 to avoid garbage value (since forces cumulates)
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
	    for (k=0; k<=7; k++)
      {
        F[i][j][k].x = 0.0;
        F[i][j][k].y = 0.0;
        F[i][j][k].z = 0.0;
      }

  // iterate all springs and update elastic and damping forces on mass points
  iterateSprings(structuralSprings,F,jello);
  iterateSprings(shearSprings,F,jello);
  iterateSprings(bendSprings,F,jello);

  // detect collisions and respond collisions with penalty method
  if (jello->incPlanePresent==1) detectPlaneCollision(F, jello);
  detectWallCollision(F, jello);

  // if a force field is present, update forces on mass point with trilinear interpolation
  if (jello->resolution != 0) interpolateForceFiled(F, jello);

  // if interactive mode is turned on by pressing 'i',
  // and if user draged the left mouse, apply equal forces
  // on all points of the jello
  if(dragCase!=-1 && dragForceMagnitude>0) handleUserDrag(F, jello);

  // at this stage, force matrix is finalized. Compute acceleration for each mass point.
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
	    for (k=0; k<=7; k++)
      {
        a[i][j][k].x = 1.0 / jello->mass * F[i][j][k].x;
        a[i][j][k].y = 1.0 / jello->mass * F[i][j][k].y;
        a[i][j][k].z = 1.0 / jello->mass * F[i][j][k].z;
      }
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello, vector<spring> * structuralSprings, 
                                 vector<spring> * shearSprings, 
                                 vector<spring> * bendSprings)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a, structuralSprings, shearSprings, bendSprings);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;
      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello, vector<spring> * structuralSprings, 
                               vector<spring> * shearSprings, 
                               vector<spring> * bendSprings)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(&buffer, a, structuralSprings, shearSprings, bendSprings);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a, structuralSprings, shearSprings, bendSprings);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a, structuralSprings, shearSprings, bendSprings);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a, structuralSprings, shearSprings, bendSprings);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }
  return;  
}
