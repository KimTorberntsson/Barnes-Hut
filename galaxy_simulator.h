/*
 * File: galaxy_simulator.h
 * --------------
 * Header file for barnes_hutgalaxy_simulator.c
 * Written by:
 * Joel Backsell, Kim Torberntsson & Alexander Bilock.
 */

#ifndef Barnes_Hut_galaxy_simulator_h
#define Barnes_Hut_galaxy_simulator_h

//Function declarations
double frand(double xmin, double xmax);
void printtime(clock_t s, clock_t e);
void display(void);
void time_step(void);
void bounce(double *x, double *y, double *u, double *v);
void update_forces(double *x, double *y, double *forceX, double *forceY);

#endif
