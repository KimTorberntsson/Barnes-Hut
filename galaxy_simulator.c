/*
 * File: galaxy_simulator.c
 * --------------
 * Implements the bruto force algorithm for n-body
 * simulation with galaxy-like initial conditions.
 * Uses OpenGl and GLUT for graphics. Written by:
 * Joel Backsell, Kim Torberntsson & Alexander Bilock.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <GLUT/glut.h>

#include "galaxy_simulator.h"
#include "graphics.h"

//Some constants and global variables
int N = 2500;
const double L = 1, W = 1, dt = 1e-3, alpha = 0.25, V = 50, epsilon = 1e-1, grav = 0.04; //grav should be 100/N
double *x, *y, *u, *v, *force_x, *force_y, *mass;
struct node_t *root;

/*
 * Function for producing a random number between two double values
 */
double frand(double xmin, double xmax) {
    return xmin + (xmax-xmin)*rand()/RAND_MAX;
}

/*
 * Prints the time between two clocks in seconds
 */
void printtime(clock_t s, clock_t e)
{
    printf("Time: %f seconds\n", (double)(e-s)/CLOCKS_PER_SEC);
}


/*
 * This function is called every time GLUT refreshes the display.
 */
void display(void) {
    //Do a time step
    time_step();
    
    //Draw points
    drawPoints(x, y, N);
}

/*
 * Updates the positions of the particles of a time step.
 */
void time_step(void) {
    //Update forces
    update_forces(x, y, force_x, force_y);
    
    //Update velocities and positions
    for(int i = 0; i < N; i++) {
        double ax = force_x[i]/mass[i];
        double ay = force_y[i]/mass[i];
        u[i] += ax*dt;
        v[i] += ay*dt;
        x[i] += u[i]*dt;
        y[i] += v[i]*dt;
        bounce(&x[i], &y[i], &u[i], &v[i]);
    }
}

/*
 * If a particle moves beyond any of the boundaries then bounce it back
 */
void bounce(double *x, double *y, double *u, double *v) {
    double W = 1.0f, H = 1.0f;
    if(*x > W) {
        *x = 2*W - *x;
        *u = -*u;
    }
    
    if(*x < 0) {
        *x = -*x;
        *u = -*u;
    }
    
    if(*y > H) {
        *y = 2*H - *y;
        *v = -*v;
    }
    
    if(*y < 0) {
        *y = -*y;
        *v = -*v;
    }
}

/*
 * Function for updating the forces of all the stars
 */
void update_forces(double *x, double *y, double *force_x, double *force_y) {
    //Set the forces to zero
    for(int i = 0; i < N; i++) {
        force_x[i] = 0;
        force_y[i] = 0;
    }
    //Update the forces
    for(int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double r = sqrt((x[i] - x[j])*(x[i] - x[j]) + (y[i] - y[j])*(y[i] - y[j]));
            double temp = -grav*mass[i]*mass[j]/((r + epsilon)*(r + epsilon)*(r + epsilon));
            force_x[i] += temp*(x[i] - x[j]);
            force_y[i] += temp*(y[i] - y[j]);
            force_x[j] -= temp*(x[i] - x[j]);
            force_y[j] -= temp*(y[i] - y[j]);
        }
    }
}


/*
 * The main method
 */
int main(int argc, char *argv[]) {
    
    //The first arguments sets if the graphics should be used
    int graphics = 1;
    if (argc > 1) {
        graphics = atoi(argv[1]);
    }
    
    //The second argument sets the number of time steps
    int time_steps = 100;
    if (argc > 2) {
        time_steps = atoi(argv[2]);
    }
    
    x = (double *)malloc(N*sizeof(double));
    y = (double *)malloc(N*sizeof(double));
    u = (double *)malloc(N*sizeof(double));
    v = (double *)malloc(N*sizeof(double));
    force_x = (double *)malloc(N*sizeof(double));
    force_y = (double *)malloc(N*sizeof(double));
    mass = (double *)malloc(N*sizeof(double));
    
    
    for(int i = 0; i < N; i++) {
        mass[i] = 1;
        double R = frand(0, L/4);
        double theta = frand(0, 2*M_PI);
        x[i] = L/2 + R*cos(theta);
        y[i] = W/2 + alpha*R*sin(theta);
        double R_prim = sqrt(pow(x[i] - L/2, 2) + pow(y[i] - W/2, 2));
        y[i] = W/2 + alpha*R*sin(theta);
        u[i] = -V*R_prim*sin(theta);
        v[i] = V*R_prim*cos(theta);
        force_x[i]=0;
        force_y[i]=0;
    }
    
    /* Run the GLUT display function if the graphics mode is on.
     * Otherwise just run the simulations without graphics
     */
    if (graphics) {
        // Initialize the graphics routines
        graphicsInit(&argc, argv, display);
        
        //Used for the display window
        glutMainLoop();
    } else {
        //Begin taking time
        long start = clock();
        
        //The main loop
        for(int i = 0; i < time_steps; i++) {
            time_step();
        }
        
        //Stop taking time and print elapsed time
        long stop = clock();
        printtime(start, stop);
    }
    
    //Free memory
    free(x);
    free(y);
    free(u);
    free(v);
    free(force_x);
    free(force_y);
    free(mass);
    
    return 0;
}
