/*
 * File: barnes_hut.c
 * --------------
 * Implements the Barnes Hut algorithm for n-body
 * simulation with galaxy-like initial conditions. 
 * Uses OpenGl and GLUT for graphics. Written by: 
 * Joel Backsell, Kim Torberntsson & Alexander Bilock.
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <GLUT/glut.h>

#include "barnes_hut.h"
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
void print_time(clock_t s, clock_t e)
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
	//Allocate memory for root
	root = malloc(sizeof(struct node_t));
	set_node(root);
	root->min_x = 0; root->max_x = 1; root->min_y = 0; root->max_y = 1;
	
	//Put particles in tree
	for(int i = 0; i < N; i++) {
		put_particle_in_tree(i, root);
	}
	
	//Calculate mass and center of mass
	calculate_mass(root);
	calculate_center_of_mass_x(root);
	calculate_center_of_mass_y(root);
	
	//Calculate forces
	update_forces();

	//Update velocities and positions
	for(int i = 0; i < N; i++) {
		double ax = force_x[i]/mass[i];
		double ay = force_y[i]/mass[i];
		u[i] += ax*dt;
		v[i] += ay*dt;
		x[i] += u[i]*dt;
		y[i] += v[i]*dt;

		/* This of course doesn't make any sense physically, 
		 * but makes sure that the particles stay within the
		 * bounds. Normally the particles won't leave the
		 * area anyway.
		*/
		bounce(&x[i], &y[i], &u[i], &v[i]);
	}

	//Free memory
	free_node(root);
	free(root);
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
 * Puts a particle recursively in the Barnes Hut quad-tree.
 */
void put_particle_in_tree(int new_particle, struct node_t *node) {
	//If no particle is assigned to the node
	if(!node->has_particle) {
		node->particle = new_particle;
		node->has_particle = 1;
	} 
	//If the node has no children
	else if(!node->has_children) {
		//Allocate and initiate children
		node->children = malloc(4*sizeof(struct node_t));
		for(int i = 0; i < 4; i++) {
			set_node(&node->children[i]);
		}

		//Set boundaries for the children
		node->children[0].min_x = node->min_x; node->children[0].max_x = (node->min_x + node->max_x)/2;
		node->children[0].min_y = node->min_y; node->children[0].max_y = (node->min_y + node->max_y)/2;
		
		node->children[1].min_x = (node->min_x + node->max_x)/2; node->children[1].max_x = node->max_x;
		node->children[1].min_y = node->min_y; node->children[1].max_y = (node->min_y + node->max_y)/2;

		node->children[2].min_x = node->min_x; node->children[2].max_x = (node->min_x + node->max_x)/2;
		node->children[2].min_y = (node->min_y + node->max_y)/2; node->children[2].max_y = node->max_y;

		node->children[3].min_x = (node->min_x + node->max_x)/2; node->children[3].max_x = node->max_x;
		node->children[3].min_y = (node->min_y + node->max_y)/2; node->children[3].max_y = node->max_y;
		
		//Put old particle into the appropriate child
		place_particle(node->particle, node);

		//Put new particle into the appropriate child
		place_particle(new_particle, node);

		//It now has children
		node->has_children = 1;
	} 
	//Add the new particle to the appropriate children
	else {
		//Put new particle into the appropriate child
		place_particle(new_particle, node);
	}
}


/*
 * Puts a particle in the right child of a node with children.
 */
void place_particle(int particle, struct node_t *node) {
	if(x[particle] <= (node->min_x + node->max_x)/2 && y[particle] <= (node->min_y + node->max_y)/2) {
		put_particle_in_tree(particle, &node->children[0]);
	} else if(x[particle] > (node->min_x + node->max_x)/2 && y[particle] < (node->min_y + node->max_y)/2) {
		put_particle_in_tree(particle, &node->children[1]);
	} else if(x[particle] < (node->min_x + node->max_x)/2 && y[particle] > (node->min_y + node->max_y)/2) {
		put_particle_in_tree(particle, &node->children[2]);
	} else {
		put_particle_in_tree(particle, &node->children[3]);
	}
}

/*
 * Sets initial values for a new node
 */
void set_node(struct node_t *node) {
	node->has_particle = 0;
	node->has_children = 0;
}

/*
 * Frees memory for a node and its children recursively.
 */
void free_node(struct node_t *node){
	if(node->has_children){
		free_node(&node->children[0]);
		free_node(&node->children[1]);
		free_node(&node->children[2]);
		free_node(&node->children[3]);
		free(node->children);
	}
}

/*
 * Displays the boundaries of the tree using OpenGL and GLUT. 
 * Used for debugging purposes.
 */
void display_tree(struct node_t *node) {
	drawRect(node->min_x, node->max_x, node->min_y, node->max_y);
	if(node->has_children == 1) {
		for(int i = 0; i < 4; i++) {
			display_tree(&node->children[i]);
		}
	}
}

/*
 * Calculates the total mass for the node. It recursively updates the mass
 * of itself and all of its children.
 */
double calculate_mass(struct node_t *node) {
	if(!node->has_particle) {
		node->total_mass = 0;
	} else if(!node->has_children) {
		node->total_mass = mass[node->particle];  
	} else {
		node->total_mass = 0;
		for(int i = 0; i < 4; i++) {
			node->total_mass += calculate_mass(&node->children[i]);
		}
	}
	return node->total_mass;
}

/*
 * Calculates the x-position of the centre of mass for the 
 * node. It recursively updates the position of itself and 
 * all of its children.
 */
double calculate_center_of_mass_x(struct node_t *node) {
	if(!node->has_children) {
		node->c_x = x[node->particle];
	} else {
		node->c_x = 0;
		double m_tot = 0; 
		for(int i = 0; i < 4; i++) {
			if(node->children[i].has_particle) {
				node->c_x += node->children[i].total_mass*calculate_center_of_mass_x(&node->children[i]);
				m_tot += node->children[i].total_mass;
			}
		}
		node->c_x /= m_tot;
	}
	return node->c_x;
}

/*
 * Calculates the y-position of the centre of mass for the 
 * node. It recursively updates the position of itself and 
 * all of its children.
 */
double calculate_center_of_mass_y(struct node_t *node) {
	if(!node->has_children) {
		node->c_y = y[node->particle];
	} else {
		node->c_y = 0;
		double m_tot = 0; 
		for(int i = 0; i < 4; i++) {
			if(node->children[i].has_particle) {
				node->c_y += node->children[i].total_mass*calculate_center_of_mass_y(&node->children[i]);
				m_tot += node->children[i].total_mass;
			}
		}
		node->c_y /= m_tot;
	}
	return node->c_y;
}

/*
 * Calculates the forces in a time step of all particles in
 * the simulation using the Barnes Hut quad tree.
 */
void update_forces(){
	for(int i = 0; i < N; i++) {
		force_x[i] = 0;
		force_y[i] = 0;
		update_forces_help(i, root);
	}
}

/*
 * Help function for calculating the forces recursively
 * using the Barnes Hut quad tree.
 */
void update_forces_help(int particle, struct node_t *node) {
    //The node is a leaf node with a particle and not the particle itself
    if(!node->has_children && node->has_particle && node->particle != particle) {
        double r = sqrt((x[particle] - node->c_x)*(x[particle] - node->c_x) + (y[particle] - node->c_y)*(y[particle] - node->c_y));
        calculate_force(particle, node, r);
    }
    //The node has children
    else if(node->has_children) {
        //Calculate r and theta
        double r = sqrt((x[particle] - node->c_x)*(x[particle] - node->c_x) + (y[particle] - node->c_y)*(y[particle] - node->c_y));
        double theta = (node->max_x - node->min_x)/r;
        
        /* If the distance to the node's centre of mass is far enough, calculate the force,
         * otherwise traverse further down the tree
         */
        if(theta < 0.5){
            calculate_force(particle, node, r);
        } else {
            update_forces_help(particle, &node->children[0]);
            update_forces_help(particle, &node->children[1]);
            update_forces_help(particle, &node->children[2]);
            update_forces_help(particle, &node->children[3]);
        }
    }
}

/*
 * Calculates and updates the force of a particle from a node.
 */
void calculate_force(int particle, struct node_t *node, double r){
    double temp = -grav*mass[particle]*node->total_mass/((r + epsilon)*(r + epsilon)*(r + epsilon));
	force_x[particle] += (x[particle] - node->c_x)*temp;
	force_y[particle] += (y[particle] - node->c_y)*temp;
}

/*
 * Main function.
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
    
	//Initiate memory for the vectors
	x = (double *)malloc(N*sizeof(double));
	y = (double *)malloc(N*sizeof(double));
	u = (double *)malloc(N*sizeof(double));
	v = (double *)malloc(N*sizeof(double));
	force_x = (double *)calloc(N, sizeof(double));
	force_y = (double *)calloc(N, sizeof(double));
    mass = (double *)malloc(N*sizeof(double));

	//Set the initial values
	for(int i = 0; i < N; i++) {
		mass[i] = 1;
		double R = frand(0, L/4);
		double theta = frand(0, 2*M_PI);
		x[i] = L/2 + R*cos(theta);
		y[i] = W/2 + alpha*R*sin(theta);
		double R_prim = sqrt(pow(x[i] - L/2, 2) + pow(y[i] - W/2, 2));
		u[i] = -V*R_prim*sin(theta);
		v[i] = V*R_prim*cos(theta);
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
		print_time(start, stop);
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
