/*
 * File: graphics.h
 * ----------------
 * Contains definitions for the graphics functions.
 *
 * CME212 Assignment 5
 * Oliver Fringer
 * Stanford University
 *
 */
#ifndef _graphics_h
#define _graphics_h

void graphicsInit(int* argc, char** argv, void* display);

void reshape(int width, int height);

void idle(void);

void display(void); 

void drawPoints(double* x, double* y, int N);

void drawRect(double min_x, double max_x, double min_y, double max_y);

#endif
