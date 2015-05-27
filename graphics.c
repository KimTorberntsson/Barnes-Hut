/*
 * File: graphics.c
 * ----------------
 * Contains some simple graphics routines.
 *
 *
 */
#include "graphics.h"
#include <GLUT/glut.h>

void drawPoints(double* x, double* y, int N)
{
    int i;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Scale and translate the scene so coordinates match those in X.
    glPushMatrix();
    glScalef(2.0,2.0,1.0);    
    glTranslatef(-0.5,-0.5,0);
    
    // Begin drawing
    glBegin(GL_POINTS);
	glColor3f(1.0, 1.0, 1.0);	
	for(i = 0; i < N; i++)
	{
        glVertex2f(x[i],y[i]);
	}
	glEnd();

    glPopMatrix(); 
  	glutSwapBuffers();
}

void drawRect(double min_x, double max_x, double min_y, double max_y)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Scale and translate the scene so coordinates match those in X.
    glPushMatrix();
    glScalef(2.0,2.0,1.0);    
    glTranslatef(-0.5,-0.5,0);
    
    // Begin drawing
    glBegin(GL_LINE_STRIP);
    glColor3f(1.0, 0, 0);   
    glVertex2f(min_x, min_y);
    glVertex2f(max_x, min_y);
    glVertex2f(max_x, max_y);
    glVertex2f(min_x, max_y);
    glVertex2f(min_x, min_y);
    glEnd();

    glPopMatrix(); 
    glutSwapBuffers();
}


void graphicsInit(int* argc, char** argv, void* display)
{
    glutInit(argc, argv);
    
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(800, 800);
    
    glutCreateWindow("Barnes Hut");
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
}    

void reshape(int width, int height)
{
    glViewport(0, 0, width, height);
}

void idle(void)
{
    glutPostRedisplay();
}

