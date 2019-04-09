#ifdef _WIN32
#	include "GL/freeglut.h"
//#elif def __APPLE__
//#	include "GLUT/glut.h"
#else 
#	include "GLUT/glut.h"
#endif

/*
	author: Mu.H., based on works from K.Winsdor.
	date: Mu, 19 - 3 - 29 and earlier

	All right reserved.*/


#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <limits>

#include "data.h"
#include "matrix.h"
#include "object.h"

#define mat dmatrix_t 

using namespace std;

void OnDisplay();
void OnKeyboard(unsigned char key, int x, int y);
void Draw();

void initGLUT(int argc, char** argv, window_t& window) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(window.width, window.height);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Mu's World");

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glShadeModel(GL_FLAT);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	glutDisplayFunc(OnDisplay);
	glutKeyboardFunc(OnKeyboard);
}

void SetCurrentColorX(unsigned int r, unsigned int g, unsigned int b) {
	foregroundColor = Colour(r,g,b);
}

void SetPixelX(window_t& window, int i, int j) {
	if (i >= window.width || j >= window.height)
		return;

	unsigned int index = 3 * (j * window.width + i);
	frame[index] = (int)(255 * foregroundColor.r);
	frame[index + 1] = (int)(255 * foregroundColor.g);
	frame[index + 2] = (int)(255 * foregroundColor.b);
}

void OnDisplay() {
	memset(frame, 255, windowH*windowW* 3);
	Draw();

	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(windowW, windowH, GL_RGB, GL_UNSIGNED_BYTE, (GLubyte*)frame);
	glutSwapBuffers();
	glFlush();
}

void QuitX() {
	exit(0);
}

Camera C;
window_t Window = window_t(windowW, windowH);
light_t light;
Colour background;

void testMatrix ( );
void testColour( );
int main(int argc, char** argv) {
	/* set the background color */

	background = Colour(0.0, 0.0, 0.0);

	/* set up light position, intensity, and color */

	dmatrix_t m;
	vec light_position = vec (Lx,Ly,Lz,1);

	Colour light_intensity(1.0, 1.0, 1.0);
	Colour light_color(1.0, 1.0, 1.0);
	light = *build_light(&light, light_position, light_color, light_intensity);

	/* build display window and synthetic camera */

	C = Camera();

	/* build a sphere */

	m = mat();
	m.translated (0.0, 0.0, 0.0);

	Colour specular_color = Colour(1.0, 1.0, 1.0);
	Colour diffuse_color = Colour(0.0, 0.0, 1.0);
	Colour ambient_color = Colour(0.0, 0.0, 1.0);

	double specular_coeff = 0.4;
	double diffuse_coeff = 0.4;
	double ambient_coeff = 0.2;

	double f = 10.0;
	double reflectivity = 0.0;

	object[nobjects] = object_t(SPHERE, m, ambient_color, diffuse_color, specular_color, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity);

	initGLUT (argc, argv, Window);
	glutMainLoop ();
	return 0;
}

void testMatrix ( )
{
	mat A = mat (0,0,0);
	A += vec ( 1, 2, 3 );
	A += vec ( 0, 1, 4 );
	A = vec ( 2, 6, 0 ) + A;
	A.print ( );

	A.inverse( ).print ( );

	//system ( "pause" );
}

void testObject( )
{ }

void testColour( )
{
	Colour c = Colour(3,2,1);
	Colour d = c;
	c += d;
	c.print();
}

void Draw() {
	double aspect = ASPECT; /* Set near plane dimensions */
	double height = Near * tan(M_PI / 180.0 * THETA / 2.0);
	double width = height * aspect;

	dmatrix_t direction;
	int i, j;
	Colour pixel;

	for (i = 0; i < Window.width; i++) {
		for (j = 0; j < Window.height; j++) {
			direction = *ray_direction(C, &Window, height, width, (double)i, (double)j);
			pixel = shade(&light, object, C.E, &direction, pixel, background, 3);

			//pixel.r = 1;
			//pixel.g = 1;
			//pixel.b = 0;

			SetCurrentColorX((int)pixel.r, (int)pixel.g, (int)pixel.b);
			SetPixelX(Window, i, Window.height - (j + 1));
			delete_dmatrix(&direction);
		}
	}
}

void OnKeyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 'q':
		QuitX();
		break;
	}
}
