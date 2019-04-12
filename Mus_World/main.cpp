#include <memory>

//
// Created by Charles on 2019-04-11.
//

#include "main.h"
#include "matrix.h"
#include "object.cpp"
#include "windows.cpp"

#include <vector>
#include <iostream>
#include <cmath>

Windows options;
std::vector<char> frame;
std::vector<std::unique_ptr<Object>> objects;
std::vector<std::unique_ptr<Light>> lights;

void render()
{
	Vec3 *framebuffer = new Vec3[options.width * options.height];
	char *pix = &frame[nChars - 1];
	double scale = tan(deg2rad(options.fov * 0.5));
	double imageAspectRatio = options.width / (double)options.height;
	Vec3 orig(0);
	for (uint32_t i = 0; i < options.height; ++i) {
		for (uint32_t j = 0; j < options.width; ++j) {
			// generate primary ray direction
			double x = (2 * (j + 0.5) / (double)options.width - 1) * imageAspectRatio * scale;
			double y = (1 - 2 * (i + 0.5) / (double)options.height) * scale;
			Vec3 dir = normalize(Vec3(x, y, -1));
			Vec3 color = castRay(orig, dir, objects, lights, options, 0);
			*(pix) = color.z*255;
			*(pix - 1) = color.y*255;
			*(pix - 2) = color.x*255;
			pix -= 3;
		}
	}
}

void OnDisplay() {
	frame.resize(nChars, 0);
	render();

	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(W, H, GL_RGB, GL_UNSIGNED_BYTE, (GLubyte*)&frame[0]);
	glutSwapBuffers();
	glFlush();
}

void initGLUT(int argc, char** argv, Windows& window) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(W, H);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Mu's World");

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glShadeModel(GL_FLAT);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	glutDisplayFunc(OnDisplay);
	glutKeyboardFunc(OnKeyboard);
}

int main(int argc, char **argv)
{
	/* Set the background color and some other colors */
	Vec3 IvoryWrite = Vec3(250.0 / 255.0, 255.0 / 255.0, 240.0 / 255.0);
	Vec3 IvoryBlack = Vec3(41.0 / 255.0, 36.0 / 255.0, 33.0 / 255.0);
	Vec3 BlockRed = Vec3(178.0/255, 34.0/255, 34.0/255);
	Vec3 Pink = Vec3(255.0/255, 192.0/255, 203.0/255);

	/* Set up light position, intensity, and color */
	Light light0 = Light(Vec3(-10,10, 5), Vec3(0.8, 0.8, 0.8));
	Light light1 = Light(Vec3(10, -10, 5), Vec3(0.4, 0.4, 0.4));
	Light light2 = Light(Vec3(10,10,5), Vec3(0.8, 0.8, 0.8));
	Light light3 = Light(Vec3(-10,-10,5), Vec3(0.4, 0.4, 0.4));
	lights.push_back(make_unique<Light>(light0));
	lights.push_back(make_unique<Light>(light1));
	lights.push_back(make_unique<Light>(light2));
	lights.push_back(make_unique<Light>(light3));

	/* Build a plane */
	double depth_mesh = -5;
	double near_mesh = 5;
	double far_mesh = -50;
	double left_mesh = -50;
	double right_mesh = 50;
	Plane *mesh = new Plane(IvoryBlack, IvoryWrite, depth_mesh, near_mesh, far_mesh, left_mesh, right_mesh);
	mesh->materialType = LAMBERT;


	Sphere *sph1 = new Sphere(Vec3(-2, depth_mesh+2, -18), 2);
	sph1->diffuseColor = IvoryWrite;
	sph1->materialType = REFLECTION;

	Sphere *sph2 = new Sphere(Vec3(3, depth_mesh+3, -14), 3);
	sph2->diffuseColor = IvoryWrite;
	sph2->materialType = REFLECTION;

	Sphere *sph3 = new Sphere(Vec3(-6, depth_mesh+3, -12), 3);
	sph3->diffuseColor = IvoryWrite;
	sph3->materialType = REFLECTION;

	objects.push_back(unique_ptr<Sphere>(sph1));
	objects.push_back(unique_ptr<Sphere>(sph2));
	objects.push_back(unique_ptr<Sphere>(sph3));

	objects.push_back(unique_ptr<Plane>(mesh));

	// setting up options
	options = Windows(W, H, IvoryWrite, 4);

	// finally, render
	initGLUT(argc, argv, options);
	glutMainLoop();

	return 0;
}