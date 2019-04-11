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
	glutCreateWindow("Assignment 3");

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glShadeModel(GL_FLAT);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	glutDisplayFunc(OnDisplay);
	glutKeyboardFunc(OnKeyboard);
}

int main(int argc, char **argv)
{
	// creating the scene (adding objects and lights)

	Sphere *sph1 = new Sphere(Vec3(-1, 0, -12), 2);
	sph1->materialType = DIFFUSE_AND_GLOSSY;
	sph1->diffuseColor = Vec3(0.6, 0, 0.8);
	Sphere *sph2 = new Sphere(Vec3(0.5, -1, -8), 1.5);
	sph2->ior = 1.5;
	sph2->materialType = REFLECTION_AND_REFRACTION;

	Sphere *sph3 = new Sphere(Vec3(1.5, 0, -10), 1);
	sph3->ior = 4;
	sph3->materialType = REFLECTION_AND_REFRACTION;

	objects.push_back(std::unique_ptr<Sphere>(sph3));
	objects.push_back(std::unique_ptr<Sphere>(sph2));
	objects.push_back(std::unique_ptr<Sphere>(sph1));

	Vec3 verts[4] = {{-5,-3,-6}, {5,-3,-6}, {5,-3,-16}, {-5,-3,-16}};
	uint32_t vertIndex[6] = {0, 1, 3, 1, 2, 3};
	Vec2 st[4] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
	MeshTriangle *mesh = new MeshTriangle(verts, vertIndex, 2, st);
	mesh->materialType = DIFFUSE_AND_GLOSSY;

	objects.push_back(std::unique_ptr<MeshTriangle>(mesh));

	//lights.push_back(std::unique_ptr<Light>(new Light(Vec3(-20, 70, 20), 0.5)));
	lights.push_back(std::make_unique<Light>(Vec3(30, 50, -12), 1));

	Vec3 IvoryWrite = Vec3(250.0 / 255.0, 255.0 / 255.0, 240.0 / 255.0);
	Vec3 Black = Vec3(0.0 / 255.0, 0.0 / 255.0, 0.0 / 255.0);
	// setting up options
	options = Windows(W, H, IvoryWrite, 10);

	// finally, render
	initGLUT(argc, argv, options);
	glutMainLoop();

	return 0;
}