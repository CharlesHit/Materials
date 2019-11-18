# Materials

A customized ray tracer, provided a high-contrast, cool-tones Rembrandt style rendering. Stylistically different with traditional ray tracer!

![Alt text](sample_output/Screen\ Shot\ 2019-04-11 at 10.58.40 PM.png)

# How to Compile

Prerequisite: add freeGLUT or GLUT/OpenGL mannually.

# Materials I did in this program

Lambert Material: Pure diffuse reflection. I use a very very low Ks and high Kd to do it. It's almost the same as which used Gouraud shading techniques (As I did in my previous one). But of course, it's more natural.

Metal: Under Phong's model, it's hard to get the effect as Monte Carlo Algorithm can do. So I simply lower the Kd and Ks together, and use it approach matel's dark surface.

Steel: Same as Metal. Set the color to a light color and using high coeffient of Kd and low for Ks then you will get it.

Grass: The grass's color is decided by reflection and refraction. That is essentially different from Lambert material. In Grass there is no specular color, and the coeffients of them is 1.

# Reference

My most work followes _An improved illumination model for shaded display_, by T.Whitted.
