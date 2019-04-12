//
// Created by Charles on 2019-04-11.
//

#include "windows.h"

Windows::Windows ( int w, int h, Vec3 backcolor, int maxRecur ) :width(w),height(h), backgroundColor(backcolor),
                                                                 maxDepth(maxRecur), bias(0.00001), fov(90)
{ }

Windows::Windows ()
{ }
