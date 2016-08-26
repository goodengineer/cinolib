/****************************************************************************
* Italian National Research Council                                         *
* Institute for Applied Mathematics and Information Technologies, Genoa     *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Author: Marco Livesu (marco.livesu@gmail.com)                             *
*                                                                           *
* Copyright(C) 2016                                                         *
* All rights reserved.                                                      *
*                                                                           *
* This file is part of CinoLib                                              *
*                                                                           *
* CinoLib is free software; you can redistribute it and/or modify           *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 3 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
****************************************************************************/
#include <cinolib/geometry/ray.h>

namespace cinolib
{

CINO_INLINE
Ray::Ray(const vec3d & p, const vec3d & dir)
{
    start     = p;
    direction = dir;
    direction.normalize();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::vector<Plane> Ray::to_planes() const
{
    vec3d n0(-direction.y(),  direction.x(),             0);
    vec3d n1(-direction.z(),              0, direction.x());
    vec3d n2(             0, -direction.z(), direction.y());

    std::vector<Plane> planes;
    if (n0.length() > 0) planes.push_back(Plane(start, n0));
    if (n1.length() > 0) planes.push_back(Plane(start, n1));
    if (n2.length() > 0) if (planes.size() < 2) planes.push_back(Plane(start, n2));
    assert(planes.size() == 2);

    return planes;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CINO_INLINE
const vec3d & Ray::dir() const
{
    return direction;
}


CINO_INLINE
const vec3d & Ray::begin() const
{
    return start;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CINO_INLINE
bool Ray::on_positive_half_space(const vec3d & p) const
{
    if ((p - start).dot(direction) >= 0) return true;
    return false;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CINO_INLINE
double Ray::dist_to_point(const vec3d & p) const
{
    vec3d u = direction;
    vec3d w = p - start;

    float cos_wv = w.dot(u);
    float cos_uu = u.dot(u);

    if (cos_wv <= 0.0) return start.dist(p);

    float b  = cos_wv / cos_uu;
    vec3d Pb = start + u*b;
    return (p-Pb).length();
}

}
