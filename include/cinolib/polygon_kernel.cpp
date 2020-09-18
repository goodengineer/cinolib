/********************************************************************************
*  This file is part of CinoLib                                                 *
*  Copyright(C) 2016: Marco Livesu                                              *
*                                                                               *
*  The MIT License                                                              *
*                                                                               *
*  Permission is hereby granted, free of charge, to any person obtaining a      *
*  copy of this software and associated documentation files (the "Software"),   *
*  to deal in the Software without restriction, including without limitation    *
*  the rights to use, copy, modify, merge, publish, distribute, sublicense,     *
*  and/or sell copies of the Software, and to permit persons to whom the        *
*  Software is furnished to do so, subject to the following conditions:         *
*                                                                               *
*  The above copyright notice and this permission notice shall be included in   *
*  all copies or substantial portions of the Software.                          *
*                                                                               *
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR   *
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
*  FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE *
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      *
*  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS *
*  IN THE SOFTWARE.                                                             *
*                                                                               *
*  Author(s):                                                                   *
*                                                                               *
*     Marco Livesu (marco.livesu@gmail.com)                                     *
*     http://pers.ge.imati.cnr.it/livesu/                                       *
*                                                                               *
*     Italian National Research Council (CNR)                                   *
*     Institute for Applied Mathematics and Information Technologies (IMATI)    *
*     Via de Marini, 6                                                          *
*     16149 Genoa,                                                              *
*     Italy                                                                     *
*********************************************************************************/
#include <cinolib/polygon_kernel.h>
#include <cinolib/min_max_inf.h>

#ifdef CINOLIB_USES_BOOST
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry.hpp>
//
typedef boost::geometry::model::d2::point_xy<float>  BoostPoint;
typedef boost::geometry::model::polygon<BoostPoint>  BoostPolygon;
#endif

namespace cinolib
{

CINO_INLINE
float polygon_kernel(const std::vector<vec3d> & poly,   // will discard z component
                            std::vector<vec3d> & kernel) // z component will be zero
{
    std::vector<vec2d> poly_2d, kernel_2d;
    for(auto p : poly) poly_2d.push_back(vec2f(p));

    float area = polygon_kernel(poly_2d, kernel_2d);

    if (area > 0)
    {
        kernel.clear();
        for(auto p : kernel_2d) kernel.push_back(vec3f(p.x(), p.y(), 0));
    }
    return area;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef CINOLIB_USES_BOOST

CINO_INLINE
float polygon_kernel(const std::vector<vec2f> & poly,
                            std::vector<vec2f> & kernel)
{
    kernel.clear();
    if (poly.empty()) return 0;

    // define 2d axis aligned bbox
    vec2f min( inf_double,  inf_double);
    vec2f max(-inf_double, -inf_double);
    for(const vec2f & p : poly)
    {
        min = min.min(p);
        max = max.max(p);
    }
    float delta = min.dist(max);

    // define half spaces
    std::vector<BoostPolygon> half_spaces;
    for(uint i=0; i<poly.size(); ++i)
    {
        vec2f A   = poly.at(i),B   = poly.at((i+1)%poly.size());
        vec2f u = (B - A).normalize();  // edge direction
        vec2f v = vec2d(-u.y(), u.x()); // direction orthogonal to u (rotated CCW)
        A -= u * delta;
        B += u * delta;
        vec2f C = B + v * delta,D = A + v * delta;
        
        BoostPolygon p;
        boost::geometry::append(p, BoostPoint(A.x(), A.y()));
        boost::geometry::append(p, BoostPoint(B.x(), B.y()));
        boost::geometry::append(p, BoostPoint(C.x(), C.y()));
        boost::geometry::append(p, BoostPoint(D.x(), D.y()));
        boost::geometry::correct(p);
        half_spaces.push_back(p);
    }

    // define kernel as intersection of half-spaces
    BoostPolygon kernel_boost = half_spaces.front();
    for(uint i=1; i<half_spaces.size(); ++i)
    {
        std::vector<BoostPolygon> res;
        boost::geometry::intersection(kernel_boost, half_spaces.at(i), res);
        if (res.empty()) return 0; // no kernel
        assert(res.size()==1);     // this should always hold. I am intersecting convex polygons...
        kernel_boost = res.front();
    }

    // convert from BoostPolygon to vec2d array
    for(uint i=0; i<kernel_boost.outer().size()-1; ++i)
    {
        kernel.push_back(vec2d(boost::geometry::get<0>(kernel_boost.outer()[i]),
                               boost::geometry::get<1>(kernel_boost.outer()[i])));
    }

    // Boost uses CW ordered vertices, so here I make the kernel CCW
    std::reverse(kernel.begin(), kernel.end());

    return boost::geometry::area(kernel_boost);
}

#else

CINO_INLINE
float polygon_kernel(const std::vector<vec2f> &,
                            std::vector<vec2f> &)
{
    std::cerr << "ERROR : Boost Polygon disabled. Recompile defining symbol CINOLIB_USES_BOOST" << std::endl;
    exit(-1);
}

#endif

}
