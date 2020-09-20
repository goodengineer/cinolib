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
#ifndef CINO_PREDICATES
#define CINO_PREDICATES

#include <cinolib/geometry/vec2.h>
#include <cinolib/geometry/vec3.h>

namespace cinolib
{
/* This file provides orient, incircle and in sphere predicates,
 * as well as additional predicates that build on top of them
 * to test point in segments/triangles/tetrahedra, and also to
 * test intersections between these entities in 2D and 3D. In the
 * default configuration, all these predicates are INEXACT, and the
 * basic orient, incircle and in sphere are basically equivalent to
 * the "fast" version of the Shewchuk predicates.
 *
 * *********************************************************************
 * IMPORTANT: to switch to EXACT PREDICATES, you must define the symbol
 * CINOLIB_USES_EXACT_PREDICATES at compilation time, and also add the
 * file <CINOLIB_HOME>/external/predicates/shewchuk.c in your project.
 * *********************************************************************
 *
 * Return values for the point_in_{segment | triangle | tet} predicates:
 * an integer flag which indicates exactly where, in the input simplex, the
 * point is located is returned. Note that a point typically belongs to
 * multiple sub-simplices. For example, a point coincident to a triangle
 * vertex belongs to a 0-dimensional simplex (the vertex), at least two
 * 2-dimensional simplices (its incident edges), and one 2-dimensional
 * simplex (the triangle). The integer flag points to the lowest dimensional
 * (sub) simplex that fully contains the point.
 *
 * WARNING: for degenerated elements such as zero length segments, zero area
 * triangles and zero volume tets, the lowest dimensional simplex containing
 * a point may not be unique. In these cases, only one of them will be returned.
 *
 * Return values for intersection tests: an integer flag which indicates
 * whether the input simplices are fully disjoint or intersect is returned.
 * Intersection can be of three types: (i) the simplices intersect only
 * at a shared sub-simplex, defining a valid simplicial complex (e.g. two
 * edges sharing a commong endpoint); (ii) the simplices intersect in a way
 * that does not define a valid simplicial complex; (iii) the simplices
 * intersect in some pathological way (e.g. two colinear triangles that
 * partially overlap, or two coplanar triangles that overlap).
 *
 * WARNING: intersections tests all assume that simplices are non degenerated.
 * If the code is compiled without defining the NDEBUG symbol, dedicated assertions
 * will make the program stop in case of zero length edges, zero area triangles,
 * or zero volume tets.
 */

// location of intersection points for point_in_{segment | triangle | tet}
// predicates. Elements' orders are compliant with the tables in:
//
//   #include <cinolib/standard_elements_tables.h>
typedef enum
{
    STRICTLY_OUTSIDE = 0,  // strictly outside the input simplex
    STRICTLY_INSIDE  = 1,  // strictly inside  the input simplex
    ON_VERT0         = 2,  // used for segs, tris and tets
    ON_VERT1         = 3,  // used for segs, tris and tets
    ON_VERT2         = 4,  // used for tris and tets
    ON_VERT3         = 5,  // used for tets
    ON_EDGE0         = 6,  // used for tris and tets
    ON_EDGE1         = 7,  // used for tris and tets
    ON_EDGE2         = 8,  // used for tris and tets
    ON_EDGE3         = 9,  // used for tets
    ON_EDGE4         = 10, // used for tets
    ON_EDGE5         = 11, // used for tets
    ON_FACE0         = 12, // used for tets
    ON_FACE1         = 13, // used for tets
    ON_FACE2         = 14, // used for tets
    ON_FACE3         = 15, // used for tets
}
PointInSimplex;

// intersection types
typedef enum
{
    DO_NOT_INTERSECT   = 0, // simplices do not intersect
    SIMPLICIAL_COMPLEX = 1, // simplices form a valid simplicial complex (i.e. they are coincident or share a sub-simplex)
    INTERSECT          = 2, // simplices intersect in a non conforming way
    OVERLAP            = 3, // for corner cases: simplices intersect and partially overlap
}                           // (e.g. colinear segments or coplanar triangles)
SimplexIntersection;

#ifdef CINOLIB_USES_EXACT_PREDICATES

/* Wrap of the popular geometric predicates described by Shewchuk in:
 *
 * Routines for Arbitrary Precision Floating-point Arithmetic and
 * Fast Robust Geometric Predicates
 *
 * WARNING: if you use these predicates, you should include in your
 * project <CINOLIB_HOME>/external/predicates/shewchuk.c and compile it,
 * otherwise the linker will not find an implementation for the methods
 * below
 */
extern "C"
{

float orient2d(const float * pa,
                const float * pb,
                const float * pc);

float orient3d(const float * pa,
                const float * pb,
                const float * pc,
                const float * pd);

float incircle(const float * pa,
                const float * pb,
                const float * pc,
                const float * pd);

float insphere(const float * pa,
                const float * pb,
                const float * pc,
                const float * pd,
                const float * pe);
}

#else

// These are equivalent to the "fast" version of Shewchuk's predicates. Hence are INEXACT
// geometric predicates solely based on the accuracy of the floating point system

CINO_INLINE
float orient2d(const float * pa,
                const float * pb,
                const float * pc);

CINO_INLINE
float orient3d(const float * pa,
                const float * pb,
                const float * pc,
                const float * pd);

CINO_INLINE
float incircle(const float * pa,
                const float * pb,
                const float * pc,
                const float * pd);

CINO_INLINE
float insphere(const float * pa,
                const float * pb,
                const float * pc,
                const float * pd,
                const float * pe);
#endif

// wrap of orient2d for cinolib points. Either exact or not depending on CINOLIB_USES_EXACT_PREDICATES
CINO_INLINE
float orient2d(const vec2f & pa,
                const vec2f & pb,
                const vec2f & pc);

// wrap of orient3d for cinolib points. Either exact or not depending on CINOLIB_USES_EXACT_PREDICATES
CINO_INLINE
float orient3d(const vec3f & pa,
                const vec3f & pb,
                const vec3f & pc,
                const vec3f & pd);

// wrap of incircle for cinolib points. Either exact or not depending on CINOLIB_USES_EXACT_PREDICATES
CINO_INLINE
float incircle(const vec2f & pa,
                const vec2f & pb,
                const vec2f & pc,
                const vec2f & pd);

// wrap of insphere for cinolib points. Either exact or not depending on CINOLIB_USES_EXACT_PREDICATES
CINO_INLINE
float insphere(const vec3f & pa,
                const vec3f & pb,
                const vec3f & pc,
                const vec3f & pd,
                const vec3f & pe);

// true if the area of the triangle p0-p1-p2 is zero
CINO_INLINE
bool points_are_colinear_2d(const vec2f & p0,
                            const vec2f & p1,
                            const vec2f & p2);

// true if the area of the triangle p0-p1-p2 is zero
CINO_INLINE
bool points_are_colinear_2d(const float * p0,
                            const float * p1,
                            const float * p2);

// true if the area of all the orthogonal 2d projections of the triangle p0-p1-p2 is zero
CINO_INLINE
bool points_are_colinear_3d(const vec3f & p0,
                           const vec3f & p1,
                           const vec3f & p2);

// true if the area of all the orthogonal 2d projections of the triangle p0-p1-p2 is zero
CINO_INLINE
bool points_are_colinear_3d(const float * p0,
                            const float * p1,
                            const float * p2);

// true if the volume of the tetrahedron p0-p1-p2-p3 is zero
CINO_INLINE
bool points_are_coplanar_3d(const vec3f & p0,
                            const vec3f & p1,
                            const vec3f & p2,
                            const vec3f & p3);

// true if the volume of the tetrahedron p0-p1-p2-p3 is zero
CINO_INLINE
bool points_are_coplanar_3d(const float * p0,
                            const float * p1,
                            const float * p2,
                            const float * p3);

// returns:
// ON_VERTi         if p coincides with the i-th vertex of s
// STRICTLY_INSIDE  if p lies inside segment s (endpoints excluded)
// STRICTLY_OUTSIDE otherwise
CINO_INLINE
PointInSimplex point_in_segment_2d(const vec2f & p,
                                   const vec2f & s0,
                                   const vec2f & s1);

// returns:
// ON_VERTi         if p coincides with the i-th vertex of s
// STRICTLY_INSIDE  if p lies inside segment s (endpoints excluded)
// STRICTLY_OUTSIDE otherwise
CINO_INLINE
PointInSimplex point_in_segment_2d(const float * p,
                                   const float * s0,
                                   const float * s1);

// returns:
// ON_VERTi         if p coincides with the i-th vertex of s
// STRICTLY_INSIDE  if p lies inside segment s (endpoints excluded)
// STRICTLY_OUTSIDE otherwise
CINO_INLINE
PointInSimplex point_in_segment_3d(const vec3f & p,
                                   const vec3f & s0,
                                   const vec3f & s1);

// returns:
// ON_VERTi         if p coincides with the i-th vertex of s
// STRICTLY_INSIDE  if p lies inside segment s (endpoints excluded)
// STRICTLY_OUTSIDE otherwise
CINO_INLINE
PointInSimplex point_in_segment_3d(const float * p,
                                   const float * s0,
                                   const float * s1);

// returns:
// ON_VERTi         if p coincides with the i-th vertex of t
// ON_EDGEj         if p lies inside the j-th edge of t (endpoints excluded)
// STRICTLY_INSIDE  if p lies inside triangle t (borders excluded)
// STRICTLY_OUTSIDE otherwise
CINO_INLINE
PointInSimplex point_in_triangle_2d(const vec2f & p,
                                    const vec2f & t0,
                                    const vec2f & t1,
                                    const vec2f & t2);

// returns:
// ON_VERTi         if p coincides with the i-th vertex of t
// ON_EDGEj         if p lies inside the j-th edge of t (endpoints excluded)
// STRICTLY_INSIDE  if p lies inside triangle t (borders excluded)
// STRICTLY_OUTSIDE otherwise
CINO_INLINE
PointInSimplex point_in_triangle_2d(const float * p,
                                    const float * t0,
                                    const float * t1,
                                    const float * t2);

// returns:
// ON_VERTi         if p coincides with the i-th vertex of t
// ON_EDGEj         if p lies inside the j-th edge of t (endpoints excluded)
// STRICTLY_INSIDE  if p lies inside triangle t (borders excluded)
// STRICTLY_OUTSIDE otherwise
CINO_INLINE
PointInSimplex point_in_triangle_3d(const vec3f & p,
                                    const vec3f & t0,
                                    const vec3f & t1,
                                    const vec3f & t2);

// returns:
// ON_VERTi         if p coincides with the i-th vertex of t
// ON_EDGEj         if p lies inside the j-th edge of t (endpoints excluded)
// STRICTLY_INSIDE  if p lies inside triangle t (borders excluded)
// STRICTLY_OUTSIDE otherwise
CINO_INLINE
PointInSimplex point_in_triangle_3d(const float * p,
                                    const float * t0,
                                    const float * t1,
                                    const float * t2);

// returns:
// ON_VERTi         if p coincides with the i-th vertex of t
// ON_EDGEj         if p lies inside the j-th edge of t (endpoints excluded)
// ON_FACEk         if p lies inside the k-th face of t (borders excluded)
// STRICTLY_INSIDE  if p lies inside tetrahedron t (borders excluded)
// STRICTLY_OUTSIDE otherwise
CINO_INLINE
PointInSimplex point_in_tet(const vec3f & p,
                            const vec3f & t0,
                            const vec3f & t1,
                            const vec3f & t2,
                            const vec3f & t3);

// returns:
// ON_VERTi         if p coincides with the i-th vertex of t
// ON_EDGEj         if p lies inside the j-th edge of t (endpoints excluded)
// ON_FACEk         if p lies inside the k-th face of t (borders excluded)
// STRICTLY_INSIDE  if p lies inside tetrahedron t (borders excluded)
// STRICTLY_OUTSIDE otherwise
CINO_INLINE
PointInSimplex point_in_tet(const float * p,
                            const float * t0,
                            const float * t1,
                            const float * t2,
                            const float * t3);

// returns:
// DO_NOT_INTERSECT     if segments are fully disjoint
// SIMPLICIAL_COMPLEX   if segments coincide or intersect at a shared vertex
// INTERSECT            if segments intersect at an inner point (for s0, s1, or both)
// OVERLAP              if segments are colinear and partially overlapped
CINO_INLINE
SimplexIntersection segment_segment_intersect_2d(const vec2f & s00,
                                                 const vec2f & s01,
                                                 const vec2f & s10,
                                                 const vec2f & s11);

// returns:
// DO_NOT_INTERSECT     if segments are fully disjoint
// SIMPLICIAL_COMPLEX   if segments coincide or intersect at a shared vertex
// INTERSECT            if segments intersect at an inner point (for s0, s1, or both)
// OVERLAP              if segments are colinear and partially overlapped
CINO_INLINE
SimplexIntersection segment_segment_intersect_2d(const float * s00,
                                                 const float * s01,
                                                 const float * s10,
                                                 const float * s11);

// returns:
// DO_NOT_INTERSECT     if segments are fully disjoint
// SIMPLICIAL_COMPLEX   if segments coincide or intersect at a shared vertex
// INTERSECT            if segments intersect at an inner point (for s0, s1, or both)
// OVERLAP              if segments are colinear and partially overlapped
CINO_INLINE
SimplexIntersection segment_segment_intersect_3d(const vec3f & s00,
                                                 const vec3f & s01,
                                                 const vec3f & s10,
                                                 const vec3f & s11);

// returns:
// DO_NOT_INTERSECT     if segments are fully disjoint
// SIMPLICIAL_COMPLEX   if segments coincide or intersect at a shared vertex
// INTERSECT            if segments intersect at an inner point (for s0, s1, or both)
// OVERLAP              if segments are colinear and partially overlapped
CINO_INLINE
SimplexIntersection segment_segment_intersect_3d(const float * s00,
                                                 const float * s01,
                                                 const float * s10,
                                                 const float * s11);

// returns:
// DO_NOT_INTERSECT     if s and t are fully disjoint
// SIMPLICIAL_COMPLEX   if s is an edge of t, or s is degenerate and coincides with a vertex of t
// INTERSECT            if s and t intersect and do not form a valid simplex
CINO_INLINE
SimplexIntersection segment_triangle_intersect_2d(const vec2f & s0,
                                                  const vec2f & s1,
                                                  const vec2f & t0,
                                                  const vec2f & t1,
                                                  const vec2f & t2);

// returns:
// DO_NOT_INTERSECT     if s and t are fully disjoint
// SIMPLICIAL_COMPLEX   if s is an edge of t, or s is degenerate and coincides with a vertex of t
// INTERSECT            if s and t intersect and do not forma a valid simplex
CINO_INLINE
SimplexIntersection segment_triangle_intersect_2d(const float * s0,
                                                  const float * s1,
                                                  const float * t0,
                                                  const float * t1,
                                                  const float * t2);

// returns:
// DO_NOT_INTERSECT     if s and t are fully disjoint
// SIMPLICIAL_COMPLEX   if s is an edge of t, or s is degenerate and coincides with a vertex of t
// INTERSECT            if s and t intersect and do not form a valid simplex
CINO_INLINE
SimplexIntersection segment_triangle_intersect_3d(const vec3f & s0,
                                                  const vec3f & s1,
                                                  const vec3f & t0,
                                                  const vec3f & t1,
                                                  const vec3f & t2);

// returns:
// DO_NOT_INTERSECT     if s and t are fully disjoint
// SIMPLICIAL_COMPLEX   if s is an edge of t, or s is degenerate and coincides with a vertex of t
// INTERSECT            if s and t intersect and do not form a valid simplex
CINO_INLINE
SimplexIntersection segment_triangle_intersect_3d(const float * s0,
                                                  const float * s1,
                                                  const float * t0,
                                                  const float * t1,
                                                  const float * t2);

// returns:
// DO_NOT_INTERSECT     if s and t are fully disjoint
// SIMPLICIAL_COMPLEX   if s is an edge of t, or s is degenerate and coincides with a vertex of t
// INTERSECT            if s and t intersect and do not form a valid simplex
CINO_INLINE
SimplexIntersection segment_tet_intersect_3d(const vec3f & s0,
                                             const vec3f & s1,
                                             const vec3f & t0,
                                             const vec3f & t1,
                                             const vec3f & t2,
                                             const vec3f & t3);

// returns:
// DO_NOT_INTERSECT     if s and t are fully disjoint
// SIMPLICIAL_COMPLEX   if s is an edge of t, or s is degenerate and coincides with a vertex of t
// INTERSECT            if s and t intersect and do not form a valid simplex
CINO_INLINE
SimplexIntersection segment_tet_intersect_3d(const float * s0,
                                             const float * s1,
                                             const float * t0,
                                             const float * t1,
                                             const float * t2,
                                             const float * t3);

// returns:
// DO_NOT_INTERSECT     if triangles are fully disjoint
// SIMPLICIAL_COMPLEX   if triangles coincide or intersect at a shared sub-simplex
// INTERSECT            if triangles intersect without making a valid simplcial complex
CINO_INLINE
SimplexIntersection triangle_triangle_intersect_2d(const vec2f & t00,
                                                   const vec2f & t01,
                                                   const vec2f & t02,
                                                   const vec2f & t10,
                                                   const vec2f & t11,
                                                   const vec2f & t12);

// returns:
// DO_NOT_INTERSECT     if triangles are fully disjoint
// SIMPLICIAL_COMPLEX   if triangles coincide or intersect at a shared sub-simplex
// INTERSECT            if triangles intersect without making a valid simplcial complex
CINO_INLINE
SimplexIntersection triangle_triangle_intersect_2d(const float * t00,
                                                   const float * t01,
                                                   const float * t02,
                                                   const float * t10,
                                                   const float * t11,
                                                   const float * t12);

// returns:
// DO_NOT_INTERSECT     if triangles are fully disjoint
// SIMPLICIAL_COMPLEX   if triangles coincide or intersect at a shared sub-simplex
// INTERSECT            if triangles intersect without making a valid simplcial complex
CINO_INLINE
SimplexIntersection triangle_triangle_intersect_3d(const vec3f & t00,
                                                   const vec3f & t01,
                                                   const vec3f & t02,
                                                   const vec3f & t10,
                                                   const vec3f & t11,
                                                   const vec3f & t12);

// returns:
// DO_NOT_INTERSECT     if triangles are fully disjoint
// SIMPLICIAL_COMPLEX   if triangles coincide or intersect at a shared sub-simplex
// INTERSECT            if triangles intersect without making a valid simplcial complex
CINO_INLINE
SimplexIntersection triangle_triangle_intersect_3d(const float * t00,
                                                   const float * t01,
                                                   const float * t02,
                                                   const float * t10,
                                                   const float * t11,
                                                   const float * t12);

// returns true if s0==s1
CINO_INLINE
bool segment_is_degenerate_2d(const vec2f & s0,
                              const vec2f & s1);

// returns true if s0==s1
CINO_INLINE
bool segment_is_degenerate_2d(const float * s0,
                              const float * s1);

// returns true if s0==s1
CINO_INLINE
bool segment_is_degenerate_3d(const vec3f & s0,
                              const vec3f & s1);

// returns true if s0==s1
CINO_INLINE
bool segment_is_degenerate_3d(const float * s0,
                              const float * s1);

// returns true if t0, t1 and t2 are colinear
CINO_INLINE
bool triangle_is_degenerate_2d(const vec2f & t0,
                               const vec2f & t1,
                               const vec2f & t2);

// returns true if t0, t1 and t2 are colinear
CINO_INLINE
bool triangle_is_degenerate_2d(const float * t0,
                               const float * t1,
                               const float * t2);

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// returns true if t0, t1 and t2 are colinear
CINO_INLINE
bool triangle_is_degenerate_3d(const vec3f & t0,
                               const vec3f & t1,
                               const vec3f & t2);

// returns true if t0, t1 and t2 are colinear
CINO_INLINE
bool triangle_is_degenerate_3d(const float * t0,
                               const float * t1,
                               const float * t2);

// returns true if t0, t1, t2 and t3 are coplanar
CINO_INLINE
bool tet_is_degenerate(const vec3f & t0,
                       const vec3f & t1,
                       const vec3f & t2,
                       const vec3f & t3);

// returns true if t0, t1, t2 and t3 are coplanar
CINO_INLINE
bool tet_is_degenerate(const float * t0,
                       const float * t1,
                       const float * t2,
                       const float * t3);

// returns true if v0 and v1 are equal
CINO_INLINE
bool vec_equals_2d(const float * v0,
                   const float * v1);

// returns true if v0 and v1 are equal
CINO_INLINE
bool vec_equals_3d(const float * v0,
                   const float * v1);
}

#ifndef  CINO_STATIC_LIB
#include "predicates.cpp"
#endif

#endif // CINO_PREDICATES
