/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/********************************************************************
 *                          IMPORTANT NOTE                          *
 *                                                                  *
 * The p4est_geometry interface will be removed shortly.            *
 * Please do NOT use this interface for newly written code.         *
 * It will be replaced with a generic transfinite blending scheme.  *
 ********************************************************************/

#ifndef P8EST_GEOMETRY_H
#define P8EST_GEOMETRY_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

/* Naming convention for different coordinate systems:
 * xyz              computational domain
 * abc              p4est vertex domain
 * rst              reference domain [-1,1]^3
 */

typedef struct p8est_geometry p8est_geometry_t;

/* Forward transformations return Jacobi determinant */
typedef void        (*p8est_geometry_X_t) (p8est_geometry_t * geom,
                                           p4est_topidx_t which_tree,
                                           const double abc[3],
                                           double xyz[3]);
typedef double      (*p8est_geometry_D_t) (p8est_geometry_t * geom,
                                           p4est_topidx_t which_tree,
                                           const double abc[3]);
typedef double      (*p8est_geometry_J_t) (p8est_geometry_t * geom,
                                           p4est_topidx_t which_tree,
                                           const double abc[3],
                                           double J[3][3]);

/* Inverse transformation to the reference element; returns -1 on error */
typedef int         (*p8est_geometry_R_t) (p8est_geometry_t * geom,
                                           p4est_topidx_t which_tree,
                                           const double txyz[3],
                                           double cabc[8][3],
                                           double abc[3], double rst[3]);

struct p8est_geometry
{
  const char         *name;     /* use prefixes to be namespace clean */
  p8est_geometry_X_t  X;
  p8est_geometry_D_t  D;
  p8est_geometry_J_t  J, Jit;   /* both return the determinant of J */
  p8est_geometry_R_t  R;        /* returns -1 on error */
};

/** Number of allowed Newton steps in p8est_geometry_R */
extern int          p8est_geometry_max_newton;

/** Compute the inverse transpose Jacobian by calling
 * the geom->J function and transpose inverting the result.
 * \return          The determinant of the Jacobian J (not of Jit).
 */
double              p8est_geometry_Jit (p8est_geometry_t * geom,
                                        p4est_topidx_t which_tree,
                                        const double abc[3],
                                        double Jit[3][3]);

/** Approximate the inverse transformation by Newton iterations.
 * The number of allowed Newton steps is p8est_geometry_max_newton.
 * \param [in] txyz Computational domain target coordinates.
 * \param [in] cabc Corners of a warped hexahedron not larger than the octree.
 * \param [out] abc Solution such that X(which_tree, abc) \approx txyz.
 * \param [out] rst Solution in reference coordinates wrt. the hexahedron cabc.
 * \return          The number of Newton iterations or -1 on failure.
 */
int                 p8est_geometry_R (p8est_geometry_t * geom,
                                      p4est_topidx_t which_tree,
                                      const double txyz[3],
                                      double cabc[8][3],
                                      double abc[3], double rst[3]);

/** The identity transformation.
 */
void                p8est_geometry_identity_X (p8est_geometry_t * geom,
                                               p4est_topidx_t which_tree,
                                               const double abc[3],
                                               double xyz[3]);

/** The Jacobi determinant of the identity transformation.
 * \return          The determinant of the Jacobian.
 */
double              p8est_geometry_identity_D (p8est_geometry_t * geom,
                                               p4est_topidx_t which_tree,
                                               const double abc[3]);

/** The Jacobian matrix of the identity transformation.
 * \return          The determinant of the Jacobian J.
 */
double              p8est_geometry_identity_J (p8est_geometry_t * geom,
                                               p4est_topidx_t which_tree,
                                               const double abc[3],
                                               double J[3][3]);

/** Create a geometry structure for the identity transformation.
 * \return          Geometry structure which must be freed with P4EST_FREE.
 */
p8est_geometry_t   *p8est_geometry_new_identity (void);

/** Create a geometry structure for the spherical shell of 24 trees.
 * This is suitable for forests obtained with p8est_connectivity_new_shell.
 * \param [in] R2   The outer radius of the shell.
 * \param [in] R1   The inner radius of the shell.
 * \return          Geometry structure which must be freed with P4EST_FREE.
 */
p8est_geometry_t   *p8est_geometry_new_shell (double R2, double R1);

/** Create a geometry structure for the solid sphere of 13 trees.
 * This is suitable for forests obtained with p8est_connectivity_new_sphere.
 * \param [in] R2   The outer radius of the sphere.
 * \param [in] R1   The outer radius of the inner shell.
 * \param [in] R0   The inner radius of the inner shell.
 * \return          Geometry structure which must be freed with P4EST_FREE.
 */
p8est_geometry_t   *p8est_geometry_new_sphere (double R2, double R1,
                                               double R0);

SC_EXTERN_C_END;

#endif /* !P8EST_GEOMETRY_H */
