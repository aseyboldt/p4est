/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Copyright (C) 2012 Carsten Burstedde
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

#ifndef P4_TO_P8
#include <p4est_algorithms.h>
#include <p4est_io.h>
#else
#include <p8est_algorithms.h>
#include <p8est_io.h>
#endif

sc_array_t        *
p4est_deflate_quadrants (p4est_t * p4est, sc_array_t **data)
{
  const size_t       qsize = sizeof (p4est_qcoord_t);
  const size_t       dsize = p4est->data_size;
  size_t             qtreez, qz;
  sc_array_t        *qarr, *darr;
  p4est_topidx_t     tt;
  p4est_tree_t      *tree;
  p4est_quadrant_t  *q;
  p4est_qcoord_t    *qap;
  void              *dap;

  qarr = sc_array_new_size (qsize,
                            (P4EST_DIM + 1) * p4est->local_num_quadrants);
  qap = (p4est_qcoord_t *) qarr->array;
  darr = NULL;
  dap = NULL;
  if (data != NULL) {
    P4EST_ASSERT (dsize > 0);
    darr = sc_array_new_size (dsize, p4est->local_num_quadrants);
    dap = darr->array;
  }
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {    
    tree = p4est_tree_array_index (p4est->trees, tt);
    qtreez = tree->quadrants.elem_count;
    for (qz = 0; qz < qtreez; ++qz) {
      q = p4est_quadrant_array_index (&tree->quadrants, qz);
      *qap++ = q->x;
      *qap++ = q->y;
#ifdef P4_TO_P8
      *qap++ = q->z;
#endif
      *qap++ = (p4est_qcoord_t) q->level;
      if (data != NULL) {
        memcpy (dap, q->p.user_data, dsize);
        dap += dsize;
      }
    }
  }
  P4EST_ASSERT ((void *) qap ==
                qarr->array + qarr->elem_size * qarr->elem_count);
  if (data != NULL) {
    P4EST_ASSERT (dap == darr->array + darr->elem_size * darr->elem_count);
    *data = darr;
  }
  return qarr;
}

#ifndef P4_TO_P8

p4est_t           *
p4est_inflate (MPI_Comm mpicomm, p4est_connectivity_t *connectivity,
               const p4est_gloidx_t *global_first_quadrant,
               const p4est_gloidx_t *pertree,
               sc_array_t *quadrants, sc_array_t *data, void *user_pointer)
{
  const p4est_gloidx_t *gfq;
  int                 mpiret;
  int                 num_procs, rank;
  p4est_topidx_t      num_trees;
  p4est_t            *p4est;
#ifdef P4EST_DEBUG
  int                 p;
#endif
 
  P4EST_GLOBAL_PRODUCTION
    ("Into " P4EST_STRING "_inflate with min quadrants\n");

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));

  /* retrieve MPI information */
  mpiret = MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* assign some data members */
  p4est = P4EST_ALLOC_ZERO (p4est_t, 1);
  p4est->mpicomm = mpicomm;
  p4est->mpisize = num_procs;
  p4est->mpirank = rank;
  p4est->data_size = data == NULL ? 0 : data->elem_size;
  p4est->user_pointer = user_pointer;
  p4est->connectivity = connectivity;
  num_trees = connectivity->num_trees;

  /* create global first quadrant offsets */
  gfq = p4est->global_first_quadrant =
    P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  memcpy (p4est->global_first_quadrant, global_first_quadrant,
          (num_procs + 1) * sizeof (p4est_gloidx_t));
#ifdef P4EST_DEBUG
  for (p = 0; p < num_procs; ++p) {
    P4EST_ASSERT (gfq[p] <= gfq[p + 1]);
  }
#endif
  p4est->global_num_quadrants = gfq[num_procs];
  p4est->local_num_quadrants = (p4est_locidx_t) (gfq[rank + 1] - gfq[rank]);

  /* allocate memory pools */
  if (p4est->data_size > 0) {
    p4est->user_data_pool = sc_mempool_new (p4est->data_size);
  }
  else {
    p4est->user_data_pool = NULL;
  }
  p4est->quadrant_pool = sc_mempool_new (sizeof (p4est_quadrant_t));


  /* print more statistics */
  P4EST_VERBOSEF ("total local quadrants %lld\n",
                  (long long) p4est->local_num_quadrants);

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_inflate\n");

  return p4est;
}

#endif /* !P4_TO_P8 */
