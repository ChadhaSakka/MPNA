#ifndef MESH_H
#define MESH_H

// For a uniform mesh on [0,1], we might store 
// the number of points N, or some routines for indexing, etc.

static inline double computeDx(int N) {
    return 1.0 / (double)N;
}

#endif

