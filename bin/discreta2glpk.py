#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import string
#import glpk
import swiglpk as glpk

"""
---------------------------------------------------------
 Read a KM-file from DISCRETA and generate a matrix
 which is readable by GLPK
---------------------------------------------------------
"""

def PrintMatrix(mat, rhs, bounds):
    """
        Print the matrix in Gurobi format
    """
    i = 0
    for row in mat:
        for k in xrange(len(row)):
            print '%d' % mat[i][k], 
        i += 1
        print ''
        
    print
    print "RHS:", " ".join(map(str, rhs))
    print "BOUNDS:", " ".join(map(str, bounds))

def createGLPKProblem(matrix, rhs, bounds):
    num_rows = len(matrix)
    num_cols = len(matrix[0])
    
    lp = glpk.glp_create_prob()
    
    glpk.glp_set_prob_name(lp, "LP relax")
    glpk.glp_set_obj_dir(lp, glpk.GLP_MAX)
    glpk.glp_add_rows(lp, num_rows)
    
    for i in range(num_rows):
        glpk.glp_set_row_bnds(lp, i + 1, glpk.GLP_FX, rhs[i], rhs[i])

	glpk.glp_add_cols(lp, num_cols);
    for i in range(len(bounds)):
        if bounds[i] == 0:
            glpk.glp_set_col_bnds(lp, i + 1, glpk.GLP_FX, 0.0, bounds[i])
        else:
            glpk.glp_set_col_bnds(lp, i + 1, glpk.GLP_DB, 0.0, bounds[i])
        
    for i in range(num_cols):
        glpk.glp_set_obj_coef(lp, i + 1, 1.0);

    ia = glpk.intArray(num_rows * num_cols + 1)
    ja = glpk.intArray(num_rows * num_cols + 1)
    ar = glpk.doubleArray(num_rows * num_cols + 1)
    
    ne = 0
    for i, row in enumerate(matrix):
        for j, el in enumerate(row):
            if el != 0:
                ne += 1
                ia[ne] = i + 1
                ja[ne] = j + 1
                ar[ne] = el
    
    glpk.glp_load_matrix(lp, ne, ia, ja, ar)
  	#glpk.glp_term_out(glpk.GLP_OFF)

    params = glpk.glp_smcp()
    glpk.glp_init_smcp(params)
    params.presolve = glpk.GLP_ON
    #params.meth = glpk.GLP_DUALP
    print >>sys.stderr, "Start simplex solver"
    glpk.glp_simplex(lp, params)

    sol = []
    for i in range(num_cols):
        sol.append(glpk.glp_get_col_prim(lp, i + 1))
    print sol 

    Z = glpk.glp_get_obj_val(lp)
    print "optimal value", Z

    if True:
        for i in range(len(bounds)):
            glpk.glp_set_col_kind(lp, i + 1, glpk.GLP_BV);
        
        params = glpk.glp_iocp()
        glpk.glp_init_iocp(params)
        params.presolve = glpk.GLP_ON
        glpk.glp_intopt(lp, params)

        sol = []
        for i in range(num_cols):
            sol.append(glpk.glp_mip_col_val(lp, i + 1))
        print sol 

        Z = glpk.glp_mip_obj_val(lp)
        print "optimal value", Z
        #glpk.glp_interior(lp, None)
    
    print "DONE"
    
    
if __name__ == '__main__':
    """
    ---------------------------------------------------------
    Find the number of rows and cols. This is the first line
    not beginning with %
    ---------------------------------------------------------
    """
    foundDims = False
    while not foundDims:
        line = sys.stdin.readline()
        line = line.strip()
        if not line[0]==r'%':
            dims = map(int,line.split())
            foundDims = True
            if len(dims)==2:
                dims.append(0)
            sys.stderr.write("rows:%d, cols=%d, flag=%d\n" % (dims[0],dims[1],dims[2]))
        else:
            pass
            '''
                ToDo:
                number of solutions, number of iterations
            '''
    
    if not foundDims:
        sys.stderr.write("Could not find matrix size -> exit\n")
        sys.exit(1)
    
    size = dims[0] * dims[1]
    matrix = [[]] * dims[0]
    if dims[2] == 0:
        rhs = [1] * dims[0]
    else:
        rhs = [0] * dims[0]
    
    r = 0
    while r < dims[0] - 1:
        row = map(int, sys.stdin.readline().strip().split())
        matrix[r] = [x for x in row[:-1]]
        rhs[r] = row[-1]
        r += 1
    sys.stderr.write("read %d lines\n" % r)

    print >>sys.stderr, "Search BOUNDS"
    foundBounds = False
    while not foundBounds:
        words = sys.stdin.readline().strip().split()
        if len(words) > 0 and words[0] == 'BOUNDS':
            foundBounds = True
            num_bounds = int(words[1])
            break
    
    print >>sys.stderr, "Read BOUNDS"
    bounds = []
    while foundBounds and len(bounds) < num_bounds:
        line = sys.stdin.readline()
        line = line.strip()
        row = map(int, line.split())
        bounds = bounds + row
        
    print >>sys.stderr, "Reading done"
    #PrintMatrix(matrix, rhs, bounds)

    createGLPKProblem(matrix, rhs, bounds)
