#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import string
#from Numeric import *
"""
---------------------------------------------------------
 Read a KM-file from DISCRETA and generate a matrix
 which is readable by Gurobi
---------------------------------------------------------
"""

def PrintMatrix(mat, rhs):
    """
        Print the matrix in Gurobi format
    """
    i = 0
    for row in mat:
        for k in xrange(len(row)):
            print '%d' % mat[i][k], 
        i += 1
        print ''
        
    

#if __name__ == '__main__':
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
    
size = dims[0]*dims[1]
matrix = [[]]*dims[0]
if dims[2]==0:
    rhs = [1]*dims[0]
else:
    rhs = [0]*dims[0]
    
r = 0
while r < dims[0]:
    line = sys.stdin.readline()
    line = line.strip()
    row = map(int,line.split())
    #for entry in row[:-1]:
    #    matrix[i] = entry
    matrix[r] = [x for x in row[:-1]]
    rhs[r] = row[-1]
    r += 1

sys.stderr.write("read %d lines\n" % r)
PrintMatrix(matrix, rhs)
