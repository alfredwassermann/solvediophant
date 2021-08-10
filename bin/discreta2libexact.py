#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import string
#from Numeric import *
"""
 Read a KM-file from DISCRETA and generate a matrix
 which is readable by libexact
"""

def PrintSystem(mat, rhs):
    """
        Print the matrix in libexact format
    """
    for i, r in enumerate(rhs):
        print "r", i+1, r

    for i in xrange(len(mat[0])):
        print "c", i+1

    i = 0
    for i, row in enumerate(mat):
        for j, e in enumerate(row):
            if e == 1:
                print "e", i+1, j+1
            elif e != 0:
                print >>sys.stderr, "Entry different from 0/1"
                sys.exit(1)


if __name__ == '__main__':
    """
     Find the number of rows and cols. This is the first line
     not beginning with %
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
    withRHS = False
    if dims[2] == 0:
        rhs = [1] * dims[0]
    else:
        withRHS = True
        rhs = [0] * dims[0]
    
    r = 0
    while r < dims[0]:
        line = sys.stdin.readline()
        row = map(int, line.strip().split())
        if withRHS:
            matrix[r] = [x for x in row[:-1]]
            rhs[r] = row[-1]
        else:
            matrix[r] = [x for x in row]
        r += 1
    
    sys.stderr.write("read %d lines\n" % r)
    PrintSystem(matrix, rhs)
