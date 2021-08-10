#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
    Implement the GLS of Theorem 4.1.(i) in 
    Direct product actions of groups on t-designs.
'''
import sys
import string

def binomial(n,k):
    if n<k or n<0 or k<0:
        print 'Input for binomial() not correct!'
        print 'n=%d k=%d' % (n,k)
        return 0L
    else:
        b = 1L
        for i in xrange(k):
            b *= (n-i)
            b /= (i+1)
        return b

if __name__ == '__main__':
    v = 9
    mat = [ [6, 1, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 5, 2, 0, 1, 1, 0, 0, 0, 0], 
            [0, 0, 4, 6, 0, 1, 2, 0, 0, 0], 
            [0, 0, 0, 0, 4, 1, 0, 2, 0, 0], 
            [0, 0, 0, 0, 0, 3, 4, 2, 4, 0], 
            [0, 0, 0, 0, 0, 0, 0, 2, 2, 6] 
          ]
    mat[0].append(2*binomial(v,5))
    mat[1].append(2*binomial(v,4)*(v-4))
    mat[2].append(2*binomial(v,3)*binomial(v-3,2))
    mat[3].append(2*binomial(v,4)*4)
    mat[4].append(2*binomial(v,3)*3*(v-3))
    mat[5].append(2*binomial(v,3)*3)
    
    print "%d %d 1" % (len(mat), len(mat[0])-1)
    for i in xrange(len(mat)):
        for j in xrange(len(mat[0])):
            print "%d " % mat[i][j],
        print ""

    print "BOUNDS %d" % (len(mat[0])-1)
    m = 1L
    for j in xrange(len(mat)):
        if mat[j][-1]>m:
            m = mat[j][-1]
    
    for j in xrange(len(mat[0])-1):
        print "%d " % m,
        
    print ""
    #print max(mat[:][0])
        
    
