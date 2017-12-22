#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import string

if __name__ == '__main__':
    mat = []
    rhs = []
    bounds = []
    num_rows = 0
    num_cols = 0

    is_first_entry = True
    for line in sys.stdin:
        c = line.strip().split(' ')
        if c[0] == 'r':
            rhs.append(int(c[2]))
            num_rows += 1
        elif c[0] == 'c':
            bounds.append(int(c[2]))
            num_cols += 1
        elif c[0] == 'e':
            if is_first_entry:
                for i in range(num_rows):
                    mat.append([0] * num_cols)
                is_first_entry = False

            mat[int(c[1]) - 1][int(c[2]) - 1] = 1
        elif c[0] == 'p':
            v = [0] * num_cols
            v[int(c[1]) - 1] = 1
            mat.append(v)
            rhs.append(1)
            num_rows += 1

    print '%\n%\n%'
    print num_rows, num_cols, 1
    for i in range(num_rows):
        print " ".join(map(str, mat[i])), ' ', rhs[i]

    print 'BOUNDS', num_cols
    print " ".join(map(str, bounds))

