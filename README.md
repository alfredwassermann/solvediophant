# Solvediophant

This program searches for solutions of a system of linear equations of the form

```math
Ax = d, \quad A\in \mathbb{Z}^{m\times n}, d\in \mathbb{Z}^m,
```

where $x$ is a vector of non-negative integers. The unknowns $x_i$ can be bounded by $0\leq x_i\leq r_i\in\mathbb{Z}$, the default is $r_i=1$, i.e. binary variables.

The algorithm is based on lattice basis reduction followed by lattice enumeration (either depth first search or limited discrepancy search). The problem is known to be NP-complete, i.e. in general, solving instances might be difficult. The success rate of the algorithm is difficult to predict. In the search of combinatorial substructures, solvediophant successfully solves systems up to 1500 unknowns. The subfolder `QOBLIB` contains instances of the market split problem, see 
[Quantum Optimization Benchmark Library - The Intractable Decathlon](https://arxiv.org/abs/2504.03832).

Please, cite the algorithm with

```
@Article{	  wassermann98,
  author	= {Wassermann, Alfred},
  title		= {Finding simple $t$-designs with enumeration techniques},
  journal	= {Journal of Combinatorial Designs},
  volume	= {6},
  number	= {2},
  year		= {1998},
  pages		= {79--90},
  doi		= {10.1002/(SICI)1520-6610}
}

@InProceedings{10.1007/978-3-030-79987-8_2,
  author	= {Wassermann, Alfred},
  editor	= {Flocchini, Paola and Moura, Lucia},
  title		= {Search for Combinatorial Objects Using Lattice Algorithms -- Revisited},
  booktitle	= {Combinatorial Algorithms},
  year		= {2021},
  publisher	= {Springer International Publishing},
  address	= {Cham},
  pages		= {20--33},
  isbn		= {978-3-030-79987-8},
  doi		= {10.1007/978-3-030-79987-8_2}
}
```

## Compilation

The algorithm is implemented in C, depends on the multiprecision integer library GMP. For Intel/AMD 64 bit architecture supporting `AVX2` it is sufficient to start compilation by typing

```
make 
```

in the root folder. The executable `sd2` will be created in the subfolder `bin`. For other architectures, it is recommended to install `openBLAS` and adapt the `Makefile` accordingly. 

## Input format

Create a problem file, e.g. `Tests/in_sd.txt` 

```
% Correct is: 25 solution
7 10 1
1 0 1 0 1 0 2 0 0 0 6
1 0 0 1 0 1 0 2 0 0 7
0 1 0 1 1 0 0 0 2 0 7
0 1 1 0 0 1 0 0 0 2 8
1 1 0 0 0 0 0 0 0 0 2
0 0 1 1 0 0 0 0 0 0 3
0 0 0 0 1 1 0 0 0 0 4
BOUNDS 10
2 2 3 3 4 4 3 3 3 4 
```

- Initial lines starting with `%` are treated as comments
- The first non-comment line contains the size of the problem, here `7 10 1`, meaning
the matrix $A$ consists of $7$ rows and $10$ columns. The `1`  is optional, meaning that there is another column
containing the right hand side vector $d$.
- Now, $A$ and $d$ are listed row-wise, i.e.

```
    --                   --       -- --
    | 1 0 1 0 1 0 2 0 0 0 |       | 6 |
    | 1 0 0 1 0 1 0 2 0 0 |       | 7 |
    | 0 1 0 1 1 0 0 0 2 0 |       | 7 |
A = | 0 1 1 0 0 1 0 0 0 2 |,  d = | 8 |
    | 1 1 0 0 0 0 0 0 0 0 |       | 2 |
    | 0 0 1 1 0 0 0 0 0 0 |       | 3 |
    | 0 0 0 0 1 1 0 0 0 0 |       | 4 |
    --                   --       -- --
```

- Optionally, upper bounds for $x$ can be supplied by the two lines

```
BOUNDS 10
2 2 3 3 4 4 3 3 3 4 
```

If no upper bounds are given, the upper bound $1$ is assumed.

## Running solvediophant

Start the program (here from the root folder):

```
bin/sd2 Tests/in_sd.txt
```

The output to __stdout__ consists of all solutions

```
1 1 1 2 2 2 1 1 1 2 
1 1 2 1 1 3 1 1 2 1 
1 1 0 3 3 1 1 1 0 3 
1 1 3 0 0 4 1 1 3 0 
1 1 2 1 3 1 0 2 1 2 
1 1 3 0 2 2 0 2 2 1 
1 1 1 2 4 0 0 2 0 3 
1 1 1 2 0 4 2 0 2 1 
1 1 0 3 1 3 2 0 1 2 
0 2 2 1 2 2 1 2 1 1 
0 2 1 2 3 1 1 2 0 2 
0 2 3 0 1 3 1 2 2 0 
0 2 1 2 1 3 2 1 1 1 
0 2 0 3 2 2 2 1 0 2 
0 2 2 1 0 4 2 1 2 0 
0 2 3 0 3 1 0 3 1 1 
0 2 2 1 4 0 0 3 0 2 
0 2 0 3 0 4 3 0 1 1 
2 0 2 1 2 2 0 1 2 2 
2 0 1 2 3 1 0 1 1 3 
2 0 3 0 1 3 0 1 3 1 
2 0 0 3 4 0 0 1 0 4 
2 0 1 2 1 3 1 0 2 2 
2 0 0 3 2 2 1 0 1 3 
2 0 2 1 0 4 1 0 3 1 
```

Additionally. the solutions are written to the output file `solutions`. This can be changed with the parameter `-o`.

The output to __stderr__ gives additional information:

```
CPU supports AVX2
No reduction was chosen.
It is set to iterate=1.
You did not supply -c* or -C*. It is set to 1099511627776=2^40.
You did not supply -maxnorm*. It is set to 1.
You did not supply -scalelastline*. It is set to 1024.
Enumeration type is 'dfs'
LLL deltas:
	 delta_low   =0.550000
	 delta_med   =0.820000
	 delta_high  =0.990000
	 delta_higher=0.999000
Num. bounded variables = 10
No SELECTEDCOLUMNS found 
>> preprocess: remove zero-forced variables
>> preprocess: remove columns with entries that are too large
>> preprocess: find rows whose entries are too small
>> preprocess: find rows whose rhs can not be reached because of gcd
Rank <= 7
Num variables after preprocess: 10
The RHS is fixed !
upper bounds found. Max=12

0 2 2 1 0 4 2 1 2 0  !!
0 2 0 3 2 2 2 1 0 2  !!
1 1 0 3 1 3 2 0 1 2  !!
1 1 1 2 2 2 1 1 1 2  !!
First reduction successful
Second reduction successful
   log(D)= 17.653234
Third reduction successful

Dimension of solution space (k): 4 compared to (columns - rank): 4
Number of nonzero entries in the last row: 1
Max bit size: 4
Use floats in enumeration
Fq: 12.000000
Fd: 1584.158400
300.000 317.547 344.325 185.420 

Dual bounds:
2.525 2.364 1.442 1.781 

First non-zero entries:
6 2 2 1 : 11
Start enumeration at level=3
Prune_cs: 0
Prune_only_zeros: 1 of 14
Prune_hoelder: 0 of 0
Prune_hoelder interval: 5
Dual bounds: 1
Loops: 57
Total number of solutions: 25

total enumeration time: 0:00:00
```

Most notably, sometimes solutions are already found during lattice basis reduction. Output in stderr is

```
0 2 2 1 0 4 2 1 2 0  !!
0 2 0 3 2 2 2 1 0 2  !!
1 1 0 3 1 3 2 0 1 2  !!
1 1 1 2 2 2 1 1 1 2  !!
```

### Parameters

Parameters of `sd2` are listed by `sd2 -h`:

```
sd2 --- multiple precision version --- 
Usage:
	sd2 options inputfile
Options:
	 inputfile: file name or '-'  for stdin
	-iterate{num} do num LLL calls with delta=delta_high
	-bkz -beta{num} do BKZ with blocksize num
	-pbkz -beta{num} do progressive BKZ with max. blocksize num
	-tours{num} maximum number of tours in a bkz call
	-c{num} scale equations by num (default=1099511627776=2**40)
	-C{num} scale equations by 2^num (default=1099511627776=2**40)
	-scalelastline{num} scale last line by num (default=1024)
	-maxnorm* ???? default=1
	-delta_low{num} delta for first LLL reduction
	-delta_med{num} delta for second reduction
	-delta_high{num} delta for second reduction and for third reduction in case of -iterate
	-delta_higher{num} delta for bkz reduction
	-o{string} write solutions to file 'string' (default='solutions')
	-time{num} stop program after num seconds
	-t{num} stop program after num solutions. Overwritten by line '% stopafter num' in input file
	-d{num} Print progress report in enumeration after that many steps
	-silent do not write solutions to stdout and solution file
	-printntl write (shortened) lattice after third reduction in NTL format to stdout
	-dump If this flag is supplied, write lattice to file 'dump_lattice.b' after the reduction phase
	-restart{string} Read dumped lattice basis from file 'string' and jump to third reduction phase
	-double Do not use floats in enumeration
	-dfs Use depth first enumeration (default)
	-lds{num} Use LDS enumeration up to num discrepancies, otherwise use dfs (default=dfs)
Signals:
	 10: print lattice, e.g. kill -10 PID
	 12 Dump lattice to file 'dump_lattice.b', e.g. kill -12 PID
```

## License

The software is licensed under a variant of the 3-clause BSD license, see the file `LICENSE`.

