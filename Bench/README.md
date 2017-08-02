#  Benchmark problems for solvediophant
## Time to first solution
* `KM_PSL_2_23_7_8.in`:
    - 0/1 problem, 2 solutions
    - `./sd2 -bkz -beta30 -p18 Bench/KM_PSL_2_23_7_8.in`
    - First solution after 4.8 sec, 30 Mio Loops
    - `sd3` with ILDS early: First solution after 0.9 sec, 2.6 Mio Loops
    - `sd3` with ILDS late: First solution after 2.8 sec, 7 Mio Loops
* `n8_282.txt`:
    - non 0/1 problem, 104 solutions
    - `./sd2 -c10000 -bkz -beta30 -p18 Bench/n8_282.txt`
    - First solution after 7 sec, 1.2 Mio loops
* `KM_C11Id2X_t3_k4.txt`:
    - `./sd2 -c10000 -bkz -beta30 -p18 Bench/KM_C11Id2X_t3_k4.txt`
    - No solution found, but there must be many
* `qdesign_2_8_4_217_q2.txt`:
    - $2-(8,4,217)_2$ design, part of a large set for $N=3$
* `arc_105_9_q13.txt`:
    - non 0/1 problem
    - at least one solution
    - ./sd2 -bkz -beta30 -p18 Bencharc_105_9_q13.txt`
    - First solution after 609 Mio loops, 38 sec
    - ILDS early: 93 Mio loops, 8 sec
* `arc_28_3_q17.txt`:
    - non 0/1 problem
