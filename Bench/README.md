#  Benchmark problems for solvediophant
## Time to first solution
### `KM_PSL_2_23_7_8.in`:
- 0/1 problem, 2 solutions
- `./sd2 -bkz -beta30 -p18 Bench/KM_PSL_2_23_7_8.in`
- First solution after 4.8 sec, 30 Mio Loops
- `sd3` with ILDS early: First solution after 0.9 sec, 900k Loops, lds_k = 4
- `sd3` with ILDS late: First solution after 2.8 sec, 7 Mio Loops

### `n8_282.txt`:
- non 0/1 problem, 104 solutions
- `./sd2 -c10000 -bkz -beta30 -p18 Bench/n8_282.txt`
- First solution after 7 sec, 1.2 Mio loops
- `./sd3 -c10000 -bkz -beta30 -dfs Bench/n8_282.txt`
- First solutions after 2.6 Mio loops, 9 sec
- `./sd3 -c10000 -bkz -beta30 -lds Bench/n8_282.txt`
- First solutions after 1.0 Mio loops, 8.6 sec

### `KM_C11Id2X_t3_k4.txt`:
- `./sd2 -c10000 -bkz -beta30 -p18 Bench/KM_C11Id2X_t3_k4.txt`
- No solution found, but there must be many
- `./sd3 -c10000 -bkz -beta40 -p18 Bench/KM_C11Id2X_t3_k4.txt`
- Progressive BKZ, ILDS early: First solution after 1min5sec, 400k Loops, lds_k=2
- `./sd3 -c10000 -bkz -beta48 -lds Bench/KM_C11Id2X_t3_k4.txt`
- First solution after 25k loops, 1min3sec
- After removing conflicting columns from the LGS:
- `./sd3 -c10000 -bkz -beta48 -lds Bench/KM_C11Id2X_t3_k4.txt`
- First solution after 150 Mio loops, lds_k=4
- bkz beta60: First solution after 6 Mio loops, lds_k=3
- bkz beta80: First solution after 758k loops, lds_k=2

### `qdesign_2_8_4_217_q2.txt`:
- $2-(8,4,217)_2$ design, part of a large set for $N=3$, , group <singer^5, frob^2>
- `./sd2 -bkz -beta30 -p18 Bench/qdesign_2_8_4_217_q2.txt`
- No solution after several days
- `./sd3 -c1000 -bkz -beta60 -p18 Bench/qdesign_2_8_4_217_q2.txt`
- solved with progressive BKZ, ILDS early, threshold 1/6: lds_k = 0, 380 Mio loops, ca 30 min
- solved with progressive BKZ, ILDS early, threshold 0: lds_k = 2, ~~5.2 Mio loops~~, 812k loops
- solved with sd2 using dump from progressive BKZ in 380 Mio loops

### `qdesign_2_8_4_217_q2_LS_prob2.txt`:
- $2-(8,4,217)_2$ design, 2nd part of a large set for $N=3$, group <singer^5, frob>
- It may have no solution
- `./sd2 -bkz -beta30 -p18 Bench/qdesign_2_8_4_217_q2_LS_prob2.txt`
- No solution after several days
- `./sd3 -c1000 -bkz -beta60 -p18 Bench/qdesign_2_8_4_217_q2_LS_prob2.txt`

### !!! `arc_105_9_q13.txt`:
- non 0/1 problem
- at least one solution
- `./sd2 -bkz -beta30 -p18 Bench/arc_105_9_q13.txt`
- First solution after 609 Mio loops, 38 sec
- `./sd3 -bkz -beta30 -p18 Bench/arc_105_9_q13.txt`
- ILDS early: 85 Mio loops, 8 sec
- ILDS early, threshold 1/6: 54 Mio loops, 5 sec

### !!! `arc_28_3_q17.txt`:
- non 0/1 problem
- `./sd2 -bkz -beta60 -p18 Bench/arc_28_3_q17.txt`
-  No solution after 4 min 40sec
- `./sd3 -bkz -beta60 -p18 Bench/arc_28_3_q17.txt`
- ILDS early, threshold 1/6, 84 Mio loops, 14 sec, lds_k = 2

### `arc_204_12_q19.txt`:
- non 0/1 problem
- `./sd2 -bkz -beta60 -p18 Bench/arc_204_12_q19.txt`
- First solution after 2456 Mio loops, 2min 26 sec
- Small problem (39 cols), but first solution has 13 discrepancies
- `./sd3 -bkz -beta30 -p18 Bench/arc_204_12_q19.txt`
- ILDS early, threshold 1/6, 776 Mio loops, 1min 4 sec (11 discrepancies)

### `4_18_9_56.txt`:
- 0/1 problem
- `../sd3 -bkz -beta60 -dfs 4_18_9_56.txt`
- Solution after 353 min, 100354  Mio loops
- Old version: first solution after several weeks (according to Vedran and Mario)
- `../sd3 -bkz -beta60 -lds10 4_18_9_56.txt`
- Solution after 103 min, 12551 Mio loops, lds_k = 5

