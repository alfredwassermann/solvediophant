# Test problems for solvediophant

* KM_Id7_t2_k3.in:
0/1 problem, 30 solutions, easy

* KM_Id9_t2_k3.in:
0/1 problem, 840 solutions, easy

* KM_PGGL_2_32_t7_k8.in:
0/1 problem, 4996426 solutions, ca. 16 min
./solvediophant -bkz -beta100 -p40 -silent Tests/KM_PGGL_2_32_t7_k8.in (18 min, 13 min mit compute_w2, 10min with btmdx2 using compute_w2)
./sd3 -bkz -beta40 -silent Tests/KM_PGGL_2_32_t7_k8.in (18 min)
All solutions are found also with lds

* KM_PSL_2_23_7_8.in:
0/1 problem, 2 solutions, ca. 1 min
./sd3 -bkz -beta40 Tests/KM_PSL_2_23_7_8.in
All solutions are found also with lds

* unitals_v1.txt:
0/1 problems, 78 solutions
./sd3 -bkz -beta32 Tests/unitals_v1.txt
first solutions: 49 sec with dfs, 43 sec with lds

* unitals_v2.txt:
0/1 problems, >= 77 solutions
./sd3 -bkz -beta40 Tests/unitals_v2.txt

* n8_279.txt:
non 0/1 problem, 1848 solutions
Takes very long time
./sd3 -c10000 -bkz -beta40 Tests/n8_279.txt

* n8_282.txt:
non 0/1 problem, 104 solutions
./sd3 -c10000 -bkz -beta40 Tests/n8_282.txt

* KM_C11Id2X_t3_k4.txt:
solvediophant does not find solutions (at least in the first 5 minutes)
but there are many
With lds solutions are output quickly (lds_k=3)
./sd3 -lds -bkz -beta36 Tests/KM_C11Id2X_t3_k4.txt

* 4-arcs.txt
32705 solution. sd3: 10 sec on btmdx2
sd3:
    dfs: 32705 solutions
    lds: 32705 solutions
    17.2.2021: bug fix sd3, 32704 solutions
    18.2.2021: bug fixed. 32708 solutions

* in_sd.txt
non 0/1 problem 25 solutions.
A buggy version found 24 solutions.
With lds: 13 solutions
