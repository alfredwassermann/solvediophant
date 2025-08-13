# Test QOBLIB market split instances

Apply solvediophant (`sd2`) to the market split instances of the
Quantum Optimization Benchmark Library (QOBLIB), see 
[Quantum Optimization Benchmark Library - The Intractable Decathlon](https://arxiv.org/abs/2504.03832).

Instances are in folder `01-marketsplit/instances`

Shell scripts (bash):

- `start_050.bash`: For all instances with $D=50$, find the first solution with LDS enumeration.
Output in folder `Results`: log file `050_lds.log` with all output, `*.sol` file for instance containing the solution.
- `start_all_sol.bash`: For *all* instances enumerate all solutions for each instances (with DFS enumeration). Output to folder `Results`: log file `all_dfs_full.log` with all output, `*.sol` file for instance containing the solutions.
- `start_parallel.bash`: For *all* instances search the first solution in parallel with LDS and DFS enumeration. As soon as one program has found a solution , the other program is terminated.  Output to folder `Results`: `parallel_run_dfs_full.log` and `parallel_run_lds_full.log` and `*.sol` file for instance containing the solution.

The bash scripts use the executable `../bin/sd2` which can be created
by typing `make` in the main folder.

