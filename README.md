# TSP Parallelisation
Travelling Salesman Problem with Simulated Annealing, Parallelisation tests.

###Usage
call `bash batch.sh` to batch compile the .cpp files in `code/`, they will be placed into `execs/`.

Tests are pre-provided, but can be regenerated `bash ./tests/gen_tests.sh`. in the `tests/` directory calling ./gen_input n > input will create a file input with n cities in it.

To run the execs from the top directory call `./execs/<exec> < ./tests/<test>` to run the chosen exec with the chosen test. I'll automate testing once I've figured out how to handle MPI/OpenMP/Xeon Phi compilations.
