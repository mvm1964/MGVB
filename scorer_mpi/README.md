How to make shared library for portable scorer_mpi:

```
gcc -c -Wall -fpic -pthread scmpi.c sqlite3.c pepMass.c item_open.c bst2023.c compute_site_prob_open_optimized.c
gcc -shared -o libscmpi.so scmpi.o sqlite3.o pepMass.o item_open.o bst2023.o compute_site_prob_open_optimized.o
```

Then move libscmpi.so to desired location (/home/ubuntu on my multipass VM).
Then:

```
mpicc -L/home/ubuntu -g scorer_mpi_portable.c -o scorer_mpi_portable -lscmpi -lmpfr -lgmp -ldl -lm
```
