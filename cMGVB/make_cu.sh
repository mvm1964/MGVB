#!/bin/bash
#gcc -g -Wall -c -fopenmp mgvb_2026_cuda.c
#nvcc -g -c -Xcompiler -fopenmp mgvb_cu.cu -o mgvb_cu.o -lm # -G option (debug device code) slows down 7x but enables cuda-gdb
#nvcc -dlink -Xcompiler -fopenmp mgvb_cu.o -o mgvb_dlink.o -lm
g++ -fopenmp -lgomp mgvb_dlink.o mgvb_2026_cuda.o scmpi.o sqlite3.o pepMass.o item_open.o bst2023.o compute_site_prob_open_optimized.o mgvb_cu.o -lmpfr -lgmp -ldl  -L/usr/local/cuda/lib64 -lcudart -lcuda -o mgvb_2026_cuda -lm


