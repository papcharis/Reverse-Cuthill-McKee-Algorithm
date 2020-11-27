# Reverse-Cuthill-McKee-Algorithm
to run sequential version:
	gcc -o seq sequential.c queue.c -O3 -g -std=c99 -lm
	./seq

to run parallel version:
	gcc -o par parallel.c queue.c -O3 -g -std=c99 -lm -fopenmp
	./par
