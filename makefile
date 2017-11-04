percolate: main.o outline.o percolate.o grid.o
	mpicc -std=gnu99 -Wall -pedantic -Werror -g -o percolate \
	main.o outline.o percolate.c grid.o -lm
	
main.o: main.c
	mpicc -std=gnu99 -Wall -pedantic -Werror -g -c main.c

percolate.o: percolate.c
	mpicc -std=gnu99 -Wall -pedantic -Werror -g -c percolate.c

grid.o: grid.c
	mpicc -std=gnu99 -Wall -pedantic -Werror -g -c grid.c

outline.o: outline.c
	mpicc -std=gnu99 -Wall -pedantic -Werror -g -c outline.c