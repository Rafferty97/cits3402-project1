percolate: main.o site2.o results.o
	gcc-7 -std=c99 -fopenmp -Wall -pedantic -Werror -o percolate \
	main.o site2.o results.o -lm -g
	
main.o: main.c
	gcc-7 -std=c99 -fopenmp -Wall -pedantic -Werror -c main.c -g
	
results.o: results.c
	gcc-7 -std=c99 -fopenmp -Wall -pedantic -Werror -c results.c -g

site2.o: site2.c
	gcc-7 -std=c99 -fopenmp -Wall -pedantic -Werror -c site2.c -g