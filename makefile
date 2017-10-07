percolate: main.o site.o bond.o clusters.o
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -o percolate \
	main.o site.o bond.o clusters.o -lm
	
main.o: main.c
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -c main.c

site.o: site.c
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -c site.c

bond.o: bond.c
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -c bond.c

clusters.o: clusters.c
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -c clusters.c
