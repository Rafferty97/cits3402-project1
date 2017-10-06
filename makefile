percolate: main.o site.o bond.o clusters.o
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -o percolate \
	main.o site.o bond.o clusters.o -lm -g
	
main.o: main.c
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -c main.c -g

site.o: site.c
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -c site.c -g

bond.o: bond.c
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -c bond.c -g

clusters.o: clusters.c
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -c clusters.c -g