percolate: main.o site.o clusters.o
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -o percolate \
	main.o site.o clusters.o -lm -g
	
main.o: main.c
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -c main.c -g

site.o: site.c
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -c site.c -g

clusters.o: clusters.c
	gcc-7 -std=gnu99 -fopenmp -Wall -pedantic -Werror -c clusters.c -g