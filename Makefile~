CFLAGS=-Wall -g 

sequential:main_sequential.c  physics.c
	gcc -std=c99 -g -o $@ $^ -lrt -lm 

LFLAGS=-lrt -lm
parallel:main_sequential.c physics.c
	mpicc $^ $(LFLAGS) $(CFLAGS) -o $@

allclean:
	rm sequential

cleanpar:
	rm parallel

