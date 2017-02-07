rndmsiistat: rndmsiistat.c Makefile cmplx.c random.c linkedlist.c datastack.c Statanas0/statanas.c misc.c
	g++ rndmsiistat.c -o rndmsiistat -lm -g -L/usr/local/lib -lgsl -lgslcblas
