CC=mpicc
CFLAGS= -g
OBJS= convection.o

convection:  $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -lm -o convection

convection.o: convection.c
	$(CC) $(CFLAGS) -c convection.c

playground: playground.o 
	gcc $(CFLAGS) playground.o -o playground

playground.o: playground.c
	gcc $(CFLAGS) -c playground.c
clean :
	rm *.o
	rm convection
