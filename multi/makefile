CC=mpicc
CFLAGS= -g
OBJS= convection.o

convection:  $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -lm -o convection

convection.o: convection.c
	$(CC) $(CFLAGS) -c convection.c

clean :
	rm *.o
	rm convection
