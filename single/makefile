CC=gcc
CFLAGS= -g
OBJS= convection.o

convection:  $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o convection

convection.o: convection.c
	$(CC) $(CFLAGS)  -c convection.c

clean :
	rm *.o
	rm convection
