OBJS = main.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

matchAdapters : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o matchAdapters

main.o : main.cpp
	$(CC) $(CFLAGS) main.cpp

clean:
	\rm *.o matchAdapters

