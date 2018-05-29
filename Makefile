
INST_DIR = $(PWD)/bin
SOURCES = heat2d.c

CC = gcc
CFLAGS = -Wall -c
LIBS = -lpmem -lomp

OBJS = $(SOURCES:%.c=%.o)

heat2d:heat2d.o
	$(CC) -o $@ heat2d.o $(LIBS)
	mkdir -p $(INST_DIR)
	mv heat2d $(INST_DIR)/heat2d

%.o:%.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -rf *~ \#_* *.o $(INST_DIR)
