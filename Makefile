program1 = spec_conv
CC = gcc
CFLAGS = -Wall -lm -O2 -pedantic

.PHONY: default all clean
.DEFAULT_GOAL:=all

all: $(program1)

$(program1): %: %.c
	$(CC) $< $(CFLAGS) -o $@ 

clean:
	\rm -f $(program1)
