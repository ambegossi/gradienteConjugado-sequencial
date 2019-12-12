IDIR =./iohb1.0

default: all

all: boeing
	gcc gradiente_conjugado.c iohb.o -Wno-pointer-to-int-cast -Wno-format-security -I$(IDIR) -O2 -o gradiente -fopenmp -pg
boeing:
	gcc -Wno-pointer-to-int-cast -Wno-format-security -I$(IDIR) -c $(IDIR)/iohb.c 
