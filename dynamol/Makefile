ROOT=..
OBJS=stats.o sysinfo.o system.o sysmanip.o
INCLUDES=-I $(ROOT)/modules
DIRS=param struc coor

include $(ROOT)/Makefile.inc

sysinfo.o : system.o output
system.o: stats.o parser output vectors files


