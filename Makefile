CPP=g++
FLAGS=-std=c++11 -Werror -O3
LINK=-lgsl -lgslcblas

corrnet: src/triangle.o 
    $(CPP) $(FLAGS) src/make_correlated_network.cpp -o corrnet $(LINK) src/triangle.o -pthread 

src/triangle.o: src/triangle.c src/triangle.h
        gcc -O -DLINUX -I/usr/X11R6/include -L/usr/X11R6/lib -DTRILIBRARY -c -o src/triangle.o src/triangle.c  

cnclean:
    rm corrnet
