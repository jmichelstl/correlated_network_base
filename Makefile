CPP=g++
FLAGS=-std=c++11 -Werror -O3
LINK=-lgsl -lgslcblas

corrnet:
    $(CPP) $(FLAGS) src/make_correlated_network.cpp -o corrnet $(LINK) -pthread

cnclean:
    rm corrnet
