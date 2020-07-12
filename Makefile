build: main.o
	g++ bin/main.o -Wall -m64 -lsfml-graphics -lsfml-window -lsfml-system -march=native -mssse3 -Ofast -flto -funroll-loops -ftree-vectorize -pthread -lsleef -o bin/mandelbrot

main.o: main.cpp
	g++ main.cpp -Wall -c -m64 -march=native -mssse3 -Ofast -flto -funroll-loops -pthread -o bin/main.o