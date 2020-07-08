build: main.o
	g++ bin/main.o -Wall -m64 -lsfml-graphics -lsfml-window -lsfml-system -mavx2 -Ofast -flto -funroll-loops -ftree-vectorize -pthread -o bin/mandelbrot

main.o: main.cpp
	g++ main.cpp -Wall -c -m64 -mavx2 -Ofast -flto -funroll-loops -pthread -o bin/main.o