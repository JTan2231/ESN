# ESN
Collection of linear algebra functions and matrix utilities for the purpose of implementing an [Echo State Network](https://en.wikipedia.org/wiki/Echo_state_network). This uses no external libraries aside from gnuplot for graphical display.

You can build and run `main.c` using gcc with the following command: `gcc main.c src/*.c -I./include/ -lm -o main`. Note that without gnuplot, graphical functions will not work. With gnuplot, a graph similar to the one below should show upon running `main`.

![alt text](https://github.com/JTan2231/ESN/blob/master/example_output.png?raw=true)
