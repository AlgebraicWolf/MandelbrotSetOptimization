#ifndef MANDELBROT_CONFIG_H_
#define MANDELBROT_CONFIG_H_

// Window/Image size in px
const unsigned int WINDOW_WIDTH = 1200;
const unsigned int WINDOW_HEIGHT = 800;

// Position of bottom left corner of the initial image
const double STARTING_POSITION_X = -2.0;
const double STARTING_POSITION_Y = -1.0;

// Initial pixel step size
const double INITIAL_STEP_SIZE = 3.0 / WINDOW_WIDTH;

// Amount of steps per one move
const unsigned int MOVE = WINDOW_WIDTH / 50;

// Zoom coefficient (step *= ZOOM_COEFFICIENT)
const double ZOOM_COEFFICIENT = 0.9;

// Rendering settings
const double DIVERGENCE_RADIUS = 2.0;  // Radius of circle with centre in (0, 0) outside of which sequences of points are considered to diverge
const unsigned int MAX_ITER = 200;      // Number of iterations after which sequence of points is considered to converge

// Benchmark settings
const unsigned int BENCHMARK_ITERATIONS = 50;  // Number mandelbrot set rendering/display iterations
#define DISPLAY_BENCHMARK_FRAME                // Set this define in case benchmark frames should be displayed

#define SIMD_OPTIMIZATION  // Set this define in order to enable parallelized functions (Requires AVX2 and FMA3 on machine)

#define MULTITHREADING  // Set this define in order to enable multithreading

#ifdef MULTITHREADING
const unsigned int THREADS = 8;                                         // Number of threads for program to run
const unsigned int THREAD_IMAGE_PART_HEIGHT = WINDOW_HEIGHT / THREADS;  // Height of chunk that is calculated in one thread
#endif

#endif  // MANDELBROT_CONFIG_H_