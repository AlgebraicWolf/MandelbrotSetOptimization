#ifndef MANDELBROT_CONFIG_H_
#define MANDELBROT_CONFIG_H_

// Window/Image size in px
const unsigned int WINDOW_WIDTH = 900;
const unsigned int WINDOW_HEIGHT = 600;

// Position of bottom left corner of the initial image
const double STARTING_POSITION_X = -2.0;
const double STARTING_POSITION_Y = -1.0;

// Initial pixel step size
const double INITIAL_STEP_SIZE = 3.0 / WINDOW_WIDTH;

// Zoom coefficient (step *= ZOOM_COEFFICIENT)
const double ZOOM_COEFFICIENT = 0.9;

// Rendering settings
const double DIVERGENCE_RADIUS = 2.0; // Radius of circle with centre in (0, 0) outside of which sequences of points are considered to diverge
const unsigned int MAX_ITER = 50; // Number of iterations after which sequence of points is considered to converge

#endif // MANDELBROT_CONFIG_H_