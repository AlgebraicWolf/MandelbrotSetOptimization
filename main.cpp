#include <SFML/Graphics.hpp>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <future>

#include "config.h"

#ifdef SIMD_OPTIMIZATION

#include <immintrin.h>
#include <sleef.h>

#endif

void renderMandelbrotSet(const double xstart, const double ystart, const double step, sf::Texture& renderTo);

void renderMandelbrotSetSmooth(const double xstart, const double ystart, const double step, sf::Texture& renderTo);

int main() {
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Mandelbrot set explorer", sf::Style::None);  // Window for the set to render in

    sf::Texture mandelbrotTexture;
    mandelbrotTexture.create(WINDOW_WIDTH, WINDOW_HEIGHT);  // Texture for the set rendering

    sf::Sprite mandelbrotSprite(mandelbrotTexture);  // Sprite with attached texture

    // Initializing coordinates and step with default values
    double x = STARTING_POSITION_X;
    double y = STARTING_POSITION_Y;
    double step = INITIAL_STEP_SIZE;

    bool update = true;  // Flag for set rerendering. Initial value is set to true
                         // In order to render image at the beginning

    bool smooth = true;  // Flag for smooth set rendering

    tm* current_date;                           // Pointer to structure for current time
    time_t current_time;                        // Variable for current time
    char filename[] = "00.00.0000:000000.png";  // String for current time

    sf::Clock clock;  // Clock for FPS measurement

    while (window.isOpen()) {  // Event loop
        sf::Event event;       // Variable for event parsing

        while (window.pollEvent(event)) {  // Loop through events
            switch (event.type) {          // Process current event
                case sf::Event::Closed:    // Upon closing the window...
                    window.close();        // Actually close it. The top level loop will not continue after this
                    break;

                case sf::Event::KeyPressed: {  // Upon button press...
                    switch (event.key.code) {  // Perform action linked to corresponding button
                        // Move up
                        case sf::Keyboard::W:
                        case sf::Keyboard::Up:
                            update = true;
                            y -= MOVE * step;
                            break;

                        // Move left
                        case sf::Keyboard::A:
                        case sf::Keyboard::Left:
                            update = true;
                            x -= MOVE * step;
                            break;

                        // Move down
                        case sf::Keyboard::S:
                        case sf::Keyboard::Down:
                            update = true;
                            y += MOVE * step;
                            break;

                        // Move right
                        case sf::Keyboard::D:
                        case sf::Keyboard::Right:
                            update = true;
                            x += MOVE * step;
                            break;

                        // Zoom in
                        case sf::Keyboard::I:
                            update = true;
                            // Perform zoom in a way that will preserve location of the current image centre
                            x += WINDOW_WIDTH * (1 - ZOOM_COEFFICIENT) * step / 2;
                            y += WINDOW_HEIGHT * (1 - ZOOM_COEFFICIENT) * step / 2;
                            step *= ZOOM_COEFFICIENT;
                            break;

                        // Zoom out
                        case sf::Keyboard::O:
                            update = true;
                            x -= WINDOW_WIDTH * (step / ZOOM_COEFFICIENT - step) / 2;
                            y -= WINDOW_HEIGHT * (step / ZOOM_COEFFICIENT - step) / 2;
                            step /= ZOOM_COEFFICIENT;
                            break;

                        // Toggle smoothing:
                        case sf::Keyboard::Space:
                            update = true;
                            smooth = !smooth;
                            break;

                        // Perform benchmark
                        case sf::Keyboard::B:
                            clock.restart();
                            for (unsigned int i = 0; i < BENCHMARK_ITERATIONS; i++) {
                                if (smooth) {  // Perform rendering with smoothing
                                    renderMandelbrotSetSmooth(x, y, step, mandelbrotTexture);
                                } else {  // Perform rendering without smoothing
                                    renderMandelbrotSet(x, y, step, mandelbrotTexture);
                                }

#ifdef DISPLAY_BENCHMARK_FRAME
                                window.clear();
                                window.draw(mandelbrotSprite);
                                window.display();
#endif
                            }
                            printf("FPS: %lf\n", static_cast<double>(BENCHMARK_ITERATIONS) / clock.getElapsedTime().asSeconds());
                            break;

                        // Save to file
                        case sf::Keyboard::Enter:
                            time(&current_time);
                            current_date = localtime(&current_time);
                            strftime(filename, sizeof(filename), "%d.%m.%Y:%H%M%S.png", current_date);
                            mandelbrotTexture.copyToImage().saveToFile(filename);
                            break;

                        // Stub for any other key
                        default:
                            break;
                    }
                } break;

                // Stub for the rest of the events
                default:
                    break;
            }
        }

        if (update) {  // In case image requires re-rendering
            update = false;
            if (smooth) {  // Perform rendering with smoothing
                renderMandelbrotSetSmooth(x, y, step, mandelbrotTexture);
            } else {  // Perform rendering without smoothing
                renderMandelbrotSet(x, y, step, mandelbrotTexture);
            }

            window.clear();
            window.draw(mandelbrotSprite);
            window.display();
        }
    }

    return 0;
}

#ifdef SIMD_OPTIMIZATION

inline __m128i colorizeVector(uint64_t iters_64, __m256d Res, __m256d Ims) {
    __m256d squaredRes = _mm256_mul_pd(Res, Res);
    __m256d Abs = _mm256_fmadd_pd(Ims, Ims, squaredRes);  // Calculate squared absolute value of pixel position

    __m128i integer_iters = _mm_setr_epi32(iters_64 & 0xFFFF, (iters_64 >> 16) & 0xFFFF, (iters_64 >> 32) & 0xFFFF, (iters_64 >> 48) & 0xFFFF);  // Load iteration counts into vector register
    __m256d iters = _mm256_cvtepi32_pd(integer_iters);                                                                                           // Convert iteration count to double representation
    __m256d two = _mm256_set1_pd(2.0);                                                                                                           // Mask for division by two
    __m256d one = _mm256_set1_pd(1.0);                                                                                                           // Mask for subtraction of one
    __m256d max = _mm256_set1_pd(255.0);                                                                                                         // Maximal value for pixel

    __m128i resultingColor = _mm_set1_epi32(0xFF000000);                                                 // Register for resulting color
    __m128i shuffleMask = _mm_setr_epi8(0, 0, 0, 0xFF, 4, 4, 4, 0xFF, 8, 8, 8, 0xFF, 12, 12, 12, 0xFF);  // Mask for building pixel colors out of values

    iters = _mm256_add_pd(iters, one);  // n + 1

    // Calculating color channel info
    __m256d argument = _mm256_div_pd(Sleef_log2d4_u35avx2(Sleef_logd4_u35avx2(Abs)), two);  // log2(log(Re * Re + Im * Im)) / 2
    argument = _mm256_add_pd(iters, argument);                                              // n + 1 + log2(log(Re * Re + Im * Im)) / 2

    __m256d cos_value = Sleef_fabsd4_avx2(Sleef_cosd4_u35avx2(argument));  // Absolute cosine value

    __m256d color_value = _mm256_mul_pd(cos_value, max);  // Calculate values for channel

    __m128i integerColorValue = _mm256_cvtpd_epi32(color_value);                                         // Cast values to uint32_t
    integerColorValue = _mm_or_si128(_mm_shuffle_epi8(integerColorValue, shuffleMask), resultingColor);  // Build pixel colors

    return integerColorValue;  // Return resulting color
}


uint32_t avgColorSSE(__m128i colors) {
    // Masks for extraction of all colors and conversion of uint8_t into uint32_t
    __m128i extractP1 = _mm_set_epi8(0x80, 0x80, 0x80, 3,
                                     0x80, 0x80, 0x80, 2,
                                     0x80, 0x80, 0x80, 1,
                                     0x80, 0x80, 0x80, 0);

    __m128i extractP2 = _mm_set_epi8(0x80, 0x80, 0x80, 7,
                                     0x80, 0x80, 0x80, 6,
                                     0x80, 0x80, 0x80, 5,
                                     0x80, 0x80, 0x80, 4);

    __m128i extractP3 = _mm_set_epi8(0x80, 0x80, 0x80, 11,
                                     0x80, 0x80, 0x80, 10,
                                     0x80, 0x80, 0x80, 9,
                                     0x80, 0x80, 0x80, 8);

    __m128i extractP4 = _mm_set_epi8(0x80, 0x80, 0x80, 15,
                                     0x80, 0x80, 0x80, 14,
                                     0x80, 0x80, 0x80, 13,
                                     0x80, 0x80, 0x80, 12);

    // Mask that merges all the uint32_t channel colors into the lowest uint32_t of SSE register; Each color is converted into uint8_t
    __m128i mergeTogether = _mm_set_epi8(0x80, 0x80, 0x80, 0x80,
                                         0x80, 0x80, 0x80, 0x80,
                                         0x80, 0x80, 0x80, 0x80,
                                         12, 8, 4, 0);

    // Extract all colors
    __m128i pixel1 = _mm_shuffle_epi8(colors, extractP1);
    __m128i pixel2 = _mm_shuffle_epi8(colors, extractP2);
    __m128i pixel3 = _mm_shuffle_epi8(colors, extractP3);
    __m128i pixel4 = _mm_shuffle_epi8(colors, extractP4);

    // Calculate sum of channels
    __m128i subsum1 = _mm_add_epi32(pixel1, pixel2);
    __m128i subsum2 = _mm_add_epi32(pixel3, pixel4);
    __m128i sum = _mm_add_epi32(subsum1, subsum2);

    // Integer division by four using binary shift
    __m128i expandedColor = _mm_srli_epi32(sum, 2);  // Divide by four
    __m128i finalColor = _mm_shuffle_epi8(expandedColor, mergeTogether);

    // Return resulting uint8_t 
    return _mm_extract_epi32(finalColor, 0);
}


__m128i getPointColorsAVX2(const double Re3, const double Re2, const double Re1, const double Re0, const double Im3, const double Im2, const double Im1, const double Im0) {
    __m256d Res = _mm256_set_pd(Re3, Re2, Re1, Re0);                                         // Initial real values of points
    __m256d Ims = _mm256_set_pd(Im3, Im2, Im1, Im0);                                         // Initial imaginary values of points
    uint64_t iters = 0;                                                                      // Initial values for iteration counter. Each 64-bit
    __m256d squaredDivergeceRadius = _mm256_set1_pd(DIVERGENCE_RADIUS * DIVERGENCE_RADIUS);  // Register used for divergence check

    __m256d curRes = _mm256_set1_pd(0);
    __m256d curIms = _mm256_set1_pd(0);  // Initial arrays of point positions

    __m256d finalRes = _mm256_set1_pd(0);
    __m256d finalIms = _mm256_set1_pd(0);  // Final positions after convergence have been checked

    __m256d two = _mm256_set1_pd(2);  // Mask for multiplication by two

    __m256d finalMask = _mm256_castsi256_pd(_mm256_set1_epi8(0xFF));  // Mask for keeping track of points for which divergence have already been proved
    __m256d copyMask;                                                 // Mask for element copying
    __m256d extractedRes;                                             // Re(z) for points diverging on current iteration
    __m256d extractedIms;                                             // Im(z) for points diverging on current iteration

    __m256d absolute;  // |z|^2
    __m256d mask;      // Comparison result

    __m256d squaredRes = _mm256_mul_pd(curRes, curRes);  // Squared Re(z) value. L: 4, TP: 0.5

    for (unsigned int i = 0; i < MAX_ITER; i++) {
        __m256d ResIms = _mm256_mul_pd(curRes, curIms);  // Re(z) * Im(z). L: 4, TP: 0.5

        curRes = _mm256_fmadd_pd(curRes, curRes, Res);  // Intermediate result for new Re(z). L: 4, TP: 0.5

        curRes = _mm256_fnmadd_pd(curIms, curIms, curRes);  // New Re(z). L: 4, TP: 0.5
        curIms = _mm256_fmadd_pd(two, ResIms, Ims);         // New Im(z). L: 4, TP: 0.5

        squaredRes = _mm256_mul_pd(curRes, curRes);  // Squared Re(z) value. L: 4, TP: 0.5

        absolute = _mm256_fmadd_pd(curIms, curIms, squaredRes);  // |z|^2. L: 4, TP: 0.5

        mask = _mm256_cmp_pd(absolute, squaredDivergeceRadius, _CMP_LE_OS);  // Perform comparison. L: 4, TP: 0.5

        uint64_t smallmask = _mm256_movemask_pd(mask);  // Translate mask into uint64_t. Usefult information is located in four lowest bytes. L: 2, TP: 1

        copyMask = _mm256_andnot_pd(mask, finalMask);     // Mask for element copying
        extractedRes = _mm256_and_pd(curRes, copyMask);   // Re(z) of elements to copy
        extractedIms = _mm256_and_pd(curIms, copyMask);   // Im(z) of elements to copy
        finalRes = _mm256_or_pd(extractedRes, finalRes);  // Copy Re(z)
        finalIms = _mm256_or_pd(extractedIms, finalIms);  // Copy Im(z)

        finalMask = _mm256_and_pd(finalMask, mask);  // Update mask of elements that have not been copied

        if (!smallmask)  // Check whether all points are diverging
            break;       // If so, exit loop

        smallmask = (smallmask & 0b0001) | ((smallmask & 0b0010) << 15) | ((smallmask & 0b0100) << 30) | ((smallmask & 0b1000) << 45);  // Restructure mask
        iters += smallmask;                                                                                                             // Perform increments where required
    }

    return colorizeVector(iters, finalRes, finalIms);  // Calculate colors of points and return
}


void renderMandelbrotSetToArrayAVX2(const double xstart, const double ystart, const double step, uint32_t* img, const unsigned int SEGMENT_HEIGHT) {
    for (unsigned int y = 0; y < SEGMENT_HEIGHT; y++) {
        for (unsigned int x = 0; x < WINDOW_WIDTH; x += 4) {  // Since one 256-bit AVX2 array can hold 4 pixels, we will iterate by four
            _mm_store_si128(reinterpret_cast<__m128i*>(img + y * WINDOW_WIDTH + x), getPointColorsAVX2(xstart + (x + 3) * step, xstart + (x + 2) * step, xstart + (x + 1) * step, xstart + x * step,
                                                                                                       ystart + y * step, ystart + y * step, ystart + y * step, ystart + y * step));
        }
    }
}


void renderMandelbrotSetToArraySmoothAVX2(const double xstart, const double ystart, const double step, uint32_t* img, const unsigned int SEGMENT_HEIGHT) {
    for (unsigned int y = 0; y < SEGMENT_HEIGHT; y++) {
        for (unsigned int x = 0; x < WINDOW_WIDTH; x++) {
            // Calculate colors of all points for supersampling using AVX2 vector extention
            __m128i colors = getPointColorsAVX2(xstart + x * step, xstart + x * step + step / 2, xstart + x * step, xstart + x * step + step / 2,
                                                ystart + y * step, ystart + y * step + step, ystart + y * step + step / 2, ystart + y * step + step / 2);

            img[y * WINDOW_WIDTH + x] = avgColorSSE(colors); // Use SSE and calculate average color. Store it.
        }
    }
}

#else

inline uint32_t colorize(unsigned int iter, double Re, double Im) {
    if (iter == MAX_ITER)
        return 0xFF000000;
    else {
        double argument = iter + 1 + log2(log(Re * Re + Im * Im)) / 2;                // Compute argument for cosine
        unsigned int color1 = static_cast<unsigned int>(fabs(cos(argument)) * 0xFF);  // Compute channel value and convert it to unsigned integer

        return 0xFF000000 | (color1 << 16) | (color1 << 8) | color1;  // Build pixel color
    }
}


inline uint32_t avgChanel(uint8_t c1, uint8_t c2, uint8_t c3, uint8_t c4) {
    return (static_cast<uint32_t>(c1) + c2 + c3 + c4) / 4; // Compute channel average value
} are using SIMD-optimized version

inline uint32_t avgColor(uint32_t RGBA1, uint32_t RGBA2, uint32_t RGBA3, uint32_t RGBA4) {
    return avgChanel(RGBA1 & 0xFF, RGBA2 & 0xFF, RGBA3 & 0xFF, RGBA4 & 0xFF) |
           (avgChanel((RGBA1 >> 8) & 0xFF, (RGBA2 >> 8) & 0xFF, (RGBA3 >> 8) & 0xFF, (RGBA4 >> 8) & 0xFF) << 8) |
           (avgChanel((RGBA1 >> 16) & 0xFF, (RGBA2 >> 16) & 0xFF, (RGBA3 >> 16) & 0xFF, (RGBA4 >> 16) & 0xFF) << 16) |
           (avgChanel((RGBA1 >> 24) & 0xFF, (RGBA2 >> 24) & 0xFF, (RGBA3 >> 24) & 0xFF, (RGBA4 >> 24) & 0xFF) << 24); // Using per-channel averaging, compute average color
}


uint32_t getPointColor(const double Re, const double Im) {
    double curRe = 0,
           curIm = 0,
           newRe = 0,
           newIm = 0;  // Initialize starting values

    unsigned int iter = 0;  // Iteration couter

    for (iter = 0; iter < MAX_ITER; iter++) {  // Convergence check loop
        newIm = curRe * curIm * 2 + Im;
        newRe = curRe * curRe - curIm * curIm + Re;  // Calculate new point positions

        curRe = newRe;
        curIm = newIm;

        if (curRe * curRe + curIm * curIm > DIVERGENCE_RADIUS * DIVERGENCE_RADIUS) {
            break;  // If point diverges, exit loop
        }
    }

    return colorize(iter, curRe, curIm);  // Calculate point color and return
}


void renderMandelbrotSetToArray(const double xstart, const double ystart, const double step, uint32_t* img, const unsigned int SEGMENT_HEIGHT) {
    for (unsigned int y = 0; y < SEGMENT_HEIGHT; y++) {
        for (unsigned int x = 0; x < WINDOW_WIDTH; x++) {
            img[y * WINDOW_WIDTH + x] = getPointColor(xstart + x * step, ystart + y * step); // Just calculating per-pixel colors and update image array
        }
    }
}


void renderMandelbrotSetToArraySmooth(const double xstart, const double ystart, const double step, uint32_t* img, const unsigned int SEGMENT_HEIGHT) {
    for (unsigned int y = 0; y < SEGMENT_HEIGHT; y++) {
        for (unsigned int x = 0; x < WINDOW_WIDTH; x++) {
            // Calculate colors of all points for supersampling
            uint32_t color1 = getPointColor(xstart + x * step, ystart + y * step);
            uint32_t color2 = getPointColor(xstart + x * step + step / 2, ystart + y * step);
            uint32_t color3 = getPointColor(xstart + x * step, ystart + y * step + step / 2);
            uint32_t color4 = getPointColor(xstart + x * step + step / 2, ystart + y * step + step / 2);
            
            img[y * WINDOW_WIDTH + x] = avgColor(color1, color2, color3, color4); // Calculate average color and store it
        }
    }
}

#endif

void renderMandelbrotSet(const double xstart, const double ystart, const double step, sf::Texture& renderTo) {
    uint8_t img[WINDOW_HEIGHT * WINDOW_WIDTH * 4]; // Array for image rendering
    uint32_t* image = reinterpret_cast<uint32_t*>(img); // uint8_t alias
#ifdef MULTITHREADING // If multithreading is ON,
    decltype(std::async(std::launch::async, renderMandelbrotSetToArrayAVX2, xstart, ystart, step, image, THREAD_IMAGE_PART_HEIGHT)) futures[THREADS];  // Declare array for std::future objects
    for (unsigned int i = 0; i < THREADS; i++) {

        // Lauch all threads with corresponding function
#ifdef SIMD_OPTIMIZATION
        futures[i] = std::async(std::launch::async, renderMandelbrotSetToArrayAVX2, xstart, ystart + i * THREAD_IMAGE_PART_HEIGHT * step, step, image + i * WINDOW_WIDTH * THREAD_IMAGE_PART_HEIGHT, THREAD_IMAGE_PART_HEIGHT);
#else
        futures[i] = std::async(std::launch::async, renderMandelbrotSetToArray, xstart, ystart + i * THREAD_IMAGE_PART_HEIGHT * step, step, image + i * WINDOW_WIDTH * THREAD_IMAGE_PART_HEIGHT, THREAD_IMAGE_PART_HEIGHT);
#endif
    }

    // Wait for threads to complete calculations
    for (unsigned int i = 0; i < THREADS; i++) {
        futures[i].wait();
    }

#else // If multithreading is OFF

    // Render whole frame in one large chunk
#ifdef SIMD_OPTIMIZATION
    renderMandelbrotSetToArrayAVX2(xstart, ystart, step, image, WINDOW_HEIGHT);
#else
    renderMandelbrotSetToArray(xstart, ystart, step, image, WINDOW_HEIGHT);
#endif

#endif
    // Update texture
    renderTo.update(img);
}

void renderMandelbrotSetSmooth(const double xstart, const double ystart, const double step, sf::Texture& renderTo) {
    uint8_t img[WINDOW_HEIGHT * WINDOW_WIDTH * 4]; // Declare array for pixels
    uint32_t* image = reinterpret_cast<uint32_t*>(img); // uint32_t alias

#ifdef MULTITHREADING // If multithreading is ON, 
    decltype(std::async(std::launch::async, renderMandelbrotSetToArraySmoothAVX2, xstart, ystart, step, image, THREAD_IMAGE_PART_HEIGHT)) futures[THREADS];  // Declare array for futures
    
    for (unsigned int i = 0; i < THREADS; i++) {
        // Lauch all tasks using proper function
#ifdef SIMD_OPTIMIZATION
        futures[i] = std::async(std::launch::async, renderMandelbrotSetToArraySmoothAVX2, xstart, ystart + i * THREAD_IMAGE_PART_HEIGHT * step, step, image + i * WINDOW_WIDTH * THREAD_IMAGE_PART_HEIGHT, THREAD_IMAGE_PART_HEIGHT);
#else
        futures[i] = std::async(std::launch::async, renderMandelbrotSetToArraySmooth, xstart, ystart + i * THREAD_IMAGE_PART_HEIGHT * step, step, image + i * WINDOW_WIDTH * THREAD_IMAGE_PART_HEIGHT, THREAD_IMAGE_PART_HEIGHT);
#endif
    }

    for (unsigned int i = 0; i < THREADS; i++) {
        futures[i].wait(); // Wait for all tasks to complete
    }
#else
    // Otherwise use proper function in order to render the whole image
#ifdef SIMD_OPTIMIZATION
    renderMandelbrotSetToArraySmoothAVX2(xstart, ystart, step, image, WINDOW_HEIGHT);
#else
    renderMandelbrotSetToArraySmooth(xstart, ystart, step, image, WINDOW_HEIGHT);
#endif

#endif

    renderTo.update(img);
}