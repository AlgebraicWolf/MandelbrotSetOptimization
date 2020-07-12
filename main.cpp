#include <immintrin.h>
#include <sleef.h>
#include <x86intrin.h>

#include <SFML/Graphics.hpp>
#include <cmath>
#include <cstdio>
#include <ctime>

#include "config.h"

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

struct pointTestResult {
    unsigned int iter;
    double Re;
    double Im;
};

inline uint32_t colorize(unsigned int iter, double Re, double Im) {
    if (iter == MAX_ITER)
        return 0xFF000000;
    else {
        double internal = iter + 1 + log2(log(Re * Re + Im * Im)) / 2;

        unsigned int color1 = static_cast<unsigned int>(fabs(cos(internal)) * 0xFF);

        return 0xFF000000 | (color1 << 16) | (color1 << 8) | color1;
    }
}

inline __m128i colorizeVector(uint64_t iters_64, __m256d Res, __m256d Ims) {
    __m256d squaredRes = _mm256_mul_pd(Res, Res);
    __m256d Abs = _mm256_fmadd_pd(Ims, Ims, squaredRes);

    __m128i integer_iters = _mm_setr_epi32(iters_64 & 0xFFFF, (iters_64 >> 16) & 0xFFFF, (iters_64 >> 32) & 0xFFFF, (iters_64 >> 48) & 0xFFFF);  // Load iteration counts into vector register
    __m256d iters = _mm256_cvtepi32_pd(integer_iters);                                                                                           // Convert iteration count to double representation
    __m256d two = _mm256_set1_pd(2.0);                                                                                                           // Mask for division by two
    __m256d one = _mm256_set1_pd(1.0);                                                                                                           // Mask for subtraction of one
    __m256d max = _mm256_set1_pd(255.0);                                                                                                         // Maximal value for pixel

    __m128i resultingColor = _mm_set1_epi32(0xFF000000);  // Register for resulting color
    __m128i shuffleMask = _mm_setr_epi8(0, 0, 0, 0xFF, 4, 4, 4, 0xFF, 8, 8, 8, 0xFF, 12, 12, 12, 0xFF);

    iters = _mm256_add_pd(iters, one);

    // Calculating color channel info
    __m256d argument = _mm256_div_pd(Sleef_log2d4_u35avx2(Sleef_logd4_u35avx2(Abs)), two);
    argument = _mm256_add_pd(iters, argument);

    __m256d cos_value = Sleef_fabsd4_avx2(Sleef_cosd4_u35avx2(argument));

    __m256d color_value = _mm256_mul_pd(cos_value, max);

    __m128i integerColorValue = _mm256_cvtpd_epi32(color_value);
    integerColorValue = _mm_or_si128(_mm_shuffle_epi8(integerColorValue, shuffleMask), resultingColor);

    return integerColorValue;
}

uint32_t getPointColor(const double Re, const double Im) {
    double curRe = 0,
           curIm = 0,
           newRe = 0,
           newIm = 0;

    unsigned int iter = 0;

    for (iter = 0; iter < MAX_ITER; iter++) {
        newIm = curRe * curIm * 2 + Im;
        newRe = curRe * curRe - curIm * curIm + Re;

        curRe = newRe;
        curIm = newIm;

        if (curRe * curRe + curIm * curIm > DIVERGENCE_RADIUS * DIVERGENCE_RADIUS) {
            break;
        }
    }

    return colorize(iter, curRe, curIm);
}

__m128i getPointColorsAVX2(const double Re3, const double Re2, const double Re1, const double Re0, const double Im3, const double Im2, const double Im1, const double Im0) {
    __m256d Res = _mm256_set_pd(Re3, Re2, Re1, Re0);                                         // Initial real values of points
    __m256d Ims = _mm256_set_pd(Im3, Im2, Im1, Im0);                                         // Initial imaginary values of points
    uint64_t iters = 0;                                                                      // Initial values for iteration counter. Each 64-bit
    __m256d squaredDivergeceRadius = _mm256_set1_pd(DIVERGENCE_RADIUS * DIVERGENCE_RADIUS);  // Register used for divergence check

    __m256d curRes = _mm256_set1_pd(0);
    __m256d curIms = _mm256_set1_pd(0);

    __m256d finalRes = _mm256_set1_pd(0);
    __m256d finalIms = _mm256_set1_pd(0);

    __m256d two = _mm256_set1_pd(2);

    __m256d finalMask = _mm256_castsi256_pd(_mm256_set1_epi8(0xFF));  // Mask for elements copying
    __m256d copyMask = _mm256_castsi256_pd(_mm256_set1_epi8(0));
    __m256d extractedRes = _mm256_castsi256_pd(_mm256_set1_epi8(0));
    __m256d extractedIms = _mm256_castsi256_pd(_mm256_set1_epi8(0));

    __m256d absolute;
    __m256d mask;

    __m256d squaredRes = _mm256_mul_pd(curRes, curRes);  // L: 4, TP: 0.5

    for (int i = 0; i < MAX_ITER; i++) {
        __m256d ResIms = _mm256_mul_pd(curRes, curIms);  // L: 4, TP: 0.5

        curRes = _mm256_fmadd_pd(curRes, curRes, Res);  // L: 4, TP: 0.5

        curRes = _mm256_fnmadd_pd(curIms, curIms, curRes);  // L: 4, TP: 0.5
        curIms = _mm256_fmadd_pd(two, ResIms, Ims);         // L: 4, TP: 0.5

        squaredRes = _mm256_mul_pd(curRes, curRes);  // L: 4, TP: 0.5

        absolute = _mm256_fmadd_pd(curIms, curIms, squaredRes);  // L: 4, TP: 0.5

        mask = _mm256_cmp_pd(absolute, squaredDivergeceRadius, _CMP_LE_OS);  // L: 4, TP: 0.5

        uint64_t smallmask = _mm256_movemask_pd(mask);  // L: 2, TP: 1

        copyMask = _mm256_andnot_pd(mask, finalMask);
        extractedRes = _mm256_and_pd(curRes, copyMask);
        extractedIms = _mm256_and_pd(curIms, copyMask);
        finalRes = _mm256_or_pd(extractedRes, finalRes);
        finalIms = _mm256_or_pd(extractedIms, finalIms);

        finalMask = _mm256_and_pd(finalMask, mask);

        if (!smallmask)  // Check whether all points are diverging
            break;

        smallmask = (smallmask & 0b0001) | ((smallmask & 0b0010) << 15) | ((smallmask & 0b0100) << 30) | ((smallmask & 0b1000) << 45);  // Restructure mask
        iters += smallmask;                                                                                                             // Perform increments where required
    }

    return colorizeVector(iters, finalRes, finalIms);
}

void renderMandelbrotSetToArrayAVX2(const double xstart, const double ystart, const double step, uint32_t* img) {
    for (unsigned int y = 0; y < WINDOW_HEIGHT; y++) {
        for (unsigned int x = 0; x < WINDOW_WIDTH; x += 4) {  // Since one 256-bit AVX2 array can hold 4 pixels
            _mm_store_si128(reinterpret_cast<__m128i*>(img + y * WINDOW_WIDTH + x), getPointColorsAVX2(xstart + (x + 3) * step, xstart + (x + 2) * step, xstart + (x + 1) * step, xstart + x * step,
                                                                                                       ystart + y * step, ystart + y * step, ystart + y * step, ystart + y * step));
        }
    }
}

void renderMandelbrotSetToArray(const double xstart, const double ystart, const double step, uint32_t* img) {
    for (unsigned int y = 0; y < WINDOW_HEIGHT; y++) {
        for (unsigned int x = 0; x < WINDOW_WIDTH; x++) {
            img[y * WINDOW_WIDTH + x] = getPointColor(xstart + x * step, ystart + y * step);
        }
    }
}

void renderMandelbrotSet(const double xstart, const double ystart, const double step, sf::Texture& renderTo) {
    uint8_t img[WINDOW_HEIGHT * WINDOW_WIDTH * 4];

    renderMandelbrotSetToArrayAVX2(xstart, ystart, step, reinterpret_cast<uint32_t*>(img));
    renderTo.update(img);
}

inline uint32_t avgChanel(uint8_t c1, uint8_t c2, uint8_t c3, uint8_t c4) {
    return (static_cast<uint32_t>(c1) + c2 + c3 + c4) / 4;
}

inline uint32_t avgColor(uint32_t RGBA1, uint32_t RGBA2, uint32_t RGBA3, uint32_t RGBA4) {
    return avgChanel(RGBA1 & 0xFF, RGBA2 & 0xFF, RGBA3 & 0xFF, RGBA4 & 0xFF) |
           (avgChanel((RGBA1 >> 8) & 0xFF, (RGBA2 >> 8) & 0xFF, (RGBA3 >> 8) & 0xFF, (RGBA4 >> 8) & 0xFF) << 8) |
           (avgChanel((RGBA1 >> 16) & 0xFF, (RGBA2 >> 16) & 0xFF, (RGBA3 >> 16) & 0xFF, (RGBA4 >> 16) & 0xFF) << 16) |
           (avgChanel((RGBA1 >> 24) & 0xFF, (RGBA2 >> 24) & 0xFF, (RGBA3 >> 24) & 0xFF, (RGBA4 >> 24) & 0xFF) << 24);
}

uint32_t avgColorSSE(__m128i colors) {
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

    __m128i mergeTogether = _mm_set_epi8(0x80, 0x80, 0x80, 0x80,
                                         0x80, 0x80, 0x80, 0x80,
                                         0x80, 0x80, 0x80, 0x80,
                                         12,   8,    4,    0);

    __m128i four = _mm_set1_epi32(4);

    __m128i pixel1 = _mm_shuffle_epi8(colors, extractP1);
    __m128i pixel2 = _mm_shuffle_epi8(colors, extractP2);
    __m128i pixel3 = _mm_shuffle_epi8(colors, extractP3);
    __m128i pixel4 = _mm_shuffle_epi8(colors, extractP4);

    __m128i subsum1 = _mm_add_epi32(pixel1, pixel2);
    __m128i subsum2 = _mm_add_epi32(pixel3, pixel4);

    __m128i sum = _mm_add_epi32(subsum1, subsum2);

    __m128i expandedColor = _mm_srli_epi32(sum, 2); // Divide by four
    __m128i finalColor = _mm_shuffle_epi8(expandedColor, mergeTogether);

    return _mm_extract_epi32(finalColor, 0);
}

void renderMandelbrotSetToArraySmooth(const double xstart, const double ystart, const double step, uint32_t* img) {
    for (unsigned int y = 0; y < WINDOW_HEIGHT; y++) {
        for (unsigned int x = 0; x < WINDOW_WIDTH; x++) {
            uint32_t color1 = getPointColor(xstart + x * step, ystart + y * step);
            uint32_t color2 = getPointColor(xstart + x * step + step / 2, ystart + y * step);
            uint32_t color3 = getPointColor(xstart + x * step, ystart + y * step + step / 2);
            uint32_t color4 = getPointColor(xstart + x * step + step / 2, ystart + y * step + step / 2);
            img[y * WINDOW_WIDTH + x] = avgColor(color1, color2, color3, color4);
        }
    }
}

void renderMandelbrotSetToArraySmoothAVX2(const double xstart, const double ystart, const double step, uint32_t* img) {
    for (unsigned int y = 0; y < WINDOW_HEIGHT; y++) {
        for (unsigned int x = 0; x < WINDOW_WIDTH; x++) {
            __m128i colors = getPointColorsAVX2(xstart + x * step, xstart + x * step + step / 2, xstart + x * step, xstart + x * step + step / 2,
                                                ystart + y * step, ystart + y * step + step, ystart + y * step + step / 2, ystart + y * step + step / 2);

            img[y * WINDOW_WIDTH + x] = avgColorSSE(colors);
        }
    }
}

void renderMandelbrotSetSmooth(const double xstart, const double ystart, const double step, sf::Texture& renderTo) {
    uint8_t img[WINDOW_HEIGHT * WINDOW_WIDTH * 4];

    renderMandelbrotSetToArraySmoothAVX2(xstart, ystart, step, reinterpret_cast<uint32_t*>(img));

    uint32_t* image = reinterpret_cast<uint32_t*>(img);

    renderTo.update(img);
}