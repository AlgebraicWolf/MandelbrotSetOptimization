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

    sf::Clock clock;                            // Clock for FPS measurement

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
        double internal = iter - 1 + log2(log(Re * Re + Im * Im)) / 2;

        unsigned int color1 = static_cast<unsigned int>(cos(internal) * 0xFF) % 0xFF;
        // unsigned int color2 = static_cast<unsigned int>(fabs(cos(0.5 * internal)) * 0xFF) % 0xFF;
        // unsigned int color3 = static_cast<unsigned int>(fabs(cos(1.5 * internal)) * 0xFF) % 0xFF;

        return 0xFF000000 | (color1 << 16) | (color1 << 8) | color1;
    }
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

void renderMandelbrotSetToArray(const double xstart, const double ystart, const double step, uint32_t* img) {
    for (unsigned int y = 0; y < WINDOW_HEIGHT; y++) {
        for (unsigned int x = 0; x < WINDOW_WIDTH; x++) {
            img[y * WINDOW_WIDTH + x] = getPointColor(xstart + x * step, ystart + y * step);
        }
    }
}

void renderMandelbrotSet(const double xstart, const double ystart, const double step, sf::Texture& renderTo) {
    uint8_t img[WINDOW_HEIGHT * WINDOW_WIDTH * 4];

    renderMandelbrotSetToArray(xstart, ystart, step, reinterpret_cast<uint32_t*>(img));
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

void renderMandelbrotSetSmooth(const double xstart, const double ystart, const double step, sf::Texture& renderTo) {
    uint8_t img[WINDOW_HEIGHT * WINDOW_WIDTH * 4];

    renderMandelbrotSetToArraySmooth(xstart, ystart, step, reinterpret_cast<uint32_t*>(img));

    uint32_t* image = reinterpret_cast<uint32_t*>(img);

    renderTo.update(img);
}