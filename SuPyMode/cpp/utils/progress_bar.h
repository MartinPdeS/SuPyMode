#pragma once

#include <vector>
#include <iostream>
#include <iomanip>


class ProgressBar {
public:
    size_t length;
    size_t barWidth;
    bool show;
    bool flush;
    size_t iteration = 0;
    double progress = 0.0;

    ProgressBar(size_t length, size_t barWidth, bool show, bool flush)
        : length(length), barWidth(barWidth), show(show), flush(flush) {}

    void show_next(double value) {
        ++iteration;
        progress = static_cast<double>(iteration) / length;

        if (!show)
            return;

        size_t pos = static_cast<size_t>(barWidth * progress);

        std::cout << "[";
        for (size_t i = 0; i < barWidth; ++i) {
            if (i < pos)
                std::cout << "-";
            else if (i == pos)
                std::cout << ">";
            else
                std::cout << " ";
        }
        std::cout << "] iteration: " << std::fixed << std::setprecision(3) << value << std::endl;

        if (flush)
            std::cout.flush();
    }
};