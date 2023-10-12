class ProgressBar{
    size_t length, barWidth;
    bool show, flush;
    size_t iteration;
    double progress;


public:
    ProgressBar(size_t &length, size_t &&barWidth, bool &show, bool &&flush)
    : length(length), barWidth(barWidth), show(show), flush(flush){iteration = 0;}

    ProgressBar(size_t &length, size_t &&barWidth, bool &&show, bool &&flush)
    : length(length), barWidth(barWidth), show(show), flush(flush){iteration = 0;}

    ProgressBar(size_t &&length, size_t &&barWidth, bool &show, bool &&flush)
    : length(length), barWidth(barWidth), show(show), flush(flush){iteration = 0;}

    void show_next(double value){
        iteration += 1;
        progress = (double) iteration / length;

        if (!show)
            return ;

        size_t pos = (size_t) (barWidth * progress);

        std::cout << "[";
        for (size_t i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "-";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << "iteration: " <<value << "\n";

        if (flush)
            std::cout.flush();
    }
};


