#include "definitions.cpp"

class Extrapolator{
  double dx;
  size_t max_accuracy = 3;
  size_t extrapolation_order = 0;

public:
    Extrapolator(double dx, size_t extrapolation_order): dx(dx), extrapolation_order(extrapolation_order){}

  std::vector<std::vector<std::vector<double>>> derivative_array = {
        {
            {+1.0, -1.0},                                                                       // order: 1, accuracy: 1, number: 2
            {+3.0 / 2.0, -2.0, +1.0 / 2.0},                                                     // order: 1, accuracy: 2, number: 3
            {+11.0 / 6.0, -3.0, +3.0 / 2.0, -1.0 / 3.0}                                         // order: 1, accuracy: 3, number: 4
        },
        {
            {+1.0, -2.0, +1.0},                                                                 // order: 2, accuracy: 1, number: 3
            {+2.0, -5.0, +4.0, -1.0},                                                           // order: 2, accuracy: 2, number: 4
            {+35.0 / 12.0 , -26.0 / 3.0, +19.0 / 2.0, -14.0 / 3.0, 11.0 / 12.0}                 // order: 2, accuracy: 3, number: 5
        },
        {
            {+1.0, -3.0, +3.0, -1.0},                                                           // order: 3, accuracy: 1, number: 4
            {+5.0 / 2.0, -9.0, +12.0, -7.0, +3.0 / 2.0},                                        // order: 3, accuracy: 2, number: 5
            {+17.0 / 4.0,  -71.0 / 4.0, +59.0 / 2.0, -49.0 / 2.0, +41.0 / 4.0, -7.0 / 4.0}      // order: 3, accuracy: 3, number: 6
        },
        {
            {+1.0, -4.0, +6.0, -4.0, +1.0},                                                     // order: 4, accuracy: 1, number: 5
            {+3.0, -14.0, +26.0, -24.0, 11.0, -2.0},                                            // order: 4, accuracy: 2, number: 6
            {+35.0 / 6.0, -31.0, +137.0 / 2.0, -242.0 / 3.0, +107.0 / 2.0, -19.0, +17.0 / 6.0}  // order: 4, accuracy: 3, number: 7
        },
    };


    std::vector<double> get_coefficients(size_t order, size_t accuracy){
        return derivative_array[order - 1][accuracy - 1];
    }

    double get_derivative(const std::vector<double> &y_list, size_t order, size_t accuracy){
        std::vector<double> coefficients = get_coefficients(order, accuracy);

        double dy = 0;
        for (size_t i=0; i<coefficients.size(); ++i){
            dy += coefficients[i] * y_list[i];
        }

        return dy / pow(dx, order);
    }

    double get_best_derivative(const std::vector<double> &y_list, size_t order){
        size_t accuracy = get_highest_possible_accuracy(y_list, order);

        return get_derivative(y_list, order, accuracy);
    }


    size_t get_highest_possible_accuracy(const std::vector<double> &y_list, size_t order){
        size_t accuracy = y_list.size() - order;

        if (accuracy > this->max_accuracy)
            return max_accuracy;
        else
            return accuracy;
    }

    double get_factorial(size_t order){
        double factorial = 1;
        for (size_t i = 1; i < order + 1; i++){
            factorial *= i;
        }
        return factorial;
    }

    double taylor_expansion(const std::vector<double> &y_list, size_t order)
    {
        if (order == 0)
            return y_list[0];

        else{
            double derivative = this->get_best_derivative(y_list, order);
            double factorial = this->get_factorial(order);
            return derivative / factorial * pow(dx, order);
        }
    }

    double extrapolate_next(const std::vector<double> &y_list){
        double next_y = 0;
        size_t local_max_order = min(y_list.size() - 1, this->extrapolation_order);

        for (size_t order = 0; order < local_max_order + 1; ++order){
            next_y += this->taylor_expansion(y_list, order);
        }

        return next_y;
    }

};

// -