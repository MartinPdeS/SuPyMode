#pragma once

#include <numeric>      // std::iota
#include <algorithm>    // std::stable_sort
#include <Eigen/Core>
#include <vector>
#include <iostream>
#include <iomanip>

template <typename T>
inline void inplace_reorder_vector(std::vector<T>& vector, std::vector<size_t>& order) {
    std::vector<T> temp(vector.size());
    for (size_t i = 0; i < order.size(); ++i)
        temp[i] = vector[order[i]];

    vector.swap(temp);
}

template <typename T>
inline std::vector<size_t> sort_indexes(const std::vector<T> &v, bool reverse = true) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  if (reverse)
      std::reverse(idx.begin(), idx.end());

  return idx;
}

inline double dydx(double x, double y)
{
    return((x - y)/2);
}

//https://www.geeksforgeeks.org/runge-kutta-4th-order-method-solve-differential-equation/
inline double rungeKutta(double x0, double y0, double x, double h)
{
    // Count number of iterations using step size or
    // step height h
    int n = (int)((x - x0) / h);

    double k1, k2, k3, k4;

    // Iterate for number of iterations
    double y = y0;
    for (int i=1; i<=n; i++)
    {
        // Apply Runge Kutta Formulas to find
        // next value of y
        k1 = h * dydx(x0, y);
        k2 = h * dydx(x0 + 0.5*h, y + 0.5*k1);
        k3 = h * dydx(x0 + 0.5*h, y + 0.5*k2);
        k4 = h * dydx(x0 + h, y + k3);

        // Update next value of y
        y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);

        // Update next value of x
        x0 = x0 + h;
    }

    return y;
}

class Extrapolator {
    public:
        double dx;
        size_t max_accuracy = 3;
        size_t extrapolation_order = 0;

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


        std::vector<double> get_coefficients(size_t order, size_t accuracy) const {
            return derivative_array[order - 1][accuracy - 1];
        }

        double get_derivative(const std::vector<double> &y_list, size_t order, size_t accuracy) const {
            std::vector<double> coefficients = get_coefficients(order, accuracy);

            double dy = 0;
            for (size_t i=0; i<coefficients.size(); ++i){
                dy += coefficients[i] * y_list[i];
            }

            return dy / pow(dx, order);
        }

        double get_best_derivative(const std::vector<double> &y_list, size_t order) const {
            size_t accuracy = get_highest_possible_accuracy(y_list, order);

            return get_derivative(y_list, order, accuracy);
        }


        size_t get_highest_possible_accuracy(const std::vector<double> &y_list, size_t order) const {
            size_t accuracy = y_list.size() - order;

            if (accuracy > this->max_accuracy)
                return max_accuracy;
            else
                return accuracy;
        }

        double get_factorial(size_t order) const {
            double factorial = 1;
            for (size_t i = 1; i < order + 1; i++){
                factorial *= i;
            }
            return factorial;
        }

        double taylor_expansion(const std::vector<double> &y_list, size_t order) const {
            if (order == 0)
                return y_list[0];

            else{
                double derivative = this->get_best_derivative(y_list, order);
                double factorial = this->get_factorial(order);
                return derivative / factorial * pow(dx, order);
            }
        }

        double extrapolate_next(const std::vector<double> &y_list) const {
            double next_y = 0;
            size_t local_max_order = std::min(y_list.size() - 1, this->extrapolation_order);

            for (size_t order = 0; order < local_max_order + 1; ++order){
                next_y += this->taylor_expansion(y_list, order);
            }

            return next_y;
        }

};



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