#pragma once

#include <numeric>      // std::iota
#include <algorithm>    // std::stable_sort
#include "eigen/Eigen"


template <typename T>
size_t get_index_of_max_value(std::vector<T> &vector)
{
  return std::max_element(vector.begin(), vector.end()) - vector.begin();
}

template <typename T>
void inplace_reorder_vector(std::vector<T>& vector, std::vector<size_t>& order)
{
    assert(vector.size() == order.size());

    for( int i = 0; i < vector.size() - 1; ++i )
    {
        while( i != order[i] )
        {
            int alt = order[i];
            std::swap( vector[i], vector[alt] );
            std::swap( order[i], order[alt] );
        }
    }
}

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v, bool reverse=true) {

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

Eigen::Matrix<double, Eigen::Dynamic, 1>
GetRandomVector(size_t size){
  Eigen::Matrix<double, Eigen::Dynamic, 1> Output = Eigen::Matrix<double, Eigen::Dynamic, 1>::Random(size);
  Output.normalize();
  return Output;
}

std::vector<size_t> get_range(size_t size){
  std::vector<size_t> output;
  output.reserve(size);

  for (size_t i=0; i<size; ++i)
    output.push_back(i);

  return output;
}

double Trapz(const Eigen::Matrix<double, Eigen::Dynamic, 1> &vector, double dx, size_t n_x, size_t n_y){

    return vector.sum();

    double sum  = 0;
    double val;

    for (size_t i=0; i<n_x; ++i){
        val = vector[i];
        if ( i % n_x == 0 || i % n_x == n_x-1 )
            val /= 2.0;
        if ( i < n_y || i > n_y * (n_x-1) )
            val /= 2.0;
        sum += val;
      }

    return sum * dx;
}

double
dydx(double x, double y)
{
    return((x - y)/2);
}

//https://www.geeksforgeeks.org/runge-kutta-4th-order-method-solve-differential-equation/
double
rungeKutta(double x0, double y0, double x, double h)
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
