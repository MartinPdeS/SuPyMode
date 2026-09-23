#include "utils.h"


// template <typename T>
// void inplace_reorder_vector(std::vector<T>& vector, std::vector<size_t>& order) {
//     std::vector<T> temp(vector.size());
//     for (size_t i = 0; i < order.size(); ++i)
//         temp[i] = vector[order[i]];

//     vector.swap(temp);
// }

// template <typename T>
// std::vector<size_t> sort_indexes(const std::vector<T> &v, bool reverse) {

//   // initialize original index locations
//   std::vector<size_t> idx(v.size());
//   std::iota(idx.begin(), idx.end(), 0);

//   // sort indexes based on comparing values in v
//   // using std::stable_sort instead of std::sort
//   // to avoid unnecessary index re-orderings
//   // when v contains elements of equal values
//   std::stable_sort(idx.begin(), idx.end(),
//        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

//   if (reverse)
//       std::reverse(idx.begin(), idx.end());

//   return idx;
// }

// double dydx(double x, double y)
// {
//     return((x - y)/2);
// }

// //https://www.geeksforgeeks.org/runge-kutta-4th-order-method-solve-differential-equation/
// double rungeKutta(double x0, double y0, double x, double h)
// {
//     // Count number of iterations using step size or
//     // step height h
//     int n = (int)((x - x0) / h);

//     double k1, k2, k3, k4;

//     // Iterate for number of iterations
//     double y = y0;
//     for (int i=1; i<=n; i++)
//     {
//         // Apply Runge Kutta Formulas to find
//         // next value of y
//         k1 = h * dydx(x0, y);
//         k2 = h * dydx(x0 + 0.5*h, y + 0.5*k1);
//         k3 = h * dydx(x0 + 0.5*h, y + 0.5*k2);
//         k4 = h * dydx(x0 + h, y + k3);

//         // Update next value of y
//         y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);

//         // Update next value of x
//         x0 = x0 + h;
//     }

//     return y;
// }
