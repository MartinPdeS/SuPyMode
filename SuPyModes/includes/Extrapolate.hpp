
double
ExtrapolateNext1(vector<VectorType>& y, ndarray& ITRList, size_t NextIter, size_t& mode){

  return y[NextIter-1][mode];
}


double
ExtrapolateNext2(vector<VectorType>& y, ndarray& X, size_t NextIter, size_t& mode){

  double * x = (double*) X.request().ptr;

  double y0 = y[NextIter-1][0];

  double d1y = +1.0 * y[NextIter-1][mode] - 1.0 * y[NextIter-2][mode];
  double dx  = +1.0 * x[NextIter-1]    - 1.0 * x[NextIter-2];

  return y0 + d1y;

}


double
ExtrapolateNext3(vector<VectorType>& y, ndarray& X, size_t& NextIter, size_t& mode){

  double * x = (double*) X.request().ptr;

  double y0 = y[NextIter-1][0];

  double d1y = +1.0 * y[NextIter-1][mode] - 1.0 * y[NextIter-2][mode];
  double d2y = -1.0 * y[NextIter-1][mode] + 2.0 * y[NextIter-2][mode] - 1.0 * y[NextIter-3][mode];
  double dx  = +1.0 * x[NextIter-1]    - 1.0 * x[NextIter-2];

  return y0 + d1y + d2y / 2.0;
}


double
ExtrapolateNext4(vector<VectorType>& y, ndarray& X, size_t& NextIter, size_t& mode){

  double * x = (double*) X.request().ptr;

  double y0 = y[NextIter-1][0];

  double d1y = +1.0 * y[NextIter-1][mode] - 1.0 * y[NextIter-2][mode];
  double d2y = -1.0 * y[NextIter-1][mode] + 2.0 * y[NextIter-2][mode] - 1.0 * y[NextIter-3][mode];
  double d3y = +1.0 * y[NextIter-1][mode] - 3.0 * y[NextIter-2][mode] + 3.0 * y[NextIter-3][mode] - 1.0 * y[NextIter-4][mode];
  double dx  = +1.0 * x[NextIter-1]    - 1.0 * x[NextIter-2];

  return y0 + d1y + d2y / 2.0 + d3y / 6.0;
}


double
ExtrapolateNext(size_t order, vector<VectorType>& y, ndarray& X, size_t& NextIter, size_t& mode){

  switch(order){

    case 1:
         return ExtrapolateNext1(y, X, NextIter, mode);
         break;

    case 2:
         if ( NextIter < 2 )      { return ExtrapolateNext1(y, X, NextIter, mode); }
         else                     { return ExtrapolateNext2(y, X, NextIter, mode); }
         break;

    case 3:
         if      ( NextIter < 2 ) { return ExtrapolateNext1(y, X, NextIter, mode); }
         else if ( NextIter < 3 ) { return ExtrapolateNext2(y, X, NextIter, mode); }
         else                     { return ExtrapolateNext3(y, X, NextIter, mode); }
         break;

    case 4:
         if      ( NextIter < 2 ) { return ExtrapolateNext1(y, X, NextIter, mode); }
         else if ( NextIter < 3 ) { return ExtrapolateNext2(y, X, NextIter, mode); }
         else if ( NextIter < 4 ) { return ExtrapolateNext3(y, X, NextIter, mode); }
         else                     { return ExtrapolateNext4(y, X, NextIter, mode); }
         break;

  };
  return 1;
}








double
_ExtrapolateNext2(vector<VectorType>& FullEigenValues, ndarray& ITRList, size_t NextIter){


  double * ITRPtr = (double*) ITRList.request().ptr;

  double y;
  double y2 = FullEigenValues[NextIter-1][0];
  double y1 = FullEigenValues[NextIter-2][0];

  double x = ITRPtr[NextIter];
  double x2 = ITRPtr[NextIter-1];
  double x1 = ITRPtr[NextIter-2];

  y = (y2 - y1) / (x2 - x1) * (x-x2) + y2;

  return y;
}
