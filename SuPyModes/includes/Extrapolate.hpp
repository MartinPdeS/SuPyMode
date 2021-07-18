
float
ExtrapolateNext1(vector<VectorXf>& y, ndarray& ITRList, size_t NextIter){

  return y[NextIter-1][0];
}


float
ExtrapolateNext2(vector<VectorXf>& y, ndarray& X, size_t NextIter){


  float * x = (float*) X.request().ptr;

  float y0 = y[NextIter-1][0];

  float d1y = +1.0 * y[NextIter-1][0] - 1.0 * y[NextIter-2][0];
  float dx  = +1.0 * x[NextIter-1]    - 1.0 * x[NextIter-2];

  return y0 + d1y;

}


float
ExtrapolateNext3(vector<VectorXf>& y, ndarray& X, size_t NextIter){

  float * x = (float*) X.request().ptr;

  float y0 = y[NextIter-1][0];

  float d1y = +1.0 * y[NextIter-1][0] - 1.0 * y[NextIter-2][0];
  float d2y = -1.0 * y[NextIter-1][0] + 2.0 * y[NextIter-2][0] - 1.0 * y[NextIter-3][0];
  float dx  = +1.0 * x[NextIter-1]    - 1.0 * x[NextIter-2];

  return y0 + d1y + d2y / 2.0;
}


float
ExtrapolateNext4(vector<VectorXf>& y, ndarray& X, size_t NextIter){

  float * x = (float*) X.request().ptr;

  float y0 = y[NextIter-1][0];

  float d1y = +1.0 * y[NextIter-1][0] - 1.0 * y[NextIter-2][0];
  float d2y = -1.0 * y[NextIter-1][0] + 2.0 * y[NextIter-2][0] - 1.0 * y[NextIter-3][0];
  float d3y = +1.0 * y[NextIter-1][0] - 3.0 * y[NextIter-2][0] + 3.0 * y[NextIter-3][0] - 1.0 * y[NextIter-4][0];
  float dx  = +1.0 * x[NextIter-1]    - 1.0 * x[NextIter-2];

  return y0 + d1y + d2y / 2.0 + d3y / 6.0;
}


float
ExtrapolateNext(size_t order, vector<VectorXf>& y, ndarray& X, size_t NextIter){

  switch(order){

    case 1: return ExtrapolateNext1(y, X, NextIter);

    case 2:
         if ( NextIter < 2 )      { return ExtrapolateNext1(y, X, NextIter); }
         else                     { return ExtrapolateNext2(y, X, NextIter); }

    case 3:
         if      ( NextIter < 2 ) { return ExtrapolateNext1(y, X, NextIter); }
         else if ( NextIter < 3 ) { return ExtrapolateNext2(y, X, NextIter); }
         else                     { return ExtrapolateNext3(y, X, NextIter); }

    case 4:
         if      ( NextIter < 2 ) { return ExtrapolateNext1(y, X, NextIter); }
         else if ( NextIter < 3 ) { return ExtrapolateNext2(y, X, NextIter); }
         else if ( NextIter < 4 ) { return ExtrapolateNext3(y, X, NextIter); }
         else                     { return ExtrapolateNext4(y, X, NextIter); }

  };
  return 1;
}








float
_ExtrapolateNext2(vector<VectorXf>& FullEigenValues, ndarray& ITRList, size_t NextIter){


  float * ITRPtr = (float*) ITRList.request().ptr;

  float NextValue;

  float y;
  float y2 = FullEigenValues[NextIter-1][0];
  float y1 = FullEigenValues[NextIter-2][0];

  float x = ITRPtr[NextIter];
  float x2 = ITRPtr[NextIter-1];
  float x1 = ITRPtr[NextIter-2];

  y = (y2 - y1) / (x2 - x1) * (x-x2) + y2;

  return y;
}
