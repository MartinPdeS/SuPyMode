
float
ExtrapolateNext1(vector<VectorXf>& FullEigenValues, ndarray& ITRList, size_t NextIter){

  return FullEigenValues[NextIter-1][0];
}


float
ExtrapolateNext2(vector<VectorXf>& FullEigenValues, ndarray& ITRList, size_t NextIter){


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



float
ExtrapolateNext3(vector<VectorXf>& FullEigenValues, ndarray& ITRList, size_t NextIter){


  float * ITRPtr = (float*) ITRList.request().ptr;

  float NextValue;

  float y;
  float y3 = FullEigenValues[NextIter-1][0];
  float y2 = FullEigenValues[NextIter-2][0];
  float y1 = FullEigenValues[NextIter-3][0];

  float x  = ITRPtr[NextIter];
  float x3 = ITRPtr[NextIter-1];
  float x2 = ITRPtr[NextIter-2];
  float x1 = ITRPtr[NextIter-3];

  y = (y2 - y1) / (x2 - x1) * (x-x2) + y2;

  return y;
}


float
ExtrapolateNext(size_t order, vector<VectorXf>& FullEigenValues, ndarray& ITRList, size_t NextIter){

  switch(order){

    case 1: return ExtrapolateNext1(FullEigenValues, ITRList, NextIter);

    case 2: return ExtrapolateNext2(FullEigenValues, ITRList, NextIter);

    case 3: return ExtrapolateNext3(FullEigenValues, ITRList, NextIter);
  };
  return 1;
}
