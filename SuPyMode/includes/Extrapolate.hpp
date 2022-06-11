

std::vector<ScalarType> CoefD1 = {+1.0, -1.0};
std::vector<ScalarType> CoefD2 = {-1.0, 2.0, -1.0};
std::vector<ScalarType> CoefD3 = {+1.0, -3.0, 3.0, -1.0};



double
ExtrapolateNext1_(vector<ScalarType>& y, std::vector<double>& x, size_t NextIter){
  return y[NextIter-1];
}


double
ExtrapolateNext2_(vector<ScalarType>& y, std::vector<double>& x, size_t NextIter){
  double d1y = CoefD1[0] * y[NextIter-1] + CoefD1[1] * y[NextIter-2];
  return ExtrapolateNext1_(y, x, NextIter) + d1y;
}


double
ExtrapolateNext3_(vector<ScalarType>& y, std::vector<double>& x, size_t& NextIter){
  double d2y = CoefD2[0] * y[NextIter-1] + CoefD2[1] * y[NextIter-2] + CoefD2[2] * y[NextIter-3];
  return ExtrapolateNext2_(y, x, NextIter) + d2y / 2.0;
}


double
ExtrapolateNext4_(vector<ScalarType>& y, std::vector<double>& x, size_t& NextIter){
  double  d3y = CoefD3[0] * y[NextIter-1] + CoefD3[1] * y[NextIter-2] + CoefD3[2] * y[NextIter-3] + CoefD3[3] * y[NextIter-4];
  return ExtrapolateNext3_(y, x, NextIter) + d3y / 6.0;
}




double
ExtrapolateNext(size_t order, vector<ScalarType>& y, std::vector<double>& x, size_t& NextIter)
{
  switch(order){

    case 1:
         return ExtrapolateNext1_(y, x, NextIter);
         break;

    case 2:
         if ( NextIter < 2 )      { return ExtrapolateNext1_(y, x, NextIter); }
         else                     { return ExtrapolateNext2_(y, x, NextIter); }
         break;

    case 3:
         if      ( NextIter < 2 ) { return ExtrapolateNext1_(y, x, NextIter); }
         else if ( NextIter < 3 ) { return ExtrapolateNext2_(y, x, NextIter); }
         else                     { return ExtrapolateNext3_(y, x, NextIter); }
         break;

    case 4:
         if      ( NextIter < 2 ) { return ExtrapolateNext1_(y, x, NextIter); }
         else if ( NextIter < 3 ) { return ExtrapolateNext2_(y, x, NextIter); }
         else if ( NextIter < 4 ) { return ExtrapolateNext3_(y, x, NextIter); }
         else                     { return ExtrapolateNext4_(y, x, NextIter); }
         break;

  };
  return 1;
}
