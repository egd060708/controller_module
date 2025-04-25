#pragma once
#include "mpcMatrix.h"
#include "QuadProg/QuadProg++.hh"

class quadprogInterface : public mpcMatrix
{
public:
    quadprogInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep);
    ~quadprogInterface() {}

private:
    quadprogpp::Matrix<double> _G, _CE, _CI;
    quadprogpp::Vector<double> _g0, _ce0, _ci0, _x;
    int n,m,q,p;
    // 重写预测函数
    Matrixr _prediction(const Matrixr &y_k, const Matrixr &x_k) override;
};