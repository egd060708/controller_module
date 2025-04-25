#pragma once

#include "mpcMatrix.h"
#include "TinyMPC/tiny_api.hpp"

class tinympcInterface : public mpcBase
{
public:
    tinympcInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, int _verbose);
    ~tinympcInterface() {}

    void setRegularisation(float _rho) { rho_value = _rho; }
    void setStateConstrain(const Matrixr &_xlb, const Matrixr &_xub)
    {
        _x_min = _xlb;
        _x_max = _xub;
    }

private:
    // 求解器
    TinySolver *solver;
    // 打印使能
    int verbose = 0;
    // QR矩阵正则化参数
    float rho_value = 1.0;
    // 输入矩阵
    tinyMatrix _Adyn, _Bdyn, _Q, _R, _x_min, _x_max, _u_min, _u_max;
    // 初始状态
    tinyMatrix _x0;
    // 重写预测函数
    Matrixr _prediction(const Matrixr &y_k, const Matrixr &x_k) override;
};