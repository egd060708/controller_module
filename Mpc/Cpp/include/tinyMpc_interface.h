#pragma once

#include "mpcMatrix.h"
#include "TinyMPC/tiny_api.hpp"

class tinympcInterface : public mpcBase
{
public:
    tinympcInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, double _rho, int _speed_up, int _verbose);
    ~tinympcInterface() {}

    void setRegularisation(double _rho) { rho_value = _rho; }

private:
    // 求解器
    TinySolver *solver;
    TinyWorkspace *work;
    // 打印使能
    int verbose = 0;
    // QR矩阵正则化参数
    double rho_value = 0.;
    // 是否使用离线拟合LQR矩阵进行加速
    int speed_up = 0.;
    // 输入矩阵
    tinyMatrix _Adyn, _Bdyn, _Q, _R, _x_min, _x_max, _u_min, _u_max, _K, _P;
    // 初始状态
    tinyMatrix _x0;
    // 重写预测函数
    Matrixr _prediction(const Matrixr &y_k, const Matrixr &x_k) override;
};