/*! @file	tinyMpc_interface.h
 *  @brief	tinyMpc接口
 *	@author	zzr
 *  @email	2231625449@qq.com
 *  @date	2025.5
 *
 *	@usage See mpcMatrix for details
 */
#pragma once

#include "mpcMatrix.h"
#include "TinyMPC/tiny_api.hpp"

class tinympcInterface : public mpcBase
{
public:
    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    tinympcInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, MPCFloat _rho = 1., int _speed_up = 1, int _verbose = 0);
    ~tinympcInterface() {}

    // 设置学习率
    void setRegularisation(MPCFloat _rho) { rho_value = _rho; }
    // 重写矩阵拷贝
    void matrixCopy() override;
private:
    // 求解器
    TinySolver *solver;
    TinyWorkspace *work;
    // 打印使能
    int verbose = 0;
    // QR矩阵正则化参数
    MPCFloat rho_value = 0.;
    // 是否使用离线拟合LQR矩阵进行加速
    int speed_up = 0.;
    // 输入矩阵
    tinyMatrix _Adyn, _Bdyn, _Q, _R, _x_min, _x_max, _u_min, _u_max, _K, _P;
    // 初始状态
    tinyMatrix _x0;
    // 重写预测函数
    Vectorr _predictionSolve(const Vectorr &y_k, const Vectorr &x_k) override;
    // 重写预测函数
    void _prediction(const Vectorr& y_k, const Vectorr& x_k) override;
    // 重写求解函数
    Vectorr _solve() override;
};