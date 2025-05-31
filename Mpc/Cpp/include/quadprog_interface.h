/*! @file	quadprog_interface.h
 *  @brief	QuadProg++接口
 *	@author	zzr
 *  @email	2231625449@qq.com
 *  @date	2025.5
 *
 *	@usage See mpcMatrix for details
 */
#pragma once
#include "mpcMatrix.h"
#include "QuadProg/QuadProg++.hh"

class quadprogInterface : public mpcMatrix
{
public:
    quadprogInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode=0);
    ~quadprogInterface() {}

private:
    // 接口输入矩阵
    quadprogpp::qpMatrix<double> _G, _CE, _CI;
    quadprogpp::qpVector<double> _g0, _ce0, _ci0, _x;
    // 问题维度
    int n,m,q,p;
    // 重写预测函数
    Matrixr _predictionSolve(const Matrixr &y_k, const Matrixr &x_k) override;
    // 重写预测函数
    void _prediction(const Matrixr& y_k, const Matrixr& x_k) override;
    // 重写矩阵拷贝
    void matrixCopy() override;
    // 重写求解函数
    Matrixr _solve() override;
};