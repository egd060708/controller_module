/*! @file	osqpE_interface.h
 *  @brief	osqp-eigen接口
 *	@author	zzr
 *  @email	2231625449@qq.com
 *  @date	2025.5
 *
 *	@usage See mpcMatrix for details
 */
#pragma once

#include "mpcMatrix.h"
#include "OsqpEigen/OsqpEigen.h"

class osqpeInterface : public mpcMatrix
{
public:
    osqpeInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode=0, bool _verbose=false);
    ~osqpeInterface(){}
private:
    // 求解器
    OsqpEigen::Solver solver;

    // mpc矩阵
    Eigen::SparseMatrix<MPCFloat> hessian;
    Eigen::Vector<MPCFloat,-1> gradient;
    Eigen::SparseMatrix<MPCFloat> linearMatrix;
    Eigen::Vector<MPCFloat,-1> lowerBound;
    Eigen::Vector<MPCFloat,-1> upperBound;
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