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
#include "mpcMatrixSparse.h"

class osqpeInterface : public mpcMatrix
{
public:
    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    osqpeInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode=0, bool _verbose=false);
    ~osqpeInterface(){}
    // 重写矩阵拷贝
    void matrixCopy() override;
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
    Vectorr _predictionSolve(const Vectorr &y_k, const Vectorr &x_k) override;
    // 重写预测函数
    void _prediction(const Vectorr& y_k, const Vectorr& x_k) override;
    // 重写求解函数
    Vectorr _solve() override;

};

class osqpeInterfaceSparse : public mpcMatrixSparse
{
public:
    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    osqpeInterfaceSparse(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, bool _verbose=false);
    ~osqpeInterfaceSparse(){}
    // 重写矩阵拷贝
    void matrixCopy() override;
private:
    // 求解器
    OsqpEigen::Solver solver;

    // mpc矩阵
    Eigen::SparseMatrix<MPCFloat> hessian;
    Eigen::Vector<MPCFloat,-1> gradient;
    Eigen::SparseMatrix<MPCFloat> linearMatrix;
    Eigen::Vector<MPCFloat,-1> lowerBound;
    Eigen::Vector<MPCFloat,-1> upperBound;
    int n;

    // 重写预测函数
    Vectorr _predictionSolve(const Vectorr &y_k, const Vectorr &x_k) override;
    // 重写预测函数
    void _prediction(const Vectorr& y_k, const Vectorr& x_k) override;
    // 重写求解函数
    Vectorr _solve() override;
};