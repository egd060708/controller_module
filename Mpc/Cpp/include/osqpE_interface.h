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
    Eigen::SparseMatrix<double> hessian;
    Eigen::VectorXd gradient;
    Eigen::SparseMatrix<double> linearMatrix;
    Eigen::VectorXd lowerBound;
    Eigen::VectorXd upperBound;
    int n,m,q,p;

    // 重写预测函数
    Matrixr _prediction(const Matrixr &y_k, const Matrixr &x_k) override;
};