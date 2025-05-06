#pragma once
#include "mpcMatrix.h"
#include "osqp.h"
#include <Eigen/Sparse>

class osqpInterface : public mpcMatrix
{
public:
    osqpInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep);
    ~osqpInterface();
    void osqpInit();
private:
    // 转换矩阵，用于从eigen稠密阵转换成稀疏矩阵
    Eigen::SparseMatrix<double> Hs;
    int* hsp = NULL;
    int* hsi = NULL;
    double* hsv = NULL;
    Eigen::MatrixXd ieqcA;
    Eigen::SparseMatrix<double> As;
    int* asp = NULL;
    int* asi = NULL;
    double* asv = NULL;
    int n,m,q,p;

    /* Exitflag */
    OSQPInt exitflag = 0;
    // OSQPInt dynamic_matrices = 0;
    // OSQPInt dynamic_settings = 0;
    /* Solver, settings, matrices */
    OSQPSolver *solver = NULL;
    OSQPSettings *settings = NULL;
    // OSQPInt qp_n = 0;// 相当于状态维度
    // OSQPInt qp_m = 0;// 相当于约束维度
    OSQPFloat *qp_q = NULL;
    OSQPFloat *qp_l = NULL;
    OSQPFloat *qp_u = NULL;
    OSQPCscMatrix *qp_P = NULL;
    OSQPCscMatrix *qp_A = NULL;

    // is problem setup
    bool isSetUp = false;

    // 重写预测函数
    Matrixr _prediction(const Matrixr &y_k, const Matrixr &x_k) override;

    // 独立矩阵处理函数
    void _matrix_transfer();
};