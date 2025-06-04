#pragma once
#include "mpcMatrix.h"
#include "osqp.h"
#include <Eigen/Sparse>
#include <vector>

class osqpInterface : public mpcMatrix
{
public:
    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    osqpInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode=0, int _verbose=0);
    ~osqpInterface();
    void osqpInit();
    // 重写矩阵拷贝
    void matrixCopy() override;
private:
    // 转换矩阵，用于从eigen稠密阵转换成稀疏矩阵
    //Eigen::SparseMatrix<MPCFloat> Hs;
    //int* hsp = NULL;
    //int* hsi = NULL;
    //MPCFloat* hsv = NULL;
    //Eigen::SparseMatrix<MPCFloat> As;
    //int* asp = NULL;
    //int* asi = NULL;
    //MPCFloat* asv = NULL;
    int n,m,q,p;
    Matrixr As;
    std::vector<OSQPFloat> Hvalues;
    std::vector<OSQPInt> Hrow_indices, Hcol_ptr;
    std::vector<OSQPFloat> Avalues;
    std::vector<OSQPInt> Arow_indices, Acol_ptr;

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

    // 重写预测求解函数
    Vectorr _predictionSolve(const Vectorr &y_k, const Vectorr &x_k) override;
    // 重写预测函数
    void _prediction(const Vectorr& y_k, const Vectorr& x_k) override;
    // 重写求解函数
    Vectorr _solve() override;

    // 独立矩阵处理函数
    void _matrix_transfer();

    // 稠密矩阵转csc稀疏矩阵
    void denseToCSC(
        const Matrixr& dense,
        std::vector<OSQPFloat>& values,
        std::vector<OSQPInt>& row_indices,
        std::vector<OSQPInt>& col_ptr,
        bool is_up_traingle
    );

    // 稠密矩阵值转稀疏矩阵向量
    void denseToCSCvel(
        const Matrixr& dense,
        std::vector<OSQPFloat>& values,
        bool is_up_traingle
    );
};