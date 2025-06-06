/*! @file	osqp_interface.cpp
 *  @brief	osqp接口
 *	@author	zzr
 *  @email	2231625449@qq.com
 *  @date	2025.5
 *
 *	@usage See mpcMatrix for details
 */
#include "../include/osqp_interface.h"
//#include "mpc_mat_workspace.h"
//#include <iostream>

/**
 * @brief osqp接口构造
 * @param _xNum 状态维度
 * @param _uNum 输入维度
 * @param _cNum 不等式约束维度
 * @param _eNum 等式约束维度
 * @param _ctrlStep 控制周期=预测周期
 * @param _flat_mode 0为不设平滑，1为预测整体与上一次平滑，2为每一步预测平滑
 * @param _verbose 是否使能打印
 */
osqpInterface::osqpInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode, int _verbose)
    : mpcMatrix(_xNum, _uNum, _cNum, _eNum, _ctrlStep, _flat_mode)
{
    n = uNum * ctrlStep;
    m = eNum * ctrlStep;
    q = cNum * ctrlStep;
    p = q + n;
    As.resize(p, n);
    As.setZero();
    As.block(0, 0, n, n).setIdentity();

    qp_q = new OSQPFloat[n];
    qp_l = new OSQPFloat[p];
    qp_u = new OSQPFloat[p];

    for (int i = 0; i < p; i++)
    {
        qp_l[i] = -1e17;
        qp_u[i] = 1e17;
    }

    settings = OSQPSettings_new();
    settings->verbose = _verbose;
//	solver = &mpc_mat_solver;

    this->osqpInit();
}

/**
 * @brief 析构函数
 * @param None
 */
osqpInterface::~osqpInterface()
{
    /* Cleanup */
    osqp_cleanup(solver);
    OSQPCscMatrix_free(qp_A);
    OSQPCscMatrix_free(qp_P);
    OSQPSettings_free(settings);
    delete qp_q;
    delete qp_l;
    delete qp_u;
}

/**
 * @brief osqp初始化
 * @param None
 */
void osqpInterface::osqpInit()
{
    this->_matrix_transfer();
    denseToCSC(H_new, Hvalues, Hrow_indices, Hcol_ptr, true);
    denseToCSC(As, Avalues, Arow_indices, Acol_ptr, false);
    qp_P = OSQPCscMatrix_new(H_new.rows(), H_new.cols(), Hvalues.size(), Hvalues.data(), Hrow_indices.data(), Hcol_ptr.data());
    qp_A = OSQPCscMatrix_new(As.rows(), As.cols(), Avalues.size(), Avalues.data(), Arow_indices.data(), Acol_ptr.data());
    settings->time_limit = this->CPU_t_static;
    settings->max_iter = this->nWSR_static;
    exitflag = osqp_setup(&solver, qp_P, qp_q, qp_A, qp_l, qp_u, p, n, settings);
//	solver->settings->max_iter = this->nWSR_static;
//	solver->settings->time_limit = this->CPU_t_static;
    if (!exitflag)
        isSetUp = true;
    /* Get the default codegen options */
    // OSQPCodegenDefines* defs = OSQPCodegenDefines_new();
    /* Sample with both vector and matrix updates */
//    defs->embedded_mode = 2;
//    exitflag = osqp_codegen(solver, "./osqp_generate/", "mpc_mat_", defs);
}

/**
 * @brief mpc预测求解
 * @param y_k 期望状态
 * @param x_k 当前轨迹
 */
Vectorr osqpInterface::_predictionSolve(const Vectorr &y_k, const Vectorr &x_k)
{
    if (isSetUp == false)
    {
        this->osqpInit();
    }
    else
    {
        this->_matrix_transfer();
        denseToCSCvel(H_new, Hvalues, true);
        denseToCSCvel(As, Avalues, false);
    }
    this->_update_qp(y_k,x_k);

    // 转换为osqp适用的向量
    for (int i = 0; i < g_new.rows(); i++)
        for (int j = 0; j < g_new.cols(); j++)
            qp_q[i * g_new.cols() + j] = g_new(i, j);

    for (int i = 0; i < q; i++)
    {
        qp_l[i] = Alb(i, 0);
        qp_u[i] = Aub(i, 0);
    }
    for (int i = 0; i < n; i++)
    {
        qp_l[q + i] = lb(i, 0);
        qp_u[q + i] = ub(i, 0);
    }

    // 更新向量
    if (!exitflag)
        exitflag = osqp_update_data_vec(solver, qp_q, qp_l, qp_u);
    // 更新矩阵
    if (!exitflag)
        exitflag = osqp_update_data_mat(solver,
                                        Hvalues.data(), OSQP_NULL, Hvalues.size(),
                                        Avalues.data(), OSQP_NULL, Avalues.size());

    // 求解
    if (!exitflag)
        exitflag = osqp_solve(solver);

    //std::cout << solver->info->run_time << std::endl;
    //std::cout << solver->info->iter << std::endl;
    // 导出结果
    Vectorr result;
    result.resize(uNum * ctrlStep);
    result.setZero();
    for (int i = 0; i < n; i++)
    {
        result(i) = solver->solution->x[i];
    }
    return result;
}

/**
 * @brief mpc问题预测
 * @param None
 */
void osqpInterface::_prediction(const Vectorr& y_k, const Vectorr& x_k)
{
    if (isSetUp == false)
    {
        this->osqpInit();
    }
    else
    {
        this->_matrix_transfer();
        
    }
    this->_update_qp(y_k, x_k);
}

/**
 * @brief 矩阵拷贝，方便加锁
 * @param None
 */
void osqpInterface::matrixCopy()
{
    denseToCSCvel(H_new, Hvalues, true);
    denseToCSCvel(As, Avalues, false);
    // 转换为osqp适用的向量
    for (int i = 0; i < g_new.rows(); i++)
        for (int j = 0; j < g_new.cols(); j++)
            qp_q[i * g_new.cols() + j] = g_new(i, j);

    for (int i = 0; i < q; i++)
    {
        qp_l[i] = Alb(i, 0);
        qp_u[i] = Aub(i, 0);
    }
    for (int i = 0; i < n; i++)
    {
        qp_l[q + i] = lb(i, 0);
        qp_u[q + i] = ub(i, 0);
    }
}

/**
 * @brief mpc问题求解
 * @param None
 */
Vectorr osqpInterface::_solve()
{
    // 更新向量
    if (!exitflag)
        exitflag = osqp_update_data_vec(solver, qp_q, qp_l, qp_u);
    // 更新矩阵
    if (!exitflag)
        exitflag = osqp_update_data_mat(solver,
            Hvalues.data(), OSQP_NULL, Hvalues.size(),
            Avalues.data(), OSQP_NULL, Avalues.size());

    // 求解
    if (!exitflag)
        exitflag = osqp_solve(solver);

    //std::cout << solver->info->run_time << std::endl;
    //std::cout << solver->info->iter << std::endl;
    // 导出结果
    Vectorr result;
    result.resize(uNum * ctrlStep);
    result.setZero();
    for (int i = 0; i < n; i++)
    {
        result(i) = solver->solution->x[i];
    }
    return result;
}



/**
 * @brief mpc矩阵转换成稀疏矩阵
 * @param None
 */
void osqpInterface::_matrix_transfer()
{
    // 生成预测矩阵
    this->_mpc_matrices();
    H_new = H + extraH;

    As.block(n, 0, q, n) = cA;
    // denseToCSC(H_new, Hvalues, Hrow_indices, Hcol_ptr, true);
    // denseToCSC(As, Avalues, Arow_indices, Acol_ptr, false);
}

/**
 * @brief eigen稠密矩阵转换成稀疏矩阵
 * @param None
 */
void osqpInterface::denseToCSC(
    const Matrixr& dense,
    std::vector<OSQPFloat>& values,
    std::vector<OSQPInt>& row_indices,
    std::vector<OSQPInt>& col_ptr,
    bool is_up_traingle
) {
    int rows = dense.rows();
    if (rows == 0) return;
    int cols = dense.cols();

    values.clear();
    row_indices.clear();
    col_ptr.clear();
    col_ptr.push_back(0); // col_ptr[0] = 0

    if (is_up_traingle)
    {
        for (int j = 0; j < cols; ++j) { // 遍历列
            for (int i = 0; i <= j; ++i) { // 遍历行
                values.push_back(dense(i, j));
                row_indices.push_back(i); // 0-based行索引
            }
            col_ptr.push_back(values.size()); // 记录当前列的结束位置
        }

    }
    else
    {
        for (int j = 0; j < cols; ++j) { // 遍历列
            for (int i = 0; i < rows; ++i) { // 遍历行
                values.push_back(dense(i, j));
                row_indices.push_back(i); // 0-based行索引
            }
            col_ptr.push_back(values.size()); // 记录当前列的结束位置
        }
    }

    // for (int j = 0; j < cols; ++j) { // 遍历列
    //     for (int i = 0; i < rows; ++i) { // 遍历行
    //         //if (dense[i][j] != 0.0) { // 非零元素
    //         if (is_up_traingle) // 上三角矩阵
    //         {
    //             if (j >= i)
    //             {
    //                 values.push_back(dense(i, j));
    //                 row_indices.push_back(i); // 0-based行索引
    //             }
    //         }
    //         else
    //         {
    //             values.push_back(dense(i, j));
    //             row_indices.push_back(i); // 0-based行索引
    //         }
    //         //}
    //     }
    //     col_ptr.push_back(values.size()); // 记录当前列的结束位置
    // }
}

void osqpInterface::denseToCSCvel(const Matrixr& dense, std::vector<OSQPFloat>& values, bool is_up_traingle)
{
    int rows = dense.rows();
    if (rows == 0) return;
    int cols = dense.cols();
    int index = 0;
    if (is_up_traingle)
    {
        for (int j = 0; j < cols; ++j) { // 遍历列
            for (int i = 0; i <= j; ++i) { // 遍历行
                values[index++] = dense(i, j);
            }
        }

    }
    else
    {
        for (int j = 0; j < cols; ++j) { // 遍历列
            for (int i = 0; i < rows; ++i) { // 遍历行
                values[index++] = dense(i, j);
            }
        }
    }
}