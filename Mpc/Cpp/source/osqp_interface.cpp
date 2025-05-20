#include "../include/osqp_interface.h"
#include <iostream>

osqpInterface::osqpInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode)
    : mpcMatrix(_xNum, _uNum, _cNum, _eNum, _ctrlStep, _flat_mode)
{
    n = uNum * ctrlStep;
    m = eNum * ctrlStep;
    q = cNum * ctrlStep;
    p = q + n;

    Hs.resize(n,n);
    As.resize(p,n);

    for (int i = 0; i < n; i++)
    {
        As.insert(q + i, i) = 1;
    }

    for (int i = 0; i < H_new.rows(); i++)
    {
        for (int j = 0; j < H_new.cols(); j++) 
        {
            //if (i <= j)
            //{
                Hs.insert(i, j) = 1e-30;
            //}
        }
    }

    for (int i = 0; i < cA.rows(); i++) 
    {
        for (int j = 0; j < cA.cols(); j++)
        {
            //if (i <= j)
            //{
                As.insert(i, j) = 1e-30;
            //}
        }
    }

    Hs.makeCompressed();
    As.makeCompressed();

    qp_q = new OSQPFloat[n];
    qp_l = new OSQPFloat[p];
    qp_u = new OSQPFloat[p];
    hsp = Hs.outerIndexPtr();
    hsi = Hs.innerIndexPtr();
    hsv = Hs.valuePtr();
    asp = As.outerIndexPtr();
    asi = As.innerIndexPtr();
    asv = As.valuePtr();
    // qp_n = static_cast<long long>(_xNum);
    // qp_m = static_cast<long long>(_uNum);
    settings = OSQPSettings_new();
    settings->alpha = 1.0;
}

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

void osqpInterface::osqpInit()
{
    this->_matrix_transfer();
     std::cout << "222222222222222222222222222222222" << std::endl;
    exitflag = osqp_setup(&solver, qp_P, qp_q, qp_A, qp_l, qp_u, p, n, settings);
    
    isSetUp = true;
}

Matrixr osqpInterface::_prediction(const Matrixr &y_k, const Matrixr &x_k)
{
    std::cout << "1111111111111111111" << std::endl;
    if (isSetUp == false)
    {
        this->osqpInit();
    }
    else
    {
        this->_matrix_transfer();
    }
    // g_new = E * x_k - L * y_k - W_bar * U_pre.block(0, 0, uNum * ctrlStep, 1) + extra_g;
    this->_update_qp(y_k,x_k);
    std::cout << "1111111111111111111" << std::endl;

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
        qp_u[q + i] = -ub(i, 0);
    }
    std::cout << "1111111111111111111" << std::endl;
    // 更新向量
    if (!exitflag)
        exitflag = osqp_update_data_vec(solver, qp_q, qp_l, qp_u);
    // 更新矩阵
    if (!exitflag)
        exitflag = osqp_update_data_mat(solver,
                                        (OSQPFloat *)hsv, OSQP_NULL, Hs.nonZeros(),
                                        (OSQPFloat *)asv, OSQP_NULL, As.nonZeros());
    std::cout << "1111111111111111111" << std::endl;
    // 求解
    if (!exitflag)
        exitflag = osqp_solve(solver);
    std::cout << "1111111111111111111" << std::endl;
    // 导出结果
    Matrixr result;
    result.resize(uNum * ctrlStep, 1);
    result.setZero();
    for (int i = 0; i < n; i++)
    {
        result(i, 0) = solver->solution->x[i];
    }
    return result;
}

void osqpInterface::_matrix_transfer()
{
    // 生成预测矩阵
    this->_mpc_matrices();
    H_new = H + extraH;
    for (int k = 0; k < Hs.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Hs, k); it; ++it)
        {
            it.valueRef() = H_new(it.row(), it.col()) + 1e-30;
        }
    }
    qp_P = OSQPCscMatrix_new(Hs.rows(), Hs.cols(), Hs.nonZeros(),
                             (OSQPFloat *)hsv, (OSQPInt *)hsi, (OSQPInt *)hsp);

    for (int l = 0; l < As.outerSize(); ++l)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(As, l); it; ++it)
        {
            if (it.row() < q)
            {
                it.valueRef() = cA(it.row(), it.col()) + 1e-30;
            }
        }
    }
    qp_A = OSQPCscMatrix_new(As.rows(), As.cols(), As.nonZeros(),
                             (OSQPFloat *)asv, (OSQPInt *)asi, (OSQPInt *)asp);
}