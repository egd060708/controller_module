/*! @file	quadprog_interface.cpp
 *  @brief	QuadProg++接口
 *	@author	zzr
 *  @email	2231625449@qq.com
 *  @date	2025.5
 *
 *	@usage See mpcMatrix for details
 */
#include "../include/quadprog_interface.h"

/**
 * @brief qp++接口构造
 * @param _xNum 状态维度
 * @param _uNum 输入维度
 * @param _cNum 不等式约束维度
 * @param _eNum 等式约束维度
 * @param _flat_mode 0为不设平滑，1为预测整体与上一次平滑，2为每一步预测平滑
 * @param _ctrlStep 控制周期=预测周期
 */
quadprogInterface::quadprogInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode)
    : mpcMatrix(_xNum, _uNum, _cNum, _eNum, _ctrlStep, _flat_mode)
{
    n = uNum * ctrlStep;
    m = eNum * ctrlStep;
    q = cNum * ctrlStep;
    p = 2 * q + 2 * n;

    _G.resize(n, n);
    _g0.resize(n);
    _CE.resize(n, m);
    _CI.resize(n, p);
    _ce0.resize(m);
    _ci0.resize(p);
    _x.resize(n);
}

/**
 * @brief mpc预测求解
 * @param y_k 期望状态
 * @param x_k 当前轨迹
 */
Vectorr quadprogInterface::_predictionSolve(const Vectorr &y_k, const Vectorr &x_k)
{
    // 生成预测矩阵
    this->_mpc_matrices();
    this->_update_qp(y_k, x_k);

    // 由于eigen库的矩阵是按列存储，因此需要手动转换为数组
    for (int i = 0; i < H_new.rows(); i++)
        for (int j = 0; j < H_new.cols(); j++)
            _G[i][j] = H_new(i, j);
    for (int i = 0; i < g_new.rows(); i++)
        for (int j = 0; j < g_new.cols(); j++)
            _g0[i * g_new.cols() + j] = g_new(i, j);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            _CE[i][j] = (cE.transpose())(i, j);
        }
    }
    for (int i = 0; i < m; i++)
    {
        _ce0[i] = Eb(i, 0);
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < q; j++)
        {
            _CI[i][j] = (cA.transpose())(i, j);
            _CI[i][q + j] = -(cA.transpose())(i, j);
        }
        for (int j = 0; j < n; j++)
        {
            _CI[i][2 * q + j] = 1;
            _CI[i][2 * q + n + j] = -1;
        }
    }
    for (int i = 0; i < q; i++)
    {
        _ci0[i] = Alb(i, 0);
        _ci0[q + i] = -Aub(i, 0);
    }
    for (int i = 0; i < n; i++)
    {
        _ci0[2 * q + i] = lb(i, 0);
        _ci0[2 * q + n + i] = -ub(i, 0);
    }
    MPCFloat value = solve_quadprog(_G, _g0, _CE, _ce0, _CI, _ci0, _x);
    Vectorr result;
    result.resize(n);
    result.setZero();
    for (int i = 0; i < n; i++)
    {
        result(i) = _x[i];
    }
    return result;
}

/**
 * @brief mpc问题预测
 * @param None
 */
void quadprogInterface::_prediction(const Vectorr& y_k, const Vectorr& x_k)
{
    // 生成预测矩阵
    this->_mpc_matrices();
    this->_update_qp(y_k, x_k);
}

/**
 * @brief 矩阵拷贝，方便加锁
 * @param None
 */
void quadprogInterface::matrixCopy()
{
    // 由于eigen库的矩阵是按列存储，因此需要手动转换为数组
    for (int i = 0; i < H_new.rows(); i++)
        for (int j = 0; j < H_new.cols(); j++)
            _G[i][j] = H_new(i, j);
    for (int i = 0; i < g_new.rows(); i++)
        for (int j = 0; j < g_new.cols(); j++)
            _g0[i * g_new.cols() + j] = g_new(i, j);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            _CE[i][j] = (cE.transpose())(i, j);
        }
    }
    for (int i = 0; i < m; i++)
    {
        _ce0[i] = Eb(i, 0);
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < q; j++)
        {
            _CI[i][j] = (cA.transpose())(i, j);
            _CI[i][q + j] = -(cA.transpose())(i, j);
        }
        for (int j = 0; j < n; j++)
        {
            _CI[i][2 * q + j] = 1;
            _CI[i][2 * q + n + j] = -1;
        }
    }
    for (int i = 0; i < q; i++)
    {
        _ci0[i] = Alb(i, 0);
        _ci0[q + i] = -Aub(i, 0);
    }
    for (int i = 0; i < n; i++)
    {
        _ci0[2 * q + i] = lb(i, 0);
        _ci0[2 * q + n + i] = -ub(i, 0);
    }
}

/**
 * @brief mpc问题求解
 * @param None
 */
Vectorr quadprogInterface::_solve()
{
    double value = solve_quadprog(_G, _g0, _CE, _ce0, _CI, _ci0, _x);
    Vectorr result;
    result.resize(n);
    result.setZero();
    for (int i = 0; i < n; i++)
    {
        result(i) = _x[i];
    }
    return result;
}