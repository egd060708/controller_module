/*! @file	tinyMpc_interface.cpp
 *  @brief	tinyMpc接口
 *	@author	zzr
 *  @email	2231625449@qq.com
 *  @date	2025.5
 *
 *	@usage See mpcMatrix for details
 */
#include "../include/tinyMpc_interface.h"
#include <iostream>

/**
 * @brief tinyMPC接口构造
 * @param _xNum 状态维度
 * @param _uNum 输入维度
 * @param _cNum 不等式约束维度
 * @param _eNum 等式约束维度
 * @param _ctrlStep 控制周期=预测周期
 * @param _rho 学习率
 * @param _speed_up 使用离线lqr加速模式
 * @param _verbose 是否使能打印
 */
tinympcInterface::tinympcInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, double _rho, int _speed_up, int _verbose)
    : mpcBase(_xNum, _uNum, _cNum, _eNum, _ctrlStep)
{
    this->_Adyn.resize(xNum, xNum);
    this->_Bdyn.resize(xNum, uNum);
    this->_Q.resize(xNum, xNum);
    this->_R.resize(uNum, uNum);
    this->_x_min.resize(xNum, ctrlStep);
    this->_x_max.resize(xNum, ctrlStep);
    this->_u_min.resize(uNum, ctrlStep-1);
    this->_u_max.resize(uNum, ctrlStep-1);
    this->_K.resize(uNum, xNum);
    this->_P.resize(xNum, xNum);
    this->_x_min.setConstant(-1e17);
    this->_x_max.setConstant(1e17);
    this->_u_min.setConstant(-1e17);
    this->_u_max.setConstant(1e17);
    this->_x0.resize(xNum, 1);
    rho_value = _rho;
    verbose = _verbose;
    speed_up = _speed_up;

    // setup
    int status = tiny_setup(&solver, _Adyn, _Bdyn, _Q, _R, rho_value,
        xNum, uNum, ctrlStep,
        _x_min, _x_max, _u_min, _u_max, verbose);

    // Update whichever settings we'd like
    solver->settings->max_iter = this->nWSR_static;
    solver->settings->en_state_bound = 0;

    // Alias solver->work for brevity
    work = solver->work;
    this->isModelUpdate = 0;
}

/**
 * @brief mpc预测求解
 * @param y_k 期望状态
 * @param x_k 当前轨迹
 */
Matrixr tinympcInterface::_prediction(const Matrixr &y_k, const Matrixr &x_k)
{
    
    _Adyn = A;
    _Bdyn = B;
    _Q = Q;
    _R = R;
    _K = K;
    _P = P;
    _u_min = lb.block(0,0,uNum*(ctrlStep-1),1).reshaped(uNum, ctrlStep-1);
    _u_max = ub.block(0, 0, uNum * (ctrlStep - 1), 1).reshaped(uNum, ctrlStep - 1);
    _x_min = xlb.reshaped(xNum, ctrlStep);
    _x_max = xub.reshaped(xNum, ctrlStep);

    //work->Q = (_Q + rho_value * tinyMatrix::Identity(xNum, xNum)).diagonal();
    //work->R = (_R + rho_value * tinyMatrix::Identity(uNum, uNum)).diagonal();
    work->Q = (_Q).diagonal();
    work->R = (_R).diagonal();
    work->Adyn = _Adyn;
    work->Bdyn = _Bdyn;

    work->x_min = _x_min;
    work->x_max = _x_max;
    work->u_min = _u_min;
    work->u_max = _u_max;

    if (speed_up)
    {
        solver->cache->rho = rho_value;
        solver->cache->Kinf = _K;
        solver->cache->Pinf = _P;
        tinyMatrix R1 = work->R.asDiagonal();
        solver->cache->Quu_inv = (R1 + _Bdyn.transpose() * _P * _Bdyn).inverse();
        solver->cache->AmBKt = (_Adyn - _Bdyn * _K).transpose();
    }
    else
    {
        int status = tiny_precompute_and_set_cache(solver->cache, _Adyn, _Bdyn, _Q, _R, xNum, uNum, rho_value, verbose);

    }

    // Initial state
    _x0 = X;

    // Reference trajectory
    work->Xref = Y_K.block(0,0,xNum * ctrlStep,1).reshaped(xNum, ctrlStep);

    // 1. Update measurement
    tiny_set_x0(solver, _x0);

    // 2. Solve MPC problem
    tiny_solve(solver);

    //std::cout << "iters: " << work->iter << std::endl;

    Matrixr result;
    result.resize(uNum * ctrlStep, 1);
    result.setZero();
    /*for (int i = 0; i < ctrlStep-1; i++)
    {
        for(int j = 0; j < uNum; j++)
        {
            result(j+i*ctrlStep,0) = work->u(j,i);
        }
    }*/
    result.block(0, 0, uNum, 1) = work->u.col(0);
    //std::cout << result.block(0, 0, uNum, 1) << std::endl;
    return result;
}