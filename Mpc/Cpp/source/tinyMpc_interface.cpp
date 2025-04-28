#include "../include/tinyMpc_interface.h"

tinympcInterface::tinympcInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, int _speed_up, int _verbose)
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
    verbose = _verbose;
    speed_up = _speed_up;

    // setup
    int status = tiny_setup(&solver, _Adyn, _Bdyn, _Q, _R, rho_value,
        xNum, uNum, ctrlStep,
        _x_min, _x_max, _u_min, _u_max, verbose);

    // Update whichever settings we'd like
    solver->settings->max_iter = this->nWSR_static;

    // Alias solver->work for brevity
    work = solver->work;
    this->isModelUpdate = 0;
}

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
    return result;
}