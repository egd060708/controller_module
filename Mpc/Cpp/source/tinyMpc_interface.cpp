#include "../include/tinyMpc_interface.h"
#include <chrono>

tinympcInterface::tinympcInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, int _verbose)
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
    this->_x_min.setConstant(-1e17);
    this->_x_max.setConstant(1e17);
    this->_u_min.setConstant(-1e17);
    this->_u_max.setConstant(1e17);
    this->_x0.resize(xNum, 1);
    verbose = _verbose;

    //// setup
    //_Adyn = A;
    //_Bdyn = B;
    //_Q = Q;
    //_R = R;
    //_u_min = lb.block(0, 0, uNum * (ctrlStep - 1), 1).reshaped(uNum, ctrlStep - 1);
    //_u_max = ub.block(0, 0, uNum * (ctrlStep - 1), 1).reshaped(uNum, ctrlStep - 1);

    //int status = tiny_setup(&solver, _Adyn, _Bdyn, _Q, _R, rho_value,
    //    xNum, uNum, ctrlStep,
    //    _x_min, _x_max, _u_min, _u_max, verbose);
}

Matrixr tinympcInterface::_prediction(const Matrixr &y_k, const Matrixr &x_k)
{
    
    _Adyn = A;
    _Bdyn = B;
    _Q = Q;
    _R = R;
    _u_min = lb.block(0,0,uNum*(ctrlStep-1),1).reshaped(uNum, ctrlStep-1);
    _u_max = ub.block(0, 0, uNum * (ctrlStep - 1), 1).reshaped(uNum, ctrlStep - 1);
    
    int status = tiny_setup(&solver, _Adyn, _Bdyn, _Q, _R, rho_value,
                            xNum, uNum, ctrlStep,
                            _x_min, _x_max, _u_min, _u_max, verbose);

    // Update whichever settings we'd like
    solver->settings->max_iter = this->nWSR_static;
        
    // Alias solver->work for brevity
    TinyWorkspace *work = solver->work;

    // Initial state
    _x0 = X;

    // Reference trajectory
    work->Xref = Y_K.block(0,0,xNum * ctrlStep,1).reshaped(xNum, ctrlStep);
    // 记录开始时间
    auto start = std::chrono::high_resolution_clock::now();
    // 1. Update measurement
    tiny_set_x0(solver, _x0);

    // 2. Solve MPC problem
    tiny_solve(solver);
    // 记录结束时间
    auto end = std::chrono::high_resolution_clock::now();
    // 计算时间差（单位：毫秒）
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    std::cout << "runtime: " << duration << " us" << std::endl;
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