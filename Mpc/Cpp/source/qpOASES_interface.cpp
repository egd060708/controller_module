/*! @file	qpOASES_interface.cpp
 *  @brief	qpOASES接口
 *	@author	zzr
 *  @email	2231625449@qq.com
 *  @date	2025.5
 *
 *	@usage See mpcMatrix for details
 */
#include "../include/qpOASES_interface.h"
 //#include <iostream>

/**
 * @brief qpoases接口构造
 * @param _xNum 状态维度
 * @param _uNum 输入维度
 * @param _cNum 不等式约束维度
 * @param _eNum 等式约束维度
 * @param _ctrlStep 控制周期=预测周期
 * @param _flat_mode 0为不设平滑，1为预测整体与上一次平滑，2为每一步预测平滑
 * @param _pl 打印等级
 */
qpoasesInterface::qpoasesInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode, PrintLevel _pl)
    : mpcMatrix(_xNum, _uNum, _cNum, _eNum, _ctrlStep, _flat_mode),
      qp_solver(_ctrlStep * _uNum, _ctrlStep * _cNum, HST_UNKNOWN)
{
    Options option;
    option.setToMPC();
    option.printLevel = _pl; // 禁用qpOASES库的打印输出
    
    //option.enableRegularisation = BT_TRUE;
    qp_solver.setOptions(option);

    // 动态分配数组内存
    H_qpOASES = new real_t[ctrlStep * uNum * ctrlStep * uNum];
    g_qpOASES = new real_t[ctrlStep * uNum];
    lb_qpOASES = new real_t[uNum * ctrlStep];
    ub_qpOASES = new real_t[uNum * ctrlStep];
    xOpt_qpOASES = new real_t[ctrlStep * uNum];
    yOpt_qpOASES = new real_t[ctrlStep * uNum + ctrlStep * cNum];
    xOpt_initialGuess = new real_t[ctrlStep * uNum];
    cA_qpOASES = new real_t[cNum * ctrlStep * uNum * ctrlStep];
    Alb_qpOASES = new real_t[cNum * ctrlStep];
    Aub_qpOASES = new real_t[cNum * ctrlStep];
    qp_out = new real_t[ctrlStep * uNum];
}

/**
 * @brief 析构函数
 * @param None
 */
qpoasesInterface::~qpoasesInterface()
{
    delete H_qpOASES;
    delete g_qpOASES;
    delete lb_qpOASES;
    delete ub_qpOASES;
    delete xOpt_qpOASES;
    delete yOpt_qpOASES;
    delete xOpt_initialGuess;
    delete cA_qpOASES;
    delete Alb_qpOASES;
    delete Aub_qpOASES;
    delete qp_out;
}

/**
 * @brief mpc预测求解
 * @param y_k 期望状态
 * @param x_k 当前轨迹
 */
Matrixr qpoasesInterface::_prediction(const Matrixr &y_k, const Matrixr &x_k)
{
    // 生成预测矩阵
    this->_mpc_matrices();
    this->_update_qp(y_k, x_k);

    // 由于eigen库的矩阵是按列存储，因此需要手动转换为数组
    //for (int i = 0; i < H_new.rows(); i++)
    //    for (int j = 0; j < H_new.cols(); j++)
    //        H_qpOASES[i * H_new.cols() + j] = H_new(i, j);
    //for (int i = 0; i < g_new.rows(); i++)
    //    for (int j = 0; j < g_new.cols(); j++)
    //        g_qpOASES[i * g_new.cols() + j] = g_new(i, j);
    //for (int i = 0; i < cA.rows(); i++)
    //    for (int j = 0; j < cA.cols(); j++)
    //        cA_qpOASES[i * cA.cols() + j] = cA(i, j);
    //for (int i = 0; i < lb.rows(); i++)
    //    for (int j = 0; j < lb.cols(); j++)
    //        lb_qpOASES[i * lb.cols() + j] = lb(i, j);
    //for (int i = 0; i < ub.rows(); i++)
    //    for (int j = 0; j < ub.cols(); j++)
    //        ub_qpOASES[i * ub.cols() + j] = ub(i, j);
    //for (int i = 0; i < Alb.rows(); i++)
    //    for (int j = 0; j < Alb.cols(); j++)
    //        Alb_qpOASES[i * Alb.cols() + j] = Alb(i, j);
    //for (int i = 0; i < Aub.rows(); i++)
    //    for (int j = 0; j < Aub.cols(); j++)
    //        Aub_qpOASES[i * Aub.cols() + j] = Aub(i, j);

    Eigen::Map<Matrixr>(H_qpOASES, H_new.rows(), H_new.cols()) = H_new;  // 直接映射，无拷贝
    Eigen::Map<Matrixr>(g_qpOASES, g_new.rows(), g_new.cols()) = g_new;
    Eigen::Map<Matrixr>(cA_qpOASES, cA.rows(), cA.cols()) = cA;
    Eigen::Map<Matrixr>(lb_qpOASES, lb.rows(), lb.cols()) = lb;
    Eigen::Map<Matrixr>(ub_qpOASES, ub.rows(), ub.cols()) = ub;
    Eigen::Map<Matrixr>(Alb_qpOASES, Alb.rows(), Alb.cols()) = Alb;
    Eigen::Map<Matrixr>(Aub_qpOASES, Aub.rows(), Aub.cols()) = Aub;

    qpOASES::returnValue ret = qpOASES::SUCCESSFUL_RETURN;
    //ret = qp_solver.init(H_qpOASES, g_qpOASES, cA_qpOASES, lb_qpOASES, ub_qpOASES, Alb_qpOASES, Aub_qpOASES, nWSR, &CPU_t);
    //std::cout << "inCPUT: " << CPU_t << std::endl;
    //std::cout << "inNWSR: " << nWSR << std::endl;
    if (isModelUpdate == 1)
    {
        // qp_solver.init(H_qpOASES, g_qpOASES, cA_qpOASES, lb_qpOASES, ub_qpOASES, Alb_qpOASES, Aub_qpOASES, nWSR, &CPU_t, xOpt_initialGuess);
        ret = qp_solver.init(H_qpOASES, g_qpOASES, cA_qpOASES, lb_qpOASES, ub_qpOASES, Alb_qpOASES, Aub_qpOASES, nWSR, &CPU_t);
        isModelUpdate = 0;
    }
    else
    {
        //ret = qp_solver.hotstart(H_qpOASES, g_qpOASES, cA_qpOASES, lb_qpOASES, ub_qpOASES, Alb_qpOASES, Aub_qpOASES, nWSR, &CPU_t, &guessedBounds, &guessedConstraints);
        ret = qp_solver.init(H_qpOASES, g_qpOASES, cA_qpOASES, lb_qpOASES, ub_qpOASES, Alb_qpOASES, Aub_qpOASES, nWSR, &CPU_t);
    }
    //std::cout << "ouCPUT: " << CPU_t << std::endl;
    //std::cout << "ouNWSR: " << nWSR << std::endl;
    nWSR = nWSR_static;
    CPU_t = CPU_t_static;
    qp_solver.getPrimalSolution(xOpt_qpOASES);
    qp_solver.getDualSolution(yOpt_qpOASES);
    //qp_solver.getBounds(guessedBounds);
    //qp_solver.getConstraints(guessedConstraints);

    Matrixr result;
    result.resize(uNum * ctrlStep, 1);
    result.setZero();
    if (ret == qpOASES::SUCCESSFUL_RETURN)
    {
        for (int i = 0; i < uNum * ctrlStep; i++)
        {
            result(i, 0) = xOpt_qpOASES[i];
        }
    }
    return result;
}





qpoasesInterfaceSimple::qpoasesInterfaceSimple(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode, PrintLevel _pl)
    : mpcMatrix(_xNum, _uNum, _cNum, _eNum, _ctrlStep, _flat_mode),
    qp_solver(_ctrlStep* _uNum, HST_UNKNOWN)
{
    Options option;
    option.setToDefault();
    option.printLevel = _pl; // 禁用qpOASES库的打印输出
    qp_solver.setOptions(option);

    // 动态分配数组内存
    H_qpOASES = new real_t[ctrlStep * uNum * ctrlStep * uNum];
    g_qpOASES = new real_t[ctrlStep * uNum];
    lb_qpOASES = new real_t[uNum * ctrlStep];
    ub_qpOASES = new real_t[uNum * ctrlStep];
    xOpt_qpOASES = new real_t[ctrlStep * uNum];
    yOpt_qpOASES = new real_t[ctrlStep * uNum + ctrlStep * cNum];
    xOpt_initialGuess = new real_t[ctrlStep * uNum];
    qp_out = new real_t[ctrlStep * uNum];
}

qpoasesInterfaceSimple::~qpoasesInterfaceSimple()
{
    delete H_qpOASES;
    delete g_qpOASES;
    delete lb_qpOASES;
    delete ub_qpOASES;
    delete xOpt_qpOASES;
    delete yOpt_qpOASES;
    delete xOpt_initialGuess;
    delete qp_out;
}

Matrixr qpoasesInterfaceSimple::_prediction(const Matrixr& y_k, const Matrixr& x_k)
{
    // 生成预测矩阵
    this->_mpc_matrices();
    g_new = E * x_k - L * y_k - W_bar * U_pre.block(0, 0, uNum * ctrlStep, 1) + extra_g;
    H_new = H + extraH;

    // 由于eigen库的矩阵是按列存储，因此需要手动转换为数组
    for (int i = 0; i < H_new.rows(); i++)
        for (int j = 0; j < H_new.cols(); j++)
            H_qpOASES[i * H_new.cols() + j] = H_new(i, j);
    for (int i = 0; i < g_new.rows(); i++)
        for (int j = 0; j < g_new.cols(); j++)
            g_qpOASES[i * g_new.cols() + j] = g_new(i, j);
    for (int i = 0; i < lb.rows(); i++)
        for (int j = 0; j < lb.cols(); j++)
            lb_qpOASES[i * lb.cols() + j] = lb(i, j);
    for (int i = 0; i < ub.rows(); i++)
        for (int j = 0; j < ub.cols(); j++)
            ub_qpOASES[i * ub.cols() + j] = ub(i, j);

    qpOASES::returnValue ret = qpOASES::SUCCESSFUL_RETURN;
    
    ret = qp_solver.init(H_qpOASES, g_qpOASES, lb_qpOASES, ub_qpOASES, nWSR, &CPU_t);

    nWSR = nWSR_static;
    CPU_t = CPU_t_static;
    qp_solver.getPrimalSolution(xOpt_qpOASES);
    qp_solver.getDualSolution(yOpt_qpOASES);
    qp_solver.getBounds(guessedBounds);

    Matrixr result;
    result.resize(uNum * ctrlStep, 1);
    result.setZero();
    if (ret == qpOASES::SUCCESSFUL_RETURN)
    {
        for (int i = 0; i < uNum * ctrlStep; i++)
        {
            result(i, 0) = xOpt_qpOASES[i];
        }
    }
    return result;
}
