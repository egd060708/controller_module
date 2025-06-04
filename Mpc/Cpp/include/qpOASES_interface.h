/*! @file	qpOASES_interface.h
 *  @brief	qpOASES接口
 *	@author	zzr
 *  @email	2231625449@qq.com
 *  @date	2025.5
 *
 *	@usage See mpcMatrix for details
 */
#pragma once

#include "mpcMatrix.h"
#include <qpOASES.hpp>
USING_NAMESPACE_QPOASES

class qpoasesInterface : public mpcMatrix
{
public:
    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // 构造函数
    qpoasesInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode=0, PrintLevel _pl=PL_NONE);
    // 析构函数
    ~qpoasesInterface();
    // 重写矩阵拷贝
    void matrixCopy() override;
private:
    // qp求解器
    SQProblem qp_solver;
    Bounds guessedBounds;
    Constraints guessedConstraints;

    // 输入qpOASES的数组
    real_t *H_qpOASES;
    real_t *g_qpOASES;
    real_t *lb_qpOASES;
    real_t *ub_qpOASES;
    real_t *xOpt_qpOASES;
    real_t *yOpt_qpOASES;
    real_t *xOpt_initialGuess;
    real_t *cA_qpOASES;
    real_t *Alb_qpOASES;
    real_t *Aub_qpOASES;
    real_t *qp_out;

    // 重写预测函数
    Vectorr _predictionSolve(const Vectorr &y_k, const Vectorr &x_k) override;
    // 重写预测函数
    void _prediction(const Vectorr& y_k, const Vectorr& x_k) override;
    // 重写求解函数
    Vectorr _solve() override;
};

class qpoasesInterfaceSimple : public mpcMatrix
{
public:
    // 构造函数
    qpoasesInterfaceSimple(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode = 0, PrintLevel _pl = PL_NONE);
    // 析构函数
    ~qpoasesInterfaceSimple();

private:
    // qp求解器
    QProblemB qp_solver;
    Bounds guessedBounds;
    Constraints guessedConstraints;

    // 输入qpOASES的数组
    real_t* H_qpOASES;
    real_t* g_qpOASES;
    real_t* lb_qpOASES;
    real_t* ub_qpOASES;
    real_t* xOpt_qpOASES;
    real_t* yOpt_qpOASES;
    real_t* xOpt_initialGuess;
    real_t* qp_out;

    // 重写预测函数
    Vectorr _predictionSolve(const Vectorr& y_k, const Vectorr& x_k) override;
    // 重写预测函数
    void _prediction(const Vectorr& y_k, const Vectorr& x_k) override;
    // 重写矩阵拷贝
    void matrixCopy() override;
    // 重写求解函数
    Vectorr _solve() override;
};