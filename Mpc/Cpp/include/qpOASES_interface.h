#pragma once

#include "mpcMatrix.h"
#include <qpOASES.hpp>
USING_NAMESPACE_QPOASES

class qpoasesInterface : public mpcMatrix
{
public:
    // 构造函数
    qpoasesInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode=0, PrintLevel _pl=PL_NONE);
    // 析构函数
    ~qpoasesInterface();

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
    Matrixr _prediction(const Matrixr &y_k, const Matrixr &x_k) override;
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
    Matrixr _prediction(const Matrixr& y_k, const Matrixr& x_k) override;
};