/*! @file	osqpE_interface.cpp
 *  @brief	osqp-eigen接口
 *	@author	zzr
 *  @email	2231625449@qq.com
 *  @date	2025.5
 *
 *	@usage See mpcMatrix for details
 */
#include "../include/osqpE_interface.h"

/**
 * @brief osqp-eigen接口构造
 * @param _xNum 状态维度
 * @param _uNum 输入维度
 * @param _cNum 不等式约束维度
 * @param _eNum 等式约束维度
 * @param _ctrlStep 控制周期=预测周期
 * @param _flat_mode 0为不设平滑，1为预测整体与上一次平滑，2为每一步预测平滑
 * @param _verbose 是否使能打印
 */
osqpeInterface::osqpeInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode, bool _verbose)
    : mpcMatrix(_xNum, _uNum, _cNum, _eNum, _ctrlStep, _flat_mode)
{
    n = uNum * ctrlStep;
    m = eNum * ctrlStep;
    q = cNum * ctrlStep;
    p = q + n;
    
    hessian.resize(n, n);
    gradient.resize(n);
    linearMatrix.resize(p,n);
    lowerBound.resize(p);
    upperBound.resize(p);

    linearMatrix.setZero();
    for (int i = 0; i < n; i++)
    {
        linearMatrix.insert(q + i, i) = 1;
    }

    for (int i = 0; i < H_new.rows(); i++)
    {
        for (int j = 0; j < H_new.cols(); j++) 
        {
            if (i <= j) // hessian矩阵是对称阵，可以只填充上三角
            {
                hessian.insert(i, j) = 1e-30;
            }
        }
    }

    for (int i = 0; i < cA.rows(); i++) 
    {
        for (int j = 0; j < cA.cols(); j++)
        {
            //if (i <= j)
            //{
                linearMatrix.insert(i, j) = 1e-30;
            //}
        }
    }

    hessian.makeCompressed();
    linearMatrix.makeCompressed();
    solver.settings()->setVerbosity(_verbose);
    solver.settings()->setWarmStart(true);// 默认使能热启动
    solver.data()->setNumberOfVariables(n);
    solver.data()->setNumberOfConstraints(p);
    
}

/**
 * @brief mpc预测求解
 * @param y_k 期望状态
 * @param x_k 当前轨迹
 */
Matrixr osqpeInterface::_predictionSolve(const Matrixr &y_k, const Matrixr &x_k)
{
    Matrixr result;
    result.resize(n, 1);
    result.setZero();
    // 生成预测矩阵
    this->_mpc_matrices();
    this->_update_qp(y_k, x_k);

    gradient = g_new;
    lowerBound.block(0,0,q,1) = Alb;
    lowerBound.block(q,0,n,1) = lb;
    upperBound.block(0,0,q,1) = Aub;
    upperBound.block(q,0,n,1) = ub;

    for (int k = 0; k < hessian.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(hessian, k); it; ++it)
        {
            it.valueRef() = H_new(it.row(), it.col()) + 1e-30;
        }
    }

    for (int l = 0; l < linearMatrix.outerSize(); ++l)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(linearMatrix, l); it; ++it)
        {
            if (it.row() < q)
            {
                it.valueRef() = cA(it.row(), it.col()) + 1e-30;
            }
        }
    }

    if (!solver.isInitialized()) {

        solver.data()->setHessianMatrix(hessian);
        solver.data()->setGradient(gradient);
        solver.data()->setLinearConstraintsMatrix(linearMatrix);
        solver.data()->setLowerBound(lowerBound);
        solver.data()->setUpperBound(upperBound);

        if (!solver.initSolver()) {
            std::cerr << "Failed to initialize solver!" << std::endl;
            return result;
        }
    }
    else {
        solver.updateHessianMatrix(hessian);
        solver.updateGradient(gradient);
        solver.updateLinearConstraintsMatrix(linearMatrix);
        solver.updateLowerBound(lowerBound);
        solver.updateUpperBound(upperBound);
        
    }

    if (solver.solveProblem() == OsqpEigen::ErrorExitFlag::NoError)
    {
        result = solver.getSolution();
    }

    return result;
}


// 重写预测函数
void osqpeInterface::_prediction(const Matrixr& y_k, const Matrixr& x_k)
{
    // 生成预测矩阵
    this->_mpc_matrices();
    this->_update_qp(y_k, x_k);

    gradient = g_new;
    lowerBound.block(0, 0, q, 1) = Alb;
    lowerBound.block(q, 0, n, 1) = lb;
    upperBound.block(0, 0, q, 1) = Aub;
    upperBound.block(q, 0, n, 1) = ub;

    for (int k = 0; k < hessian.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(hessian, k); it; ++it)
        {
            it.valueRef() = H_new(it.row(), it.col()) + 1e-30;
        }
    }

    for (int l = 0; l < linearMatrix.outerSize(); ++l)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(linearMatrix, l); it; ++it)
        {
            if (it.row() < q)
            {
                it.valueRef() = cA(it.row(), it.col()) + 1e-30;
            }
        }
    }
}

// 重写矩阵拷贝
void osqpeInterface::matrixCopy()
{
    if (!solver.isInitialized()) {

        solver.data()->setHessianMatrix(hessian);
        solver.data()->setGradient(gradient);
        solver.data()->setLinearConstraintsMatrix(linearMatrix);
        solver.data()->setLowerBound(lowerBound);
        solver.data()->setUpperBound(upperBound);

        if (!solver.initSolver()) {
            std::cerr << "Failed to initialize solver!" << std::endl;
            return;
        }
    }
    else {
        solver.updateHessianMatrix(hessian);
        solver.updateGradient(gradient);
        solver.updateLinearConstraintsMatrix(linearMatrix);
        solver.updateLowerBound(lowerBound);
        solver.updateUpperBound(upperBound);

    }
}

// 重写求解函数
Matrixr osqpeInterface::_solve()
{
    Matrixr result;
    result.resize(n, 1);
    result.setZero();
    if (solver.solveProblem() == OsqpEigen::ErrorExitFlag::NoError)
    {
        result = solver.getSolution();
    }

    return result;
}

