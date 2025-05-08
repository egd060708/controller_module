#include "../include/osqpE_interface.h"


osqpeInterface::osqpeInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, bool _verbose)
    : mpcMatrix(_xNum, _uNum, _cNum, _eNum, _ctrlStep)
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
            //if (i <= j)
            //{
                hessian.insert(i, j) = 1e-30;
            //}
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

Matrixr osqpeInterface::_prediction(const Matrixr &y_k, const Matrixr &x_k)
{
    Matrixr result;
    result.resize(n, 1);
    result.setZero();
    // 生成预测矩阵
    this->_mpc_matrices();
    H_new = H + extraH;
    g_new = E * x_k - L * y_k - W_bar * U_pre.block(0, 0, uNum * ctrlStep, 1) + extra_g;

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
