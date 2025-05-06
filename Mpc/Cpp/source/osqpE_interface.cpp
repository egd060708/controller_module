#include "../include/osqpE_interface.h"


osqpeInterface::osqpeInterface(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep)
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

    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);
    solver.data()->setNumberOfVariables(n);
    solver.data()->setNumberOfConstraints(p);

    //solver.data()->setHessianMatrix(hessian);
    //solver.data()->setGradient(gradient);
    //solver.data()->setLinearConstraintsMatrix(linearMatrix);
    //solver.data()->setLowerBound(lowerBound);
    //solver.data()->setUpperBound(upperBound);
    /*solver.initSolver();*/
}

void osqpeInterface::_matrix_transfer()
{
    
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
    //std::cout << "444444444444444444444" << std::endl;
    for (int i = 0; i < H_new.rows(); i++)
        for (int j = 0; j < H_new.cols(); j++)
            if (H_new(i, j) != 0)
            {
                hessian.coeffRef(i, j) = H_new(i, j);
            }
            

    for (int i = 0; i < cA.rows(); i++)
        for (int j = 0; j < cA.cols(); j++)
            if (cA(i, j) != 0)
            {
                linearMatrix.coeffRef(i, j) = cA(i, j);
            }
            
    //std::cout << "333333333333333333333" << std::endl;
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
        //solver.updateHessianMatrix(hessian);
        solver.updateGradient(gradient);
        solver.updateLinearConstraintsMatrix(linearMatrix);
        solver.updateLowerBound(lowerBound);
        solver.updateUpperBound(upperBound);
    }
    //std::cout << "111111111111111" << std::endl;
    
    
    if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError)
    {
        result = solver.getSolution();
    }
    //std::cout << "222222222222222" << std::endl;
    return result;
}
