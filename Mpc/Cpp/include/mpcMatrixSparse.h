#pragma once
#include "mpcMatrix.h"
#include <Eigen/Sparse>

class mpcMatrixSparse : public mpcBase
{
public:
    // 构造函数
    mpcMatrixSparse(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep);
    // 重写输入约束配置
    void setInputConstrain(const Vectorr &_lb, const Vectorr &_ub) override;
    // 重写状态约束配置
    void setStateConstrain(const Vectorr &_xlb, const Vectorr &_xub) override;
    // 重写权重更新
    void setWeightParams(const Vectorr &_Q, const Matrixr &_F, const Vectorr &_R, const Vectorr &_W) override;

protected:
    // 构建稀疏形式的mpc矩阵
    Eigen::SparseMatrix<MPCFloat> Ps, Acs;
    Vectorr g, u, l;
    // 更新qp矩阵
    void _update_qp(const Vectorr &y_k, const Vectorr &x);
};