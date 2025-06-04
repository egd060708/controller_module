#include "../include/mpcMatrixSparse.h"
#include <iostream>

/*************************************************** Sparse *******************************************************/

/**
 * @brief mpc矩阵稀疏构造
 * @param _xNum 状态维度
 * @param _uNum 输入维度
 * @param _cNum 不等式约束维度
 * @param _eNum 等式约束维度
 * @param _ctrlStep 控制周期=预测周期
 */
mpcMatrixSparse::mpcMatrixSparse(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep)
    : mpcBase(_xNum, _uNum, _cNum, _eNum, _ctrlStep)
{
    this->Ps.resize(xNum * (ctrlStep + 1) + uNum * ctrlStep, xNum * (ctrlStep + 1) + uNum * ctrlStep);
    this->Acs.resize(2 * xNum * (ctrlStep + 1) + uNum * ctrlStep, xNum * (ctrlStep + 1) + uNum * ctrlStep);
    this->g.resize(xNum * (ctrlStep + 1) + uNum * ctrlStep);
    this->g.setZero();
    this->l.resize(2 * xNum * (ctrlStep + 1) + uNum * ctrlStep);
    this->l.setZero();
    this->u.resize(2 * xNum * (ctrlStep + 1) + uNum * ctrlStep);
    this->u.setZero();
    for (int i = 0; i < xNum * (ctrlStep + 1) + uNum * ctrlStep; i++)
    {
        this->Ps.insert(i, i) = 1e-30;                      // 对角线元素为1
        this->Acs.insert(i + xNum * (ctrlStep + 1), i) = 1; // 对角线元素为1
    }
    for (int i = 0; i < xNum * (ctrlStep + 1); i++)
    {
        this->Acs.insert(i, i) = -1; // 对角线元素为-1
    }
    for (int i = 0; i < ctrlStep; i++)
        for (int j = 0; j < xNum; j++)
            for (int k = 0; k < xNum; k++)
            {
                this->Acs.insert(xNum * (i + 1) + j, xNum * i + k) = 1e-30;
            }
    for (int i = 0; i < ctrlStep; i++)
        for (int j = 0; j < xNum; j++)
            for (int k = 0; k < uNum; k++)
            {
                this->Acs.insert(xNum * (i + 1) + j, uNum * i + k + xNum * (ctrlStep + 1)) = 1e-30;
            }
    this->Ps.makeCompressed(); // 压缩稀疏矩阵
    this->Acs.makeCompressed(); // 压缩稀疏矩阵
}

/**
 * @brief 设置输入约束
 * @param _lb 输入约束下界
 * @param _ub 输入约束上界
 */
void mpcMatrixSparse::setInputConstrain(const Vectorr &_lb, const Vectorr &_ub)
{
    for(int i=0;i<ctrlStep;i++)
    {
        this->l.block(2*xNum*(ctrlStep+1)+uNum*i,0,uNum,1)= _lb;
        this->u.block(2*xNum*(ctrlStep+1)+uNum*i,0,uNum,1)= _ub;
    }
}

/**
 * @brief 设置状态约束
 * @param _lb 状态约束下界
 * @param _ub 状态约束上界
 */
void mpcMatrixSparse::setStateConstrain(const Vectorr &_xlb, const Vectorr &_xub)
{
    for(int i=0;i<ctrlStep+1;i++)
    {
        this->l.block(xNum*(ctrlStep+1)+xNum*i,0,xNum,1)= _xlb;
        this->u.block(xNum*(ctrlStep+1)+xNum*i,0,xNum,1)= _xub;
    }
}

/**
 * @brief 独立更新权重参数
 * @param _Q 状态权重向量
 * @param _F 终端补偿权重向量
 * @param _R 输入权重向量
 * @param _W 输入平滑权重向量
 */
void mpcMatrixSparse::setWeightParams(const Vectorr &_Q, const Matrixr &_F, const Vectorr &_R, const Vectorr &_W)
{
    this->Q = _Q.asDiagonal();
    this->F = _F;
    this->R = _R.asDiagonal();
    this->W = _W.asDiagonal();
    // 高速拼接向量并赋值
    Vectorr comb(_Q.size()*ctrlStep + _Q.size() + _R.size()*ctrlStep);
    comb.topRows(xNum * ctrlStep) = _Q.replicate(ctrlStep, 1); // 重复Q ctrlStep次
    comb.bottomRows(uNum * ctrlStep) = _R.replicate(ctrlStep, 1); // 重复R ctrlStep次
    comb.middleRows(xNum * ctrlStep, xNum) = _Q; // F只需要一次
    this->Ps.diagonal() = comb.sparseView(); // 设置对角线元素为Q,F,R，并且避免隐式转换开销
}

/**
 * @brief 更新qp矩阵
 * @param y_k 期望状态
 * @param x_k 当前轨迹
 */
void mpcMatrixSparse::_update_qp(const Vectorr &y, const Vectorr &x)
{
    this->g.topRows(xNum * ctrlStep) = (-this->Q * y).replicate(ctrlStep, 1); // 更新状态部分
    this->g.middleRows(xNum * ctrlStep, xNum) = -this->Q * y; // 更新终端补偿部分
    for (int i = 0; i < ctrlStep; i++)
        for (int j = 0; j < xNum; j++)
            for (int k = 0; k < xNum; k++)
            {
                this->Acs.coeffRef(xNum * (i + 1) + j, xNum * i + k) = A(j,k);
            }
    for (int i = 0; i < ctrlStep; i++)
        for (int j = 0; j < xNum; j++)
            for (int k = 0; k < uNum; k++)
            {
                this->Acs.coeffRef(xNum * (i + 1) + j, uNum * i + k + xNum * (ctrlStep + 1)) = B(j,k);
            }
    this->l.topRows(xNum) = -x;
    this->u.topRows(xNum) = -x;
}