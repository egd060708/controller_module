/*! @file	mpcMatrix.h
 *  @brief	MPC问题矩阵封装基类
 *	@author	zzr
 *  @email	2231625449@qq.com
 *  @date	2025.5
 *
 *	@usage
 *      1、mpcInit=setStateSpace+setWeightParams
 *      2、setStateSpace(options)
 *      3、setWeightParams(options)
 *      4、setLqrFeedback(options)
 *	    5、setInputConstrain
 *      6、setIeqConstrain(options)
 *      7、setEqConstrain(options)
 *      8、mpcUpdate
 *      9、mpcSolve
 *      10、getOutput
 */
#pragma once

#include <Eigen/Dense>
#include "../config/mpcInterfaceCfg.h"

// 定义常规矩阵类型(静态)
#define MatrixSr(r, c) Eigen::Matrix<MPCFloat, r, c, Eigen::ColMajor>
// 定义方阵类型（Square）
#define MatrixSsr(d) Eigen::Matrix<MPCFloat, d, d, Eigen::ColMajor>
// 定义常规矩阵类型(动态)
#define Matrixr Eigen::Matrix<MPCFloat, -1, -1, Eigen::ColMajor>
// 定义常规向量类型(动态)
#define Vectorr Eigen::Vector<MPCFloat, -1>

class mpcBase
{
public:
    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // 基类构造
    mpcBase(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep);
    // 矩阵维度
    int xNum, uNum, cNum, eNum, ctrlStep;
    // 离散状态空间方程
    Matrixr A, B;
    // 状态权重矩阵Q，状态末端补偿矩阵R，输入权重矩阵R，平滑矩阵W
    Matrixr Q, F, R, W;
    // LQR反馈增益矩阵K，代价矩阵P
    Matrixr K, P;
    // 状态向量
    Vectorr Y, X, U; // 目标向量、当前状态、输出向量
    // 预测结果
    Vectorr Y_K;   // qp求解给定
    Vectorr X_K;   // qp求解状态
    Vectorr U_K;   // qp求解输出
    Vectorr U_pre; // 上一次qp求解的输出
    // 约束矩阵
    // 输入约束
    Vectorr lb, ub;
    // 状态约束
    Vectorr xlb, xub;
    // 不等式约束
    Matrixr cA;
    Vectorr Alb, Aub;
    // 等式约束
    Matrixr cE;
    Vectorr Eb;
    // qp求解执行参数
    int nWSR_static = 1000; // 最大单轮qp迭代次数
    int nWSR = 1000;
    MPCFloat CPU_t_static = 1.; // 最长CPU使用时间
    MPCFloat CPU_t = 1.;
    uint8_t isModelUpdate = 1; // 系统模型是否更新

    // 初始化
    virtual void mpcInit(const Matrixr &_A, const Matrixr &_B, const Vectorr &_Q, const Matrixr &_F, const Vectorr &_R, const Vectorr &_W, MPCFloat _Ts = 0);
    // 独立更新状态空间
    virtual void setStateSpace(const Matrixr &_A, const Matrixr &_B, MPCFloat _Ts = 0);
    // 独立更新权重参数
    virtual void setWeightParams(const Vectorr &_Q, const Matrixr &_F, const Vectorr &_R, const Vectorr &_W);
    // 状态更新
    virtual void mpcUpdate(const Vectorr &_Y, const Vectorr &_X, int _nWSR, MPCFloat _cpu_t);
    // mpc问题预测+求解
    virtual void mpcPredictionSolve();
    // mpc问题预测
    virtual void mpcPrediction();
    // 矩阵拷贝(prediction->solve)，用于多线程/进程
    virtual void matrixCopy() = 0;
    // mpc问题求解
    virtual void mpcSolve();
    // 设置输入约束
    virtual void setInputConstrain(const Vectorr &_lb, const Vectorr &_ub);
    // 设置状态约束
    virtual void setStateConstrain(const Vectorr &_xlb, const Vectorr &_xub);
    // 设置不等式约束
    virtual void setIeqConstrain(const Matrixr &_cA, const Vectorr &_Alb, const Vectorr &_Aub);
    // 设置等式约束
    virtual void setEqConstrain(const Matrixr &_cE, const Vectorr &_Eb);
    // 设置离线LQR反馈增益
    virtual void setLqrFeedback(const Matrixr &_K, const Matrixr &_P);
    // 手动更新上一次求解的输出
    virtual void updateLastU(const Matrixr &_U_pre, const int &_start, const int &_length);
    // 获取输出值
    Vectorr getOutput()
    {
        return U;
    }
    // 获取预测向量
    Vectorr getPreState()
    {
        return X_K;
    }
    // 获取输出向量
    Vectorr getPreCtrl()
    {
        return U_K;
    }

protected:
    // 执行预测并求解qp问题
    virtual Vectorr _predictionSolve(const Vectorr &y_k, const Vectorr &x) = 0;
    // 执行预测
    virtual void _prediction(const Vectorr &y_k, const Vectorr &x) = 0;
    // 执行求解
    virtual Vectorr _solve() = 0;
};

class mpcMatrix : public mpcBase
{
public:
    // 构造函数
    mpcMatrix(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode = 0);
    // 重写权重更新
    void setWeightParams(const Vectorr &_Q, const Matrixr &_F, const Vectorr &_R, const Vectorr &_W) override;
    // 设置额外代价
    void setExtraCost(const Matrixr &_extraH, const Matrixr &_extra_g);

protected:
    // 二次型描述矩阵
    Matrixr G, E, L, H;
    Matrixr extraH, extra_g;
    // 预测矩阵
    Matrixr M, C;
    // 参数矩阵
    Matrixr Q_bar, R_bar, W_bar, g_new, H_new;
    // 用于第二类输入平滑的矩阵
    Matrixr Iup, Idown, Wup, Wn;
    // 输入平滑的方式
    uint8_t flat_mode = 0; // 0为不设平滑，1为预测整体与上一次平滑，2为每一步预测平滑
    // mpc矩阵生成
    void _mpc_matrices();
    // 更新qp矩阵
    void _update_qp(const Vectorr &y_k, const Vectorr &x);
};