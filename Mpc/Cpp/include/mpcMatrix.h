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

class mpcBase
{
    public:
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
        Matrixr Y, X, U; // 目标向量、当前状态、输出向量
        // 预测结果
        Matrixr Y_K; // qp求解给定
        Matrixr X_K;              // qp求解状态
        Matrixr U_K;                      // qp求解输出
        Matrixr U_pre;                    // 上一次qp求解的输出
        // 约束矩阵
        // 输入约束
        Matrixr lb, ub;
        // 状态约束
        Matrixr xlb, xub;
        // 不等式约束
        Matrixr cA, Alb, Aub;
        // 等式约束
        Matrixr cE, Eb;
        // qp求解执行参数
        int nWSR_static = 1000; // 最大单轮qp迭代次数
        int nWSR = 1000;
        MPCFloat CPU_t_static = 1.; // 最长CPU使用时间
        MPCFloat CPU_t = 1.;
        uint8_t isModelUpdate = 1; // 系统模型是否更新

        // 初始化
        virtual void mpcInit(const Matrixr& _A,const Matrixr& _B,const Matrixr& _Q,const Matrixr& _F,const Matrixr& _R,const Matrixr& _W, MPCFloat _Ts = 0);
        // 独立更新状态空间
        virtual void setStateSpace(const Matrixr& _A,const Matrixr& _B, MPCFloat _Ts = 0);
        // 独立更新权重参数
        virtual void setWeightParams(const Matrixr& _Q,const Matrixr& _F,const Matrixr& _R,const Matrixr& _W);
        // 状态更新
        void mpcUpdate(const Matrixr& _Y,const Matrixr& _X, int _nWSR, MPCFloat _cpu_t);
        // mpc问题预测+求解
        void mpcPredictionSolve();
        // mpc问题预测
        void mpcPrediction();
        // 矩阵拷贝(prediction->solve)，用于多线程/进程
        virtual void matrixCopy() = 0;
        // mpc问题求解
        void mpcSolve();
        // 设置输入约束
        void setInputConstrain(const Matrixr& _lb,const Matrixr& _ub);
        // 设置状态约束
        void setStateConstrain(const Matrixr& _xlb, const Matrixr& _xub);
        // 设置不等式约束
        void setIeqConstrain(const Matrixr& _cA,const Matrixr& _Alb,const Matrixr& _Aub);
        // 设置等式约束
        void setEqConstrain(const Matrixr& _cE,const Matrixr& _Eb);
        // 设置离线LQR反馈增益
        void setLqrFeedback(const Matrixr& _K, const Matrixr& _P);
        // 手动更新上一次求解的输出
        void updateLastU(const Matrixr& _U_pre, const int& _start, const int& _length);
        // 获取输出值
        Matrixr getOutput(){
            return U;
        }
        // 获取预测向量
        Matrixr getPreState(){
            return X_K;
        }
        // 获取输出向量
        Matrixr getPreCtrl(){
            return U_K;
        }
    protected:
        // 执行预测并求解qp问题
        virtual Matrixr _predictionSolve(const Matrixr& y_k, const Matrixr& x_k) = 0;
        // 执行预测
        virtual void _prediction(const Matrixr& y_k, const Matrixr& x_k) = 0;
        // 执行求解
        virtual Matrixr _solve() = 0;
};

class mpcMatrix : public mpcBase
{
    public:
        // 构造函数
        mpcMatrix(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode = 0);
        // 二次型描述矩阵
        Matrixr G, E, L, H;
        Matrixr extraH, extra_g;
        // 预测矩阵
        Matrixr M, C;
        // 参数矩阵
        Matrixr Q_bar, R_bar, W_bar, g_new, H_new;
        // 用于第二类输入平滑的矩阵
        Matrixr Iup, Idown, Wup, Wn;

        // 重写权重更新
        // void mpcInit(const Matrixr& _A, const Matrixr& _B, const Matrixr& _Q, const Matrixr& _F, const Matrixr& _R, const Matrixr& _W, MPCFloat _Ts = 0) override;
        void setWeightParams(const Matrixr& _Q,const Matrixr& _F,const Matrixr& _R,const Matrixr& _W) override;
        // 设置额外代价
        void setExtraCost(const Matrixr& _extraH,const Matrixr& _extra_g);
    protected:
        // 输入平滑的方式
        uint8_t flat_mode = 0;// 0为不设平滑，1为预测整体与上一次平滑，2为每一步预测平滑
        // mpc矩阵生成
        void _mpc_matrices();
        // 更新qp矩阵
        void _update_qp(const Matrixr& y_k, const Matrixr& x_k);
};