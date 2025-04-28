#pragma once

#include <Eigen/Dense>

// 定义常规矩阵类型(静态)
#define MatrixSr(r, c) Eigen::Matrix<double, r, c>
// 定义方阵类型（Square）
#define MatrixSsr(d) Eigen::Matrix<double, d, d>
// 定义常规矩阵类型(动态)
#define Matrixr Eigen::MatrixXd

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
        // 不等式约束
        Matrixr cA, Alb, Aub;
        // 等式约束
        Matrixr cE, Eb;
        // qp求解执行参数
        int nWSR_static = 1000; // 最大单轮qp迭代次数
        int nWSR = 1000;
        double CPU_t_static = 1.; // 最长CPU使用时间
        double CPU_t = 1.;
        uint8_t isModelUpdate = 1; // 系统模型是否更新

        // 初始化
        void mpcInit(const Matrixr& _A,const Matrixr& _B,const Matrixr& _Q,const Matrixr& _F,const Matrixr& _R,const Matrixr& _W, double _Ts);
        // 状态更新
        void mpcUpdate(const Matrixr& _Y,const Matrixr& _X, int _nWSR, double _cpu_t);
        // mpc问题求解
        void mpcSolve();
        // 设置输入约束
        void setConstrain(const Matrixr& _lb,const Matrixr& _ub);
        // 设置不等式约束
        void setIeqConstrain(const Matrixr& _cA,const Matrixr& _Alb,const Matrixr& _Aub);
        // 设置等式约束
        void setEqConstrain(const Matrixr& _cE,const Matrixr& _Eb);
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
        virtual Matrixr _prediction(const Matrixr& y_k, const Matrixr& x_k) = 0;
};

class mpcMatrix : public mpcBase
{
    public:
        // 构造函数
        mpcMatrix(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep);
        // 二次型描述矩阵
        Matrixr G, E, L, H;
        Matrixr extraH, extra_g;
        // 预测矩阵
        Matrixr M, C;
        Matrixr Q_bar, R_bar, W_bar, g_new, H_new;

        // 设置额外代价
        void setExtraCost(const Matrixr& _extraH,const Matrixr& _extra_g);
    protected:
        // mpc矩阵生成
        void _mpc_matrices();

};