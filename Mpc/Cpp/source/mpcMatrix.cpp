/*! @file	mpcMatrix.cpp
 *  @brief	MPC问题矩阵封装基类
 *	@author	zzr
 *  @email	2231625449@qq.com
 *  @date	2025.5
 *
 *	@usage
 *      1、mpcInit
 *      2、setLqrFeedback(options)
 *	    3、setInputConstrain
 *      4、setIeqConstrain(options)
 *      5、setEqConstrain(options)
 *      6、mpcUpdate
 *      7、mpcSolve
 *      8、getOutput
 */
#include "../include/mpcMatrix.h"

/********************* Base *********************/

/**
 * @brief 基类构造
 * @param _xNum 状态维度
 * @param _uNum 输入维度
 * @param _cNum 不等式约束维度
 * @param _eNum 等式约束维度
 * @param _ctrlStep 控制周期=预测周期
 */
mpcBase::mpcBase(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep)
    : xNum(_xNum), uNum(_uNum), cNum(_cNum), eNum(_eNum), ctrlStep(_ctrlStep)
{
    // 离散状态空间方程
    A.resize(_xNum, _xNum);
    B.resize(_xNum, _uNum);
    // 状态权重矩阵Q，状态末端补偿矩阵R，输入权重矩阵R，平滑矩阵W
    Q.resize(_xNum, _xNum);
    F.resize(_xNum, _xNum);
    R.resize(_uNum, _uNum);
    W.resize(_uNum, _uNum);
    // LQR导出反馈增益K，P
    K.resize(_uNum, _xNum);
    P.resize(_xNum, _xNum);
    // 状态向量
    Y.resize(_xNum, 1); // 目标向量
    X.resize(_xNum, 1); // 当前状态
    U.resize(_uNum, 1); // 输出向量
    // 预测结果
    Y_K.resize(_xNum * _ctrlStep, 1); // qp求解给定
    X_K.resize(_xNum, 1);        // qp求解状态
    U_K.resize(_uNum * _ctrlStep, 1); // qp求解输出
    U_pre.resize(_uNum * _ctrlStep, 1); // 上一次qp求解的输出
    // 约束矩阵
    // 输入约束
    lb.resize(_uNum * _ctrlStep, 1);
    ub.resize(_uNum * _ctrlStep, 1);
    // 状态约束
    xlb.resize(_xNum * _ctrlStep, 1);
    xub.resize(_xNum * _ctrlStep, 1);
    // 不等式约束
    cA.resize(_cNum * _ctrlStep, _uNum * _ctrlStep);
    Alb.resize(_cNum * _ctrlStep, 1);
    Aub.resize(_cNum * _ctrlStep, 1);
    // 等式约束
    cE.resize(_eNum * _ctrlStep, _uNum * _ctrlStep);
    Eb.resize(_eNum * _ctrlStep, 1);

    A.setZero();
    B.setZero();
    Q.setIdentity();
    F.setIdentity();
    R.setZero();
    W.setZero();
    Y.setZero();
    X.setZero();
    U.setZero();
    Y_K.setZero();
    X_K.setZero();
    U_K.setZero();
    U_pre.setZero();
    lb.setConstant(-1e17);
    ub.setConstant(1e17);
    xlb.setConstant(-1e17);
    xub.setConstant(1e17);
    cA.setZero();
    Alb.setConstant(-1e17);
    Aub.setConstant(1e17);
}

/**
 * @brief 初始化
 * @param _A 设置状态空间方程A矩阵
 * @param _B 设置状态空间方程B矩阵
 * @param _Q 状态权重矩阵
 * @param _F 终端补偿权重矩阵
 * @param _R 输入权重矩阵
 * @param _W 输入平滑权重矩阵
 * @param _Ts 离散周期，若为0，则默认输入的_A,_B已经为离散状态空间方程，否则，做一阶线性化处理
 */
void mpcBase::mpcInit(const Matrixr& _A,const Matrixr& _B,const Matrixr& _Q,const Matrixr& _F,const Matrixr& _R,const Matrixr& _W, double _Ts)
{
    if (_Ts <= 0)
    {
        // 输入的是离散
        this->A = _A;
        this->B = _B;
    }
    else
    {
        // 输入的是连续，需要做离散化(使用简单的一阶离散)
        Eigen::MatrixXd AI;
        AI.resize(xNum, xNum);
        AI.setIdentity();
        this->A = AI + _Ts * _A;
        this->B = _Ts * _B;
        // TODO：做更精准的离散化
    }
    this->Q = _Q;
    this->F = _F;
    this->R = _R;
    this->W = _W;
}

/**
 * @brief 状态更新
 * @param _Y 期望状态
 * @param _X 当前状态
 * @param _nWSR 最大迭代次数
 * @param _cpu_t 最大迭代时间
 */
void mpcBase::mpcUpdate(const Matrixr& _Y,const Matrixr& _X, int _nWSR, double _cpu_t)
{
    this->Y = _Y;
    this->X = _X;
    this->nWSR_static = _nWSR;
    this->CPU_t_static = _cpu_t;
    this->X_K.block(0, 0, xNum, 1) = this->X;
    for (int i = 0; i < ctrlStep; i++)
    {
        this->Y_K.block(i * xNum, 0, xNum, 1) = this->Y;
    }
}

/**
 * @brief mpc问题求解
 * @param None
 */
void mpcBase::mpcSolve()
{
    // 执行预测
    this->U_K.block(0, 0, uNum * ctrlStep, 1) = _prediction(this->Y_K, this->X);     // 把预测输出记录下来
    this->U = this->U_K.block(0, 0, uNum, 1); // 选择第一个周期的输出作为最后输出
    this->U_pre = this->U_K; // 更新上一次求解的输出
}

/**
 * @brief 设置输入约束
 * @param _lb 输入约束下界
 * @param _ub 输入约束上界
 */
void mpcBase::setInputConstrain(const Matrixr& _lb,const Matrixr& _ub)
{
    for (int i = 0; i < ctrlStep; i++)
    {
        this->lb.block(i * uNum, 0, uNum, 1) = _lb;
        this->ub.block(i * uNum, 0, uNum, 1) = _ub;
    }
}

/**
 * @brief 设置状态约束
 * @param _lb 状态约束下界
 * @param _ub 状态约束上界
 */
void mpcBase::setStateConstrain(const Matrixr& _xlb, const Matrixr& _xub)
{
    for (int i = 0; i < ctrlStep; i++)
    {
        this->xlb.block(i * xNum, 0, xNum, 1) = _xlb;
        this->xub.block(i * xNum, 0, xNum, 1) = _xub;
    }
}

/**
 * @brief 设置不等式约束
 * @param _cA 不等式约束左乘矩阵
 * @param _Alb 不等式约束下界
 * @param _Aub 不等式约束上界
 */
void mpcBase::setIeqConstrain(const Matrixr& _cA,const Matrixr& _Alb,const Matrixr& _Aub)
{
    this->cA.setZero();
    for (int i = 0; i < ctrlStep; i++)
    {
        this->cA.block(i * cNum, i * uNum, cNum, uNum) = _cA;
        this->Alb.block(i * cNum, 0, cNum, 1) = _Alb;
        this->Aub.block(i * cNum, 0, cNum, 1) = _Aub;
    }
}

/**
 * @brief 设置等式约束
 * @param _cE 等式约束左乘矩阵
 * @param _Eb 等式约束边界
 */
void mpcBase::setEqConstrain(const Matrixr& _cE,const Matrixr& _Eb)
{
    this->cE.setZero();
    for(int i = 0; i < ctrlStep; i++)
    {
        this->cE.block(i * eNum, i * uNum, eNum, uNum) = _cE;
        this->Eb.block(i * eNum, 0, eNum, 1) = _Eb;
    }
}

/**
 * @brief 设置离线LQR反馈增益
 * @param _K 反馈增益矩阵
 * @param _P 代价矩阵
 */
void mpcBase::setLqrFeedback(const Matrixr& _K, const Matrixr& _P)
{
    this->K = _K;
    this->P = _P;
}

/**
 * @brief 手动更新上一次求解的输入
 * @param _U_pre 更新向量
 * @param _start 更新的位置
 * @param _length 更新向量的长度
 */
void mpcBase::updateLastU(const Matrixr& _U_pre, const int& _start, const int& _length)
{
    for (int i = 0; i < ctrlStep; i++)
    {
        this->U_pre.block(i * uNum + _start, 0, _length, 1) = _U_pre;
    }
}

/********************* Matrix *********************/
/**
 * @brief mpc矩阵构造
 * @param _xNum 状态维度
 * @param _uNum 输入维度
 * @param _cNum 不等式约束维度
 * @param _eNum 等式约束维度
 * @param _ctrlStep 控制周期=预测周期
 * @param _flat_mode 0为不设平滑，1为预测整体与上一次平滑，2为每一步预测平滑
 */
mpcMatrix::mpcMatrix(int _xNum, int _uNum, int _cNum, int _eNum, int _ctrlStep, uint8_t _flat_mode)
    : mpcBase(_xNum, _uNum, _cNum, _eNum, _ctrlStep)
{
    this->flat_mode = _flat_mode;
    // 中间变量
    G.resize(_xNum, _xNum);
    E.resize(_ctrlStep * _uNum, _xNum);
    L.resize(_uNum * _ctrlStep, _ctrlStep * _xNum);
    H.resize(_ctrlStep * _uNum, _ctrlStep * _uNum);
    extraH.resize(_ctrlStep * _uNum, _ctrlStep * _uNum);
    extra_g.resize(_ctrlStep * _uNum, 1);
    // 过程变量
    M.resize(_ctrlStep * _xNum, _xNum);
    C.resize(_ctrlStep * _xNum, _ctrlStep * _uNum);
    Q_bar.resize(_ctrlStep * _xNum, _ctrlStep * _xNum);
    R_bar.resize(_ctrlStep * _uNum, _ctrlStep * _uNum);
    W_bar.resize(_ctrlStep * _uNum, _ctrlStep * _uNum);
    g_new.resize(_ctrlStep * _uNum, 1);
    H_new.resize(_ctrlStep * _uNum, _ctrlStep * _uNum);
    Iup.resize(_ctrlStep * _uNum, _ctrlStep * _uNum);
    Idown.resize(_ctrlStep * _uNum, _ctrlStep * _uNum);
    Wup.resize(_ctrlStep * _uNum, _ctrlStep * _uNum);
    Wn.resize(_ctrlStep * _uNum, _ctrlStep * _uNum);

    G.setZero();
    E.setZero();
    L.setZero();
    H.setZero();
    extraH.setZero();
    extra_g.setZero();
    M.setZero();
    C.setZero();
    Q_bar.setZero();
    R_bar.setZero();
    W_bar.setZero();
    g_new.setZero();
    H_new.setZero();
    Iup.setZero();
    Idown.setZero();
    Wup.setZero();
    Wn.setZero();
}

/**
 * @brief 设置额外的代价函数
 * @param _extraH 额外代价H矩阵
 * @param _extra_g 额外代价g矩阵
 */
void mpcMatrix::setExtraCost(const Matrixr& _extraH,const Matrixr& _extra_g)
{
    for (int i = 0; i < ctrlStep; i++)
    {
        this->extraH.block(i * uNum, i * uNum, uNum, uNum) = _extraH;
        this->extra_g.block(i * uNum, 0, uNum, 1) = _extra_g;
    }
}

/**
 * @brief 初始化
 * @param _A 设置状态空间方程A矩阵
 * @param _B 设置状态空间方程B矩阵
 * @param _Q 状态权重矩阵
 * @param _F 终端补偿权重矩阵
 * @param _R 输入权重矩阵
 * @param _W 输入平滑权重矩阵
 * @param _Ts 离散周期，若为0，则默认输入的_A,_B已经为离散状态空间方程，否则，做一阶线性化处理
 */
void mpcMatrix::mpcInit(const Matrixr& _A, const Matrixr& _B, const Matrixr& _Q, const Matrixr& _F, const Matrixr& _R, const Matrixr& _W, double _Ts)
{
    if (_Ts <= 0)
    {
        // 输入的是离散
        this->A = _A;
        this->B = _B;
    }
    else
    {
        // 输入的是连续，需要做离散化(使用简单的一阶离散)
        Eigen::MatrixXd AI;
        AI.resize(xNum, xNum);
        AI.setIdentity();
        this->A = AI + _Ts * _A;
        this->B = _Ts * _B;
        // TODO：做更精准的离散化
    }
    this->Q = _Q;
    this->F = _F;
    this->R = _R;
    this->W = _W;

    // 构建kron积
    for (int i = 0; i < ctrlStep; i++)
    {
        this->Q_bar.block(i * xNum, i * xNum, xNum, xNum) = this->Q;
        this->R_bar.block(i * uNum, i * uNum, uNum, uNum) = this->R;
        if (this->flat_mode != 0)
        {
            this->W_bar.block(i * uNum, i * uNum, uNum, uNum) = this->W;
        }
        if (this->flat_mode==2 && i < (ctrlStep - 1))
        {
            this->Iup.block(i * uNum, i * uNum + uNum, uNum, uNum).setIdentity();
            this->Idown.block(i * uNum + uNum, i * uNum, uNum, uNum).setIdentity();
            this->Wup.block(i * uNum, i * uNum + uNum, uNum, uNum) = this->W;
        }
    }
    this->Q_bar.block((ctrlStep-1) * xNum, (ctrlStep-1) * xNum, xNum, xNum) = this->F;
    if (this->flat_mode == 2)
    {
        this->Wn = this->W_bar + this->Iup * this->W_bar * this->Idown - 2 * this->Wup;
    }
}

/**
 * @brief mpc矩阵生成
 * @param None
 */
void mpcMatrix::_mpc_matrices()
{
    this->M.block(0, 0, xNum, xNum).setIdentity();

    Matrixr tmp;
    tmp.resize(xNum,xNum);
    tmp.setIdentity();
    // 填充C矩阵和M矩阵
    for (int i = 1; i < ctrlStep; i++)
    {
        int rowStart = i * xNum;
        this->C.block(rowStart, 0, xNum, uNum) = tmp * B;
        if (rowStart > xNum)
        {
            this->C.block(rowStart, uNum, xNum, C.cols() - uNum) = C.block(rowStart - xNum, 0, xNum, C.cols() - uNum);
        }
        tmp = this->A * tmp;
        this->M.block(rowStart, 0, xNum, xNum) = tmp;
    }

    this->G = M.transpose() * Q_bar * M;         // G: n x n
    this->L = C.transpose() * Q_bar;             // F: NP x n
    this->E = L * M;                             // E: NP x n
    if (this->flat_mode == 0)
    {
        this->H = C.transpose() * Q_bar * C + R_bar;
    }
    else if (this->flat_mode == 1)
    {
        this->H = C.transpose() * Q_bar * C + R_bar + W_bar; // NP x NP
    }
    else if (this->flat_mode == 2)
    {
        this->H = C.transpose() * Q_bar * C + R_bar + Wn;
    }
}

/**
 * @brief 更新qp矩阵
 * @param y_k 期望状态
 * @param x_k 当前轨迹
 */
void mpcMatrix::_update_qp(const Matrixr& y_k, const Matrixr& x_k)
{
    H_new = H + extraH;
    if (this->flat_mode == 0)
    {
        g_new = E * x_k - L * y_k + extra_g;
    }
    else if (this->flat_mode == 1)
    {
        g_new = E * x_k - L * y_k - W_bar * U_pre.block(0, 0, uNum * ctrlStep, 1) + extra_g;
    }
    else if (this->flat_mode == 2)
    {
        g_new = E * x_k - L * y_k + extra_g;
        g_new.block(0, 0, uNum, 1) -= this->W * U_pre.block(0, 0, uNum, 1);
    }
}
