/*! @file	model_fit.h
 *  @brief	模型矩阵生成类型
 *	@author	zzr
 *  @email	2231625449@qq.com
 *  @date	2023.9.12
 *
 *	@usage
 *  	1、setFunctions
 *		2、modelGenerate
 */
#ifndef _MODEL_FIT_H_
#define _MODEL_FIT_H_

#include <Eigen/Dense>

/* 模板参数：矩阵模型的行数，矩阵模型的列数，拟合函数的阶数 */
template<int rows,int cols,int order>
class modelFit {
public:
	//Eigen::Matrix<MPCFloat, rows*cols,order+1> model;//行数区分阶次参数，列数区分矩阵位置参数
	// Eigen::Matrix<MPCFloat,-1,-1> model;
	Eigen::Matrix<MPCFloat, order + 1, rows* cols, Eigen::ColMajor> model;

	/* 递归调用求解方程结果: 自变量，方程系数，方程阶数 */
	MPCFloat functionSolve(const MPCFloat _x, const Eigen::Vector<MPCFloat,-1> _para, int _orderNum)
	{
		if (_orderNum)
		{
			return _para(_orderNum) * powf(_x, _orderNum) + functionSolve(_x, _para, _orderNum - 1);
		}
		else
		{
			return _para(_orderNum);
		}
	}
public:
	modelFit()
	{
		//model.setZero();
	}

	/* 模型生成函数	形参：自变量 */
	Eigen::Matrix<MPCFloat,-1,-1> modelGenerate(MPCFloat _x) {
		Eigen::Matrix<MPCFloat,-1,-1> result;
		result.resize(rows, cols);
		for(int i=0;i < rows;i++)
			for (int j = 0; j < cols; j++) {
				Eigen::Vector<MPCFloat, -1> tmp = model.col(i * cols + j);
				result(i, j) = functionSolve(_x,tmp, order);
			}
		return result;
	}
	/* 模型生成函数(矩阵) 形参：自变量*/
	Eigen::Matrix<MPCFloat,-1,-1> modelGenerateMat(MPCFloat _x) {
		Eigen::Matrix<MPCFloat, 1, order + 1> xMat;
		for (int i = 0; i < order + 1; i++) {
			xMat(i) = pow(_x, i);
		}
		// Eigen::Matrix<MPCFloat,-1,-1> out = xMat * model;
		// return Eigen::Map<Eigen::Matrix<MPCFloat,-1,-1>>(out.data(), cols, rows).transpose();
		//return out.reshaped(cols, rows).transpose();

		Eigen::Matrix<MPCFloat,-1,-1> result;
		result.resize(rows, cols);
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++) {
				result(i, j) = xMat*model.col(i * cols + j);
			}
		return result;
	}
	/* 设置拟合方程参数，自动按照方程个数进行循环写入 */
	void setFunctions(MPCFloat _functions[rows*cols*(order + 1)])
	{
		// 已经隐式转换为列优先
		/*model = Eigen::Map<Eigen::MatrixXd>(Eigen::Map<Eigen::Matrix<MPCFloat, rows* cols, order + 1, Eigen::RowMajor>>(_functions).data(), rows * cols, order + 1);*/
		model = Eigen::Map<Eigen::Matrix<MPCFloat, order + 1, rows* cols, Eigen::ColMajor>>(_functions);
	}

};

#endif
