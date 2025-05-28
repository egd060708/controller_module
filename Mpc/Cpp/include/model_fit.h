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
	Eigen::Matrix<double, rows*cols,order+1> model;//行数区分阶次参数，列数区分矩阵位置参数

	/* 递归调用求解方程结果: 自变量，方程系数，方程阶数 */
	double functionSolve(const double _x, const Eigen::VectorXd _para, int _orderNum)
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
		model.setZero();
	}

	/* 模型生成函数	形参：自变量 */
	Eigen::MatrixXd modelGenerate(double _x) {
		Eigen::MatrixXd result;
		result.resize(rows, cols);
		for(int i=0;i < rows;i++)
			for (int j = 0; j < cols; j++) {
				Eigen::VectorXd tmp = model.row(i * cols + j);
				result(i, j) = functionSolve(_x,tmp, order);
			}
		return result;
	}

	/* 模型生成函数(矩阵) 形参：自变量*/
	Eigen::MatrixXd modelGenerateMat(double _x){
		Eigen::Matrix<double, order+1, rows*cols> xMat;
		for (int i = 0; i < order + 1; i++) {
			xMat.row(i) = Eigen::VectorXd::Constant(rows * cols, powf(_x, i));
		}
		return (model * xMat).view(rows, cols);
	}

	/* 设置拟合方程参数，自动按照方程个数进行循环写入 */
	void setFunctions(double _functions[rows*cols*(order + 1)])
	{
		model = Eigen::Map<Eigen::Matrix<double, rows* cols, order + 1, Eigen::RowMajor>>(_functions);
	}

};

#endif
