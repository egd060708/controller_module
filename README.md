# controller_module
controller code library

该仓库提供以下代码和接口：

1. Pid模块
2. Mpc模块
   - 基于Eigen的Mpc问题转换为QP问题计算库mpcMatrix
   - 基于qpOASES求解库的接口qpOASES_interface
   - 基于quadprog++求解库的接口quadprog_interface
   - 基于tinyMpc求解库的接口tinyMpc_interface
   - 基于osqp-eigen求解库的接口osqpE_interface
