本算法是按照论文<<Large-scale L-BFGS using MapReduce>>实现的，在计算lbfgs中可以采用mapreduce计算
双循环结构中的信息，可以解决海量特征的优化问题。本算法中完成了本机模拟mapreduce过程，稍后可以改成MPI版本
