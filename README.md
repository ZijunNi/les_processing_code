# les_processing_code
主要包括对LES/fDNS/DNS程序生成的数据进行后处理的代码，以下是对各代码文件或函数的简要说明。

## part_debug_data_processing.m
*于2025.5.26建立*
本代码处理“particle_solver”中les_part_debug_mpi函数的输出结果，包含一个函数`part_data_read`，该函数以文件名为索引读取并合并les_part_debug_mpi函数的输出结果，最终得到按粒子编号排序的数据文件（矩阵的第一列为粒子编号，后续列为对应的数据）。

## part_data_info.m
*于2025.5.29建立*
本代码生成了一个结构体，其内容是“particle_solver”中保存的所有粒子变量（以`part_`开头）的相关信息，包括读取时的变量名、分量数等，最终保存到文件`part_var_info.mat`中。

## scheme_verify.m
*于2025.6.6建立*
本代码主要对LES下求解粒子中的随机变量的数值格式进行正定性验证。主要是P_22和P_33两个存在开根号的量。