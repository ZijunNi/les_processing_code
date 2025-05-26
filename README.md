# les_processing_code
主要包括对LES/fDNS/DNS程序生成的数据进行后处理的代码，以下是对各代码文件或函数的简要说明。

## part_debug_data_processing.m
*于2025.5.26建立*
本代码主要处理“particle_solver”中les_part_debug_mpi函数的输出结果，包含一个函数`part_data_read`，该函数以文件名为索引读取并合并les_part_debug_mpi函数的输出结果，最终得到按粒子编号排序的数据文件（矩阵的第一列为粒子编号，后续列为对应的数据）。