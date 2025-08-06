# les_processing_code
主要包括对LES/fDNS/DNS程序生成的数据进行后处理的代码，以下是对各代码文件或函数的简要说明。
## 数据及文件夹部分
默认的处理后的数据文件放置在`processed_data`文件夹下，绘制的相关图像文件保存在`figure`文件夹下。主目录下的`grid_test.mat`文件保存有$Re_\tau = 180$的槽道湍流网格信息。

## 代码部分
### part_debug_data_processing.m
*于2025.5.26建立*
本代码处理**particle_solver**中les_part_debug_mpi函数的输出结果，包含一个函数`part_data_read`，该函数以文件名为索引读取并合并les_part_debug_mpi函数的输出结果，最终得到按粒子编号排序的数据文件（矩阵的第一列为粒子编号，后续列为对应的数据）。

### part_data_info.m
*于2025.5.29建立*
本代码生成了一个结构体，其内容是**particle_solver**中保存的所有粒子变量（以`part_`开头）的相关信息，包括读取时的变量名、分量数等，最终保存到文件`part_var_info.mat`中。

### scheme_verify.m
*于2025.6.6建立*
本代码主要对LES下求解粒子中的随机变量的数值格式进行正定性验证。主要是P_22和P_33两个存在开根号的量。

### particle_pdf_get.m 
*于2025.6.10建立*
这是一个函数文件，该函数处理particle_read.m脚本读取并保存的粒子数据文件，并根据数据文件绘制粒子PDF在垂向上的分布。

### particle_position_display.m
*于2025.6.18建立*
本代码主要对particle_read.m脚本读取并保存的粒子数据文件进行可视化处理，显示每个粒子在流场中的位置。

### lagrange_data_process.m
*于2025.7.29建立*
本代码主要处理**particle_solver**中write_data_to_file输出的文件（目前主要用于Lagrange数据处理）。代码中包括两个函数：
- `load_lagrange_data`：读取`lagrange_data`文件夹下的Lagrange粒子数据；
- `verify_lagrange_data`：对Lagrange粒子进行物理验证，验证其速度和轨迹是否物理。
- `plot_lagrange_velocity`：绘制Lagrange粒子的轨迹并用颜色表示指定速度分量的大小。对于单个粒子输入的情况，绘制的是点线连接而成的轨迹，点的颜色表示指定速度分量的大小；对于多个粒子输入的情况，绘制的所有粒子的位置分布散点图，同时颜色表征指定速度分量的大小。
- `calculate_lagrange_timescale`：输入单个粒子的时序信息，计算其三个分量的速度自相关函数以及对应的Lagrange积分时间尺度。自相关函数计算调用Matlab库函数`xcorr`实现，Lagrange积分时间尺度则采用积分到第一个零点后截断的方式进行计算。具体的介绍参见`particle_solver`项目中`doc`文件夹下的`fluid_lagrange_readme`文件。
- `calculate_mean_velocity`：计算全体流体Lagrange粒子的平均速度剖面和速度均方根剖面。