% 本代码生成了一个结构体，保存所有粒子变量的相关信息

%%%%%%partticle_solver程序中的变量表%%%%%
    % real,save,allocatable:: part_vel_relative   (:,:)!uf-up，共四个分量，第0个分量存储相对速度大小，
    % real,save,allocatable:: part_tau_p          (:,:)
    % real,save,allocatable:: part_tau_r          (:,:)
    % real,save,allocatable:: part_k_sgs          (:,:)
    % real,save,allocatable:: part_sgs_dissipation(:,:)!epsilon
    % real,save,allocatable:: part_t_sgs          (:,:)
    % real,save,allocatable:: part_t_l            (:,:)
    % real,save,allocatable:: part_b_i            (:,:)
    % real,save,allocatable:: part_c_i            (:,:)
    % real,save,allocatable:: part_c_star         (:,:)
    % real,save,allocatable:: part_b_ii           (:,:)
%%%%%变量表结束%%%%%

part_var_info(1).name = 'vel_relative';
part_var_info(1).length = 4;

part_var_info(2).name = 'tau_p';
part_var_info(2).length = 1;

part_var_info(3).name = 'k_sgs';
part_var_info(3).length = 1;

part_var_info(4).name = 'tau_r';
part_var_info(4).length = 6;

part_var_info(5).name = 'sgs_dissipation';
part_var_info(5).length = 1;

part_var_info(6).name = 'sgs_diss';
part_var_info(6).length = 1;

part_var_info(7).name = 't_sgs';
part_var_info(7).length = 1;

part_var_info(8).name = 't_l';
part_var_info(8).length = 3;

part_var_info(9).name = 'b_i';
part_var_info(9).length = 3;

part_var_info(10).name = 'c_i';
part_var_info(10).length = 3;

part_var_info(11).name = 'c_star';
part_var_info(11).length = 3;

part_var_info(12).name = 'b_ii';
part_var_info(12).length = 3;



save('part_var_info.mat',"part_var_info")
