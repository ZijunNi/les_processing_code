% 本代码功能是读取“particle_solver”中les_part_debug_mpi函数的输出结果，并进行合并等处理
format long
clear,clc

target_folder = './data/390debug';% 数据文件夹路径

step = 389;% 待读取的时间步

%% T_l和tau_p的对比验证

T_l_i = part_data_read(target_folder,'t_l_scatter',step);% 读取T_l_i

tau_p = part_data_read(target_folder,'tau_p_scatter',step);% 读取'tau_p'



% 检查获取数据的编号是否一致
check = sum(abs(tau_p(:,1)-T_l_i(:,1)));
if(check)
    error('待对比的数据中的粒子编号不一致');
end
error_points = 0;
% 绘制散点图：
scatter(T_l_i(:,2),tau_p(:,2),'g.',DisplayName='$T_{L,1}$');
error_points = error_points + length(find(T_l_i(:,2)>tau_p(:,2)));
hold on
scatter(T_l_i(:,3),tau_p(:,2),'r.',DisplayName='$T_{L,2}=T_{L,3}$');
% scatter(T_l_i(:,4),tau_p(:,2),'.',DisplayName='$T_{L,3}$');
error_points = error_points + length(find(T_l_i(:,3)>tau_p(:,2)));

plot([1e-2,max(tau_p(:,2))],[1e-2,max(tau_p(:,2))],'b-.',DisplayName='X = Y',LineWidth=2);
hold off
set(gca,'XScale','log')
set(gca,'YScale','log')
legend('Interpreter','latex');
xlabel('$T_{L,i}$','Interpreter','latex');
ylabel('$\tau_p$','Interpreter','latex');


% 寻找tau_p<T_l的粒子
% a = find(T_l_i(:,2)>tau_p(:,2));
% a = [a;find(T_l_i(:,3)>tau_p(:,2))];
% a = [a;find(T_l_i(:,4)>tau_p(:,2))];
% a = unique(a);
% error_pos = tau_p(a,3:5);
% scatter3(error_pos(:,1),error_pos(:,2),error_pos(:,3),'.')% 绘制其在流场中的位置
% xlim([0 4*pi]);
% ylim([0 2]);
% zlim([0 2*pi]);

%% 物理量验证
figure;
for i = 10:10:391
    step = i;
    delta_x = part_data_read(target_folder,'uvwfnew',step);% 读取粒子位移量
    data_length = 3;%数据向量长度

    info = find(delta_x(:,1)==490023);

    delta_x(info,2:2+data_length-1)
    
    % 绘制散点图，观察各分量与位置的关系
    axis = 1;% 轴数xyz=123
    size = abs(delta_x(:,1+axis));
    max_size = max(size);%用于归一化“气泡”大小
    min_size = min(size);%用于归一化“气泡”大小
    scatter(delta_x(:,1+data_length+1),1 - abs(1-delta_x(:,1+data_length+2)),1e+2*size/max_size);%,'filled');%最大“气泡”的半径取为100，将上半槽道的粒子映射到下半
    colorbar;
    % scatter(delta_x(:,5),1 - abs(1-delta_x(:,6)),log(size/min_size*1.01),'filled');%最小“气泡”的半径取为log(1.01)，将上半槽道的粒子映射到下半
    title(['Max Value = ',num2str(max(size)),',Min Value = ',num2str(min_size)])
    hold on
    xlim([0 4*pi]);
    ylim([1e-4 1.01]);
    set(gca,'YScale','log');
    box on
    pause
    clf
end
hold off
close all

%% 函数验证：对两种函数插值的效果进行对比
clc;
step = 5;
target_folder = './data/test_functions';% 数据文件夹路径
second_lagrange = part_data_read(target_folder,'uvwfnew',step);% 读取数据
linear = part_data_read(target_folder,'uvwfnew_linear',step);% 读取数据

error = abs(second_lagrange(:,2:4)-linear(:,2:4))./abs(second_lagrange(:,2:4));

error = [second_lagrange(:,1),error,second_lagrange(:,5:7)];

function result = part_data_read(target_folder,file_mid_name,step)
    % 获取符合文件名格式的所有文件
    % file_mid_name - 读取文件的中间名
    % step - 待读取的时间步
    % target_folder - 数据文件夹路径

    file_list = dir(fullfile(target_folder, ['step_',num2str(step),'_part_',file_mid_name,'_file_*.txt']));
    if isempty(file_list)
        error('未找到符合格式的文件');
    end
    
    
    % 读取所有数据并合并到一个矩阵all_data
    all_data = [];
    
    for i = 1:length(file_list)
        file_path = fullfile(target_folder, file_list(i).name);
        data = load(file_path);
        % all_ids = [all_ids; data(:,1)];
        all_data = [all_data;data];
    end
    result = sortrows(all_data, 1);% 输出的矩阵均按粒子编号排序
end