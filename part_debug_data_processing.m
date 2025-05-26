% 本代码功能是读取“particle_solver”中les_part_debug_mpi函数的输出结果，并进行合并等处理

clear,clc

target_folder = './data/debug';% 数据文件夹路径

step = 1;% 待读取的时间步

T_l_i = part_data_read(target_folder,'T_l_i',step);% 读取T_l_i

tau_p = part_data_read(target_folder,'tau_p',step);% 读取'tau_p'



% 检查获取数据的编号是否一致
check = sum(abs(tau_p(:,1)-T_l_i(:,1)));
if(check)
    error('待对比的数据中的粒子编号不一致');
end

% 绘制散点图：
scatter(T_l_i(:,2),tau_p(:,2),'.',DisplayName='$T_{L,1}$');
hold on
scatter(T_l_i(:,3),tau_p(:,2),'.',DisplayName='$T_{L,2}$');
scatter(T_l_i(:,4),tau_p(:,2),'.',DisplayName='$T_{L,3}$');
plot([1e-2,max(tau_p(:,2))],[1e-2,max(tau_p(:,2))],DisplayName='X = Y');
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

%% 位移验证

% delta_x = part_data_read(target_folder,'delta_x',step);% 读取粒子位移量
% 
% % 绘制散点图，观察流向位移与位置的关系
% scatter(delta_x(:,5),delta_x(:,6),5e4*abs(delta_x(:,2)));
% xlim([0 4*pi]);
% ylim([0 2]);
% box on


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