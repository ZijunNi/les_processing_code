% 本代码尝试从JHTDB中读取Lagrange数据
clear,clc

time = 10.0;
points = [3,0.1,3];

velocity = getVelocity3D(time, points);



%%
% Lagrange粒子追踪 - MATLAB脚本版本
clear; clc;

% 设置追踪参数
initial_time = 0.0;
initial_position = [1.0, 0.3, 1.0]; % [x, y, z]
dt = 0.0013;
num_steps = 20000;

% 周期边界条件参数
x_period = 8 * pi;
z_period = 3 * pi;

% 预分配结果数组以提高效率
time_results = zeros(num_steps + 1, 1);
position_results = zeros(num_steps + 1, 3);
velocity_results = zeros(num_steps + 1, 3);

% 初始状态
current_time = initial_time;
current_position = initial_position;

% 获取初始速度
current_velocity = getVelocity3D(current_time, current_position);

% 保存初始状态
time_results(1) = current_time;
position_results(1, :) = current_position;
velocity_results(1, :) = current_velocity;

% 开始追踪循环
for step = 1:num_steps
    % 获取当前速度
    current_velocity = getVelocity3D(current_time, current_position);
    
    % 使用前向欧拉法推进粒子位置
    new_position = current_position + current_velocity * dt;
    new_time = current_time + dt;
    
    % 处理周期边界条件
    % x方向: 0 ~ 8pi
    if new_position(1) < 0
        new_position(1) = new_position(1) + x_period;
    elseif new_position(1) >= x_period
        new_position(1) = new_position(1) - x_period;
    end
    
    % z方向: 0 ~ 3pi
    if new_position(3) < 0
        new_position(3) = new_position(3) + z_period;
    elseif new_position(3) >= z_period
        new_position(3) = new_position(3) - z_period;
    end
    
    % 获取新位置的速度
    new_velocity = getVelocity3D(new_time, new_position);
    
    % 保存结果
    time_results(step + 1) = new_time;
    position_results(step + 1, :) = new_position;
    velocity_results(step + 1, :) = new_velocity;
    
    % 更新当前状态
    current_time = new_time;
    current_position = new_position;
end

% 显示前几步结果
% fprintf('Lagrange粒子追踪结果（前10步）:\n');
% fprintf('步数 |   时间   |     位置(x,y,z)     |     速度(vx,vy,vz)\n');
% fprintf('-----|----------|---------------------|----------------------\n');
% for i = 1:min(10, num_steps + 1)
%     pos = position_results(i, :);
%     vel = velocity_results(i, :);
%     fprintf('%3d  | %8.3f | (%6.3f,%6.3f,%6.3f) | (%6.3f,%6.3f,%6.3f)\n', ...
%             i-1, time_results(i), pos(1), pos(2), pos(3), vel(1), vel(2), vel(3));
% end

%%
% 可视化结果
figure;

% 3D轨迹图
subplot(2, 2, 1);
plot3(position_results(:,1), position_results(:,2), position_results(:,3), 'b-', 'LineWidth', 2);
hold on;
plot3(position_results(1,1), position_results(1,2), position_results(1,3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot3(position_results(end,1), position_results(end,2), position_results(end,3), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('粒子3D运动轨迹');
legend('轨迹', '起点', '终点', 'Location', 'best');
grid on;

% 各方向位置随时间变化
subplot(2, 2, 2);
plot(time_results, position_results(:,1), 'r-', 'LineWidth', 2); hold on;
plot(time_results, position_results(:,2), 'g-', 'LineWidth', 2);
plot(time_results, position_results(:,3), 'b-', 'LineWidth', 2);
xlabel('时间'); ylabel('位置');
title('位置随时间变化');
legend('X', 'Y', 'Z', 'Location', 'best');
grid on;

% 各方向速度随时间变化
subplot(2, 2, 3);
plot(time_results, velocity_results(:,1), 'r-', 'LineWidth', 2); hold on;
plot(time_results, velocity_results(:,2), 'g-', 'LineWidth', 2);
plot(time_results, velocity_results(:,3), 'b-', 'LineWidth', 2);
xlabel('时间'); ylabel('速度');
title('速度随时间变化');
legend('Vx', 'Vy', 'Vz', 'Location', 'best');
grid on;

% X-Z平面投影（考虑周期边界）
subplot(2, 2, 4);
% 处理周期边界以显示连续轨迹
x_mod = mod(position_results(:,1), x_period);
z_mod = mod(position_results(:,3), z_period);
plot(x_mod, z_mod, 'b-', 'LineWidth', 2);
hold on;
plot(x_mod(1), z_mod(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(x_mod(end), z_mod(end), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
xlabel('X (mod 8π)'); ylabel('Z (mod 3π)');
title('X-Z平面投影（周期边界）');
legend('轨迹', '起点', '终点', 'Location', 'best');
grid on;
xlim([0, x_period]);
ylim([0, z_period]);

save('JHTDB_particle_tracking_results.mat', 'time_results', 'position_results', 'velocity_results', ...
     'initial_time', 'initial_position', 'dt', 'num_steps', 'x_period', 'z_period');

fprintf('\n结果已保存到 particle_tracking_results.mat\n');
fprintf('总共追踪了 %d 步，从 t=%.2f 到 t=%.2f\n', num_steps, initial_time, time_results(end));

%%
clear,clc
load('JHTDB_particle_tracking_results.mat');

all_time = time_results(1:end-5);

all_data = velocity_results(1:end-5,:);
all_loc = position_results(1:end-5,:);

num_steps = length(all_time);

save('JHTDB_lagrange_data.mat', 'all_time', 'all_data',"all_loc");



%%


function velocity_components = getVelocity3D(time, points, authkey, dataset)
%GETVELOCITY3D 获取三维速度分量
%   根据输入的时间和三维坐标，输出对应的三维速度分量
%
%   输入参数：
%   time - 时间点
%   points - 三维坐标点数组，n_points x 3矩阵，每行包含(x,y,z)坐标
%   authkey - 认证密钥（可选，默认为'cn.edu.pku.stu.nizijun-955f3869'）
%   dataset - 数据集名称（可选，默认为'channel'）
%
%   输出参数：
%   velocity_components - 三维速度分量，n_points x 3矩阵，每行包含(u,v,w)速度分量

    % 设置默认参数
    if nargin < 3
        authkey = 'cn.edu.pku.stu.nizijun-955f3869';
    end
    if nargin < 4
        dataset = 'channel';
    end
    
    % 设置固定参数
    variable = 'velocity';
    temporal_method = 'pchip'; 
    spatial_method = 'lag4';
    spatial_operator = 'field';

    % 获取数据点数
    n_points = size(points, 1);
    
    % ---- GetData ----
    tic
    fprintf('\nRequesting %s at %i points for time = %.1f...\n', variable, n_points, time);
    result = getData_JHTDB(authkey, dataset, variable, time, temporal_method, spatial_method, spatial_operator, points);
    toc
    
    % 提取速度分量（假设result的第三列包含速度值）
    % 注意：根据getData函数的实际返回结构，可能需要调整列索引
    velocity_components = result; % 如果返回的是标量速度值
    
    % 如果getData返回的是三个速度分量，可能需要调整为：
    % velocity_components = result(:, 3:5); % 假设3-5列分别是u,v,w分量
    
    fprintf('Data retrieval completed. Retrieved %i velocity values.\n', n_points);
end

% 可选：添加一个辅助函数用于生成规则网格点
function points = generateGridPoints(nx, ny, nz, x_range, y_range, z_range)
%GENERATEGRIDPOINTS 生成规则的三维网格点
%
%   输入参数：
%   nx, ny, nz - 各方向的网格点数
%   x_range - x方向范围 [x_min, x_max]
%   y_range - y方向范围 [y_min, y_max]  
%   z_range - z方向范围 [z_min, z_range]
%
%   输出参数：
%   points - 生成的三维坐标点数组

    n_points = nx * ny * nz;
    points = zeros(n_points, 3);
    
    x_points = linspace(x_range(1), x_range(2), nx);
    y_points = linspace(y_range(1), y_range(2), ny);
    z_points = linspace(z_range(1), z_range(2), nz);
    
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                points((i - 1) * ny * nz + (j - 1) * nz + k, 1) = x_points(i);  
                points((i - 1) * ny * nz + (j - 1) * nz + k, 2) = y_points(j);
                points((i - 1) * ny * nz + (j - 1) * nz + k, 3) = z_points(k);
            end
        end
    end
end

% 使用示例：
% % 示例1：使用默认参数
% time = 10.0;
% nx = 100; ny = 5; nz = 30;
% points = generateGridPoints(nx, ny, nz, [0, 8*pi], [-1, -1+5/1000], [0.2*pi, 0.4*pi]);
% velocity = getVelocity3D(time, points);
%
% % 示例2：自定义参数
% time = 15.0;
% custom_points = [1.0, -0.5, 1.0; 2.0, -0.3, 1.2; 3.0, -0.1, 1.5]; % 3个自定义点
% custom_authkey = 'your-auth-key-here';
% custom_dataset = 'your-dataset-name';
% velocity = getVelocity3D(time, custom_points, custom_authkey, custom_dataset);