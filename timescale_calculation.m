%%%%%%%%%%%% 数据读取部分 %%%%%%%%%%%%%%
    % ----------- JHTDB单点时序数据
    clear,clc;
    filename = 'JHTDB_euler.mat';
    load(filename);
    if(contains(filename,'euler'))
        % loc = all_loc;
        loc = [10.33,0.9,4.6];% 手动指定
    elseif(contains(filename,'lagrange'))
        loc(1) = mean(all_loc(:,1));
        loc(2) = mean(all_loc(:,2));
        loc(3) = mean(all_loc(:,3));
    end
    % ----------- JHTDB单点时序数据
 %%
    clear,clc
    % ----------- 自有程序数据
    filename = 'fluid_euler_2025-08-12_17-41-21.mat';
    foldername = 'processed_data';
    load(fullfile(foldername,filename));
    num_sample = 1;
    loc = particleData{num_sample}.data(1,2:4);
    all_data = particleData{num_sample}.data(:,5:7);
    all_time = particleData{num_sample}.data(:,1);
    % ----------- 自有程序数据
%%%%%%%%%%%% 数据读取部分 %%%%%%%%%%%%%%

%%

% ------- 自有程序数据处理：提取特定位置的采样结果
j = 1;
for i = 1:length(particleData)
    location = particleData{i}.data(1,2:4);
    if(location(2)==0.000547238044500000)
        test_data{j} = particleData{i};
        j = j + 1;
    end
end
% ------- 自有程序数据处理：提取特定位置的采样结果

%% JHTDB 单点数据处理

    tic
    x = all_data(:,1);% 流向分量
    mean(x)
    x = x - mean(x);
    t = all_time;
    [rho,lags] = my_xcorr(x,t,0000);
    [x_inte,inte,zero_cross_timescale,zero_cross_point] = integrate_from_start(lags,rho);
    num_sample = 1;
    
    [~,t_l,r_squared,x_fit,y_fit] = lesutils.fit_exp_model(lags, rho, 0); % 对自相关函数进行拟合
    t_l = -1/t_l;% 拟合系数为负倒数

    % 保存结果
    all_loc(num_sample,1:3) = loc;
    all_t_l(num_sample) = t_l;
    all_r_squared(num_sample) = r_squared;
    % 保存结果

    if(1)
        figure;
        sgtitle(['JHTDB Timeseries Data at (',num2str(loc(1)),',',num2str(loc(2)),',',num2str(loc(3)),')'])
        % 第一个子图 
        subplot(1, 2, 1);
        plot(lags, rho,DisplayName='Autocorrelation Function');
        hold on
        plot(x_fit,y_fit,'r-.',DisplayName=['Exponential Fit $\rho(s)=\exp(s/T_L)$' newline ' with $T_L = $',num2str(t_l)]);
        grid on;
        legend('Interpreter','latex')
        title('Autocorrelation Funtion')
        
        % 第二个子图 
        subplot(1, 2, 2);
        plot(x_inte,inte);
        hold on
        line([0.5 x_inte(end)],[t_l,t_l],'linestyle','--');
        line([0.5 x_inte(end)],[zero_cross_timescale,zero_cross_timescale],'linestyle','--');
        text(1,t_l+0.005,['$T_L = $', num2str(t_l)],'Interpreter','latex')
        text(1,zero_cross_timescale+0.005,['$T_L = $', num2str(zero_cross_timescale)],'Interpreter','latex')
        grid on;
        title('Integration of Autocorrelation Funtion')
    end
    toc
% end

save(['timescale_',filename],"all_r_squared","all_t_l","all_loc")


%% 自有数据处理：提取特定位置范围内的粒子数据

tic
num_sample = 669;
j = 1;
% for num_sample = 1:length(particleData)

    if(contains(filename,'euler'))
        loc = particleData{num_sample}.data(1,2:4);
    elseif(contains(filename,'lagrange'))
        loc(1) = mean(particleData{num_sample}.data(:,2));
        loc(2) = mean(particleData{num_sample}.data(:,3));
        loc(3) = mean(particleData{num_sample}.data(:,4));
    end
    if(abs(loc(2)-0.1)<0.01)
        filtered_particleData{j} = particleData{num_sample};
        j = j+1;
    end
% end


%%
plot(0:length(matrix_sequence)-1, autocorr);
xlabel('延迟 \tau');
ylabel('自相关 R(\tau)');
grid on;
%% 自有数据时间尺度计算：单个粒子的计算

    %%%%%------原始代码，对单个粒子求解
        all_data = particleData{num_sample}.data(:,5:7);
        all_time = particleData{num_sample}.data(:,1);
        
        
        x = all_data(:,1);% 流向分量
        x = x - mean(x);
        t = all_time;
        [rho,lags] = my_xcorr(x,t,00000);
        [x_inte,inte] = integrate_from_start(lags,rho);
        
        [~,t_l,r_squared,x_fit,y_fit] = lesutils.fit_exp_model(lags, rho, 0); % 对自相关函数进行拟合
        t_l = -1/t_l;% 拟合系数为负倒数
    
        % 保存结果
        all_loc(num_sample,1:3) = loc;
        all_t_l(num_sample) = t_l;
        all_r_squared(num_sample) = r_squared;
        % 保存结果
    %%%%%------原始代码，对单个粒子求解

    if(1)
        figure;
        sgtitle(['DNS Timeseries Data at (',num2str(loc(1)),',',num2str(loc(2)),',',num2str(loc(3)),')'])
        % 第一个子图 
        subplot(1, 2, 1);
        plot(lags, rho,DisplayName='Autocorrelation Function');
        hold on
        plot(x_fit,y_fit,'r-.',DisplayName=['Exponential Fit $\rho(s)=\exp(s/T_L)$' newline ' with $T_L = $',num2str(t_l)]);
        grid on;
        legend('Interpreter','latex')
        title('Autocorrelation Funtion')
        
        % 第二个子图 
        subplot(1, 2, 2);
        plot(x_inte,inte);
        hold on
        line([0.5 x_inte(end)],[t_l,t_l],'linestyle','--');
        text(1,t_l+0.005,['$T_L = $', num2str(t_l)],'Interpreter','latex')
        grid on;
        title('Integration of Autocorrelation Funtion')
    end
    toc
% end

save(['timescale_',filename],"all_r_squared","all_t_l","all_loc")

%% 绘图部分
y_plus_loc = all_loc(:,2)*180;
scatter(y_plus_loc,all_t_l,all_r_squared/min(all_r_squared),DisplayName='$T_L$ Values');
% set(gca,'XScale','log')
hold on

if(contains(filename,'euler'))
    uniq_loc = unique(y_plus_loc);
    mean_t_l = zeros(size(uniq_loc));
    
    % 遍历每个唯一值
    for i = 1:length(uniq_loc)
        % 找到当前唯一值在A中的所有索引
        indices = (y_plus_loc == uniq_loc(i));
        % 提取对应B中的值
        valuesFromB = all_t_l(indices);
        % 计算平均值
        mean_t_l(i) = mean(valuesFromB);
    end
elseif(contains(filename,'lagrange'))
    temp = load("grid_test.mat");
    uniq_loc = temp.grid_test.data(:,1);
    uniq_loc = [0;uniq_loc(1:2:end);2]*180;
    
    % 确保A是单调递增的
    uniq_loc = sort(uniq_loc);
    
    % 初始化结果数组
    mean_t_l = zeros(length(uniq_loc)-1, 1); % 存储每个区间的平均值
    midpoints = zeros(length(uniq_loc)-1, 1);      % 存储每个区间的中点
    
    % 遍历每个区间
    for i = 1:length(uniq_loc)-1
        % 获取当前区间的上下界
        lower_bound = uniq_loc(i);
        upper_bound = uniq_loc(i+1);
        
        % 计算区间中点
        midpoints(i) = (lower_bound + upper_bound) / 2;
        
        % 找到B值在当前区间内的索引
        indices = find(y_plus_loc >= lower_bound & y_plus_loc < upper_bound);
        
        % 如果找到符合条件的值
        if ~isempty(indices)
            % 提取对应的C值
            c_values = all_t_l(indices);
            
            % 计算平均值
            mean_t_l(i) = mean(c_values);
        else
            % 如果没有找到符合条件的值，设置为NaN
            mean_t_l(i) = NaN;
        end
    end
    uniq_loc = midpoints;

end

plot(uniq_loc,mean_t_l,DisplayName='Mean $T_L$',LineWidth=2)
xlim([0 360])
xlabel('$y^+$',Interpreter='latex')
ylabel('Lagrange Timescale $T_L$',Interpreter='latex')
legend(Interpreter='latex')
