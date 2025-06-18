function [pdf_particle,exp_y_plus,total_num_part] = particle_pdf_get(N,begin_step,end_step,filename,Re_tau)
% 本函数处理particle_read.m脚本读取并保存的粒子数据文件，并根据数据文件绘制粒子PDF在垂向上的分布
% 输入：
% N - 统计区间数目，函数中可切换为对数等距区间或等距区间
% begin_step & end_step - PDF统计的起止范围，用数据集的个数表示（1开始间隔为1）
% filename - 读取的粒子数据文件名，默认在./data文件夹下
% Re_tau -  读取的粒子所在流场的摩擦雷诺数

% 输出
% pdf_particle - 粒子数密度PDF在垂向的分布
% exp_y_plus - pdf_particle对应的垂向坐标
% total_num_part - 参与统计的粒子总数

filename = fullfile("data",filename);

load(filename)

delta_visc = 1/Re_tau;

counter = zeros(1,N-1);


total_number = length(data_particle{1}(:,1));
exp_y = linspace(-2.5,0,N);exp_y = 10.^(exp_y);%对数等距区间
% exp_y = linspace(1/180,1,N);%exp_y(1)=0;%等距区间


for i = begin_step:end_step

    num_particle(i) = length(data_particle{i});
    y_particle = data_particle{i}(:,2);

    for j = 1:length(y_particle)% 将上半槽道映射到下半
        if y_particle(j)>0%>0:翻转槽道；>1:上半下半叠加
            y_particle(j) = 2-y_particle(j);
        end
    end

    
    for j = 1:N-1
        number = length(find(y_particle>=exp_y(j)&y_particle<exp_y(j+1)));
        if i == begin_step
            delta_block(j) = (exp_y(j+1)-exp_y(j));
            exp_y_plus(j) = (exp_y(j)+exp_y(j+1))/(2*delta_visc);
        end
        counter(j) = counter(j)+number/delta_block(j);
    end
end

total_num_part = sum(num_particle);
for i=1:N-1
    pdf_particle(i) = counter(i)/(delta_block(i));
end
pdf_particle=counter/length(data_particle);
pdf_in = sum(num_particle)/length(data_particle);
pdf_particle=pdf_particle/pdf_in;

end
