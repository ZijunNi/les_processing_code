clc,clear

load('./data/20250606_st_1.mat');
number =  (ending - begin)/step +1;
jump = 1;%每隔jump个数据样本绘制一帧

figure; 
xlim([0,4*pi])
ylim([0,2*pi])

    for i = 1:jump:number
    % 绘制当前帧的粒子
        number_particle = data_particle{i}(:,1);
        x_particle = data_particle{i}(:,1);
        y_particle = data_particle{i}(:,2);
        z_particle = data_particle{i}(:,3);
        %x_particle_draw = x_particle(y_particle<y_range(2)/360&y_particle>y_range(1)/360);
        %y_particle_draw = y_particle(y_particle<y_range(2)/360&y_particle>y_range(1)/360);
        %z_particle_draw = z_particle(y_1particle<y_range(2)/360&y_particle>y_range(1)/360);
        y_range = [00,30];
        re_tau = 180;%re_tau=1/2时上述y_range需为物理范围，而非y^+范围
        particle_draw1 = x_particle;%(z_particle<(1+0.01)*pi & z_particle>(1-0.01)*pi);%(y_particle<y_range(2)/(2*re_tau)&y_particle>y_range(1)/(2*re_tau));%(z_particle<(1+0.01)*pi & z_particle>(1-0.01)*pi);%
        particle_draw2 = y_particle;%(z_particle<(1+0.01)*pi & z_particle>(1-0.01)*pi);%(y_particle<y_range(2)/(2*re_tau)&y_particle>y_range(1)/(2*re_tau));%(z_particle<(1+0.01)*pi & z_particle>(1-0.01)*pi);%
        clf
        plot(particle_draw1,particle_draw2,'ro','MarkerFaceColor',[102 102 102]/225,'MarkerSize',1);
        %axis([0 pi/20 0 0.2]);    
        num(i) = length(particle_draw1);
        grid on
        daspect([1,1,1])
        title([num2str(i-1),'共有粒子',num2str(num(i)),'个'])
        pause(0.25)
    end



