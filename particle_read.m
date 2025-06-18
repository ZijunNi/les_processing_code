function data_particle = particle_read(begin,step,ending,path,save_name,post_rank)
% 本函数主要读取'particle_solver'程序中得到的经后处理后的粒子计算结果
%-----输入参数：
% begin - 开始步数
% step - 'particle_solver'中导出数据的间隔
% ending - 结束步数
% path - 数据文件的绝对路径
% save_name - 读取数据后保存.mat文件的文件名，默认保存在./data文件夹下
% post_rank - 进行后处理时的并行块数（与计算的并行块不一定相同）

i = 1;

save_name = fullfile("data",save_name);

if exist(save_name,'file')
    txt = input('待保存文件已存在，请检查。若要覆盖之，请输入Y，若要增加之，请输入A \n','s');
    if ~strcmp(txt, 'Y')
        return
    end
end


for num=begin:step:ending
    tic
    % 转换为固定长度的字符串，长度为7，前面填充零
    formattedStr = sprintf('%07d', num);
    
    %filename_particle =[path,'Particle',formattedStr,'.0000.dat'];
    if(post_rank==0)
        filename_particle = [path,'Particle',formattedStr,'.0000.dat'];
        particle_data = readmatrix(filename_particle);
        % write(60+rank,'(100E18.6)')xyzp(:,i),uvwpnew(:,i),Auvwpnew(:,i),ee(0:3,i),lna(:,i),uvwfnew(:,i),real(idtype),real(idrank),real(idnum)
        data_particle{i} = particle_data;
    elseif(post_rank~=0)
        particle_data = [];
        %scatter-complete-2000000-0000.dat
        for j = 1:post_rank
            particle_data_mid = [];
            PostRank = sprintf('%04d', j-1);
            filename_particle = [path,'scatter-complete-',formattedStr,'-',PostRank,'.dat'];
            wrk_mid = readmatrix(filename_particle);
            particle_data_mid(:,5:7) = wrk_mid(:,10:12);%第10到12列为压力梯度脉动项
            particle_data_mid(:,4) = wrk_mid(:,1);%第四列为编号
            particle_data_mid(:,1:3) = wrk_mid(:,2:4);%第一至三列为位置坐标，以上安排为了兼容旧格式的数据
            particle_data = [particle_data;particle_data_mid];
        end
        data_particle{i} = particle_data;
    end
    i = i+1;
    disp(['已完成',num2str((num-begin)/(ending-begin)*100),'%'])
    toc
end
% delete(p);
save(save_name,'data_particle','begin','step','ending');

end