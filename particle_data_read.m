% 本代码主要读取'particle_solver'程序中得到的经后处理后的粒子计算结果，并进行处理


%% 读取部分 
clear,clc;
dns_ref = particle_read(00000,2000,1000000,'F:\本机文档\data\data_DNS_st_1\post\postout\', ...
        'dns_st_1_ref.mat',16);
fdns_ref = particle_read(00000,2000,1300000,'F:\本机文档\data\data_20250606_st_1\post\postout\', ...
        '20250606_st_1.mat',16);

%% 读取尾部部分
clear,clc;
dns_ref = particle_read(900000,1000,1000000,'F:\本机文档\data\data_DNS_st_1\post\postout\', ...
        'dns_st_1_ref_tail.mat',16);
fdns_ref = particle_read(900000,1000,1000000,'F:\本机文档\data\data_20250606_st_1\post\postout\', ...
        '20250606_st_1_tail.mat',16);

%% 收敛性检查部分

N=40;
Re_tau = 180;

figure(3)

tic
parfor i =1:500
    [dns_pdf_particle,~,~] = particle_pdf_get(N,i,i+1,"dns_st_1_ref.mat",Re_tau);
    dns_max_vel(i) = max(dns_pdf_particle);
    steps(i) = (i-1)*2000;
end

parfor i = 1:650
    [fdns_pdf_particle,~,~] = particle_pdf_get(N,i,i+1,"20250606_st_1.mat",Re_tau);
    fdns_max_vel(i) = max(fdns_pdf_particle);
    steps_fdns(i) = (i-1)*2000;
end
toc

plot(steps,dns_max_vel,'-',DisplayName='DNS')
hold on
plot(steps_fdns,fdns_max_vel,'-',DisplayName='fDNS, CF = (2,2)')
hold off
legend();
xlabel('Steps')
ylabel('Maxinum of Particle PDF')


%% PDF绘制部分
figure;
N = 50;
begin_step = 490;
end_step = 500;
Re_tau = 180;

[dns_pdf_particle,~,dns_number] = particle_pdf_get(N,begin_step,end_step,'dns_st_1_ref.mat',Re_tau);
[fdns_pdf_particle,exp_y_plus,fdns_number] = particle_pdf_get(N,begin_step,end_step,'20250606_st_1.mat',Re_tau);

semilogx(exp_y_plus,fdns_pdf_particle,'-s','LineWidth',2,'DisplayName','New Test, fDNS, CF = (2,2)')
hold on
semilogx(exp_y_plus,dns_pdf_particle,'-o','LineWidth',2,'DisplayName','DNS Ref.')
hold off

title([ 'Total Number of Particles = ',num2str(dns_number)]);
xlabel('y^+')
ylabel('PDF')

legend()

%% 粒子流向速度剖面检查
data_fold = 'processed_data';
%post输出：int(lna(1,i)),xyzp(1:3,i),lna(2:3,i),uvwpnew(:,i),ee(1:3,i)
data_file = 'dns_st_1_ref.mat';
data_loc = fullfile(data_fold,data_file);
data = load(data_loc);
dns_data = data.data_particle; % DNS 参考数据

data_file = '20250606_st_1.mat';
data_loc = fullfile(data_fold,data_file);
data = load(data_loc);
fdns_data = data.data_particle; % fDNS 数据

clearvars -except fdns_data dns_data


grid = readmatrix('./data/grid_test.dat');
grid = grid(:,2)/180;
grid = [0;grid];


single_data_fdns = [];
single_data_dns = [];

for j = 100:130 % 数据所在时间步
    single_data_fdns = [single_data_fdns;fdns_data{j}];
    single_data_dns = [single_data_dns;dns_data{j}];
end



for i = 1:length(grid)-1
    %%%% fDNS
    loc1 = find(single_data_fdns(:,3) > grid(i));
    loc2 = find(single_data_fdns(loc1,3) < grid(i+1));
    loc = loc1(loc2);
    
    fdns_mean_u_p(i) = mean(single_data_fdns(loc, 7));
    fdns_mean_u_p_wall(i) = mean(single_data_fdns(loc, 9));
    fdns_rms_u_p(i) = sqrt(mean((single_data_fdns(loc, 7) - fdns_mean_u_p(i)).^2));
    y_pos(i) = 180 * (grid(i) + grid(i+1)) / 2;

    %%%% DNS
    loc1 = find(single_data_dns(:,3) > grid(i));
    loc2 = find(single_data_dns(loc1,3) < grid(i+1));
    loc = loc1(loc2);
    
    dns_mean_u_p(i) = mean(single_data_dns(loc, 7));
    dns_mean_u_p_wall(i) = mean(single_data_fdns(loc, 9));
    dns_rms_u_p(i) = sqrt(mean((single_data_dns(loc, 7) - dns_mean_u_p(i)).^2));
    y_pos(i) = 180 * (grid(i) + grid(i+1)) / 2;
end
%%
figure;
subplot(2,1,1)
semilogx(y_pos, fdns_mean_u_p,DisplayName='fDNS');
hold on
semilogx(y_pos, dns_mean_u_p,DisplayName='DNS');
hold off
legend('Location','northwest',Interpreter='latex');
xlabel('$y^+$',Interpreter='latex')
ylabel('$\overline{u_{p,x}}$',Interpreter='latex')
xlim([0.1 180])
ylim([0 1.4])
% daspect([1/(1.4) 1/(180-0.1) 1])

subplot(2,1,2)
plot(y_pos, fdns_rms_u_p,DisplayName='fDNS');
hold on
plot(y_pos, dns_rms_u_p,DisplayName='DNS');
hold off
legend('Location','northwest',Interpreter='latex');
xlabel('$y^+$',Interpreter='latex')
ylabel('$\sqrt{\overline{(u^\prime_{p,x})^2}}$',Interpreter='latex')
xlim([0.1 180])
% ylim([0 1.4])
% daspect([1/(1.4) 1/(180-0.1) 1])
%%
figure;
plot(y_pos, fdns_mean_u_p_wall,DisplayName='fDNS');
hold on
plot(y_pos, dns_mean_u_p_wall,DisplayName='DNS');
hold off
legend('Location','northwest',Interpreter='latex');
xlabel('$y^+$',Interpreter='latex')
ylabel('$\overline{u_{p,z}}$',Interpreter='latex')
% xlim([0.1 180])
% ylim([0 1.4])