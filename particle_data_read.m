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
