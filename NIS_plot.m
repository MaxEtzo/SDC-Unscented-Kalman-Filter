f_l = fopen('lidar_NIS','r');
eps_lidar = fscanf(f_l,'%f');
figure(1)
plot(eps_lidar)
hold on
chi2 = ones(1,length(eps_lidar))*5.99;
plot(chi2)

f_r = fopen('radar_NIS','r');
eps_radar = fscanf(f_r,'%f');

figure(3)
plot(eps_radar)
hold on
chi3 = ones(1,length(eps_radar))*7.81;
plot(chi3)