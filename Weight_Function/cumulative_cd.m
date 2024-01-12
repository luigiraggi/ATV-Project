clear all
clc
close all

%% ESTRAZIONE DATI
cumulative_coeffs_z = table2array(importfile_bin_bike('binfield_bike_z.txt'));
z_coord = table2array(importfile_zcoord('z_coordinates.txt'));
cumulative_coeffs_y = table2array(importfile_bin_bike_y('binfield_bike_y.txt'));
y_coord = table2array(importfile_zcoord('y_coordinates.txt'));
cumulative_coeffs_x = table2array(importfile_bin_bike_x('binfield_bike_x.txt'));
x_coord = table2array(importfile_zcoord('x_coordinates.txt'));

delta_z = z_coord(20);
delta_y = abs(y_coord(1));

%% MEAN Cd_z 
for i = 1:size(cumulative_coeffs_z,1)
    counter = 1;
    for j = 1:20
        Cd_z(j,i) = cumulative_coeffs_z(i,counter);
        counter = counter + 12;
    end
end
Cd_z = mean(Cd_z, 2);

%% MEAN Cd_y 
for i = 1:size(cumulative_coeffs_y,1)
    counter = 1;
    for j = 1:20
        Cd_y(j,i) = cumulative_coeffs_y(i,counter);
        counter = counter + 12;
    end
end
Cd_y = mean(Cd_y(:,101:501), 2);

%% MEAN Cd_x
for i = 1:size(cumulative_coeffs_x,1)
    counter = 1;
    for j = 1:20
        Cd_x(j,i) = cumulative_coeffs_x(i,counter);
        counter = counter + 12;
    end
end
Cd_x = mean(Cd_x, 2);

%% WEIGHTING COEFF z
mean_Cd_slope_z = Cd_z(end,1)/1;

syms z
weight_z = 0;
for i = 1:20
    z2 = z_coord(i);
    Cd_z2 = Cd_z(i);
    if i == 1
        z1 = 0;
        Cd_z1 = 0;
    else
        z1 = z_coord(i-1);  
        Cd_z1 = Cd_z(i-1);
    end
    spacing = (z2-z1)/delta_z;
    Cd_slope_z = (Cd_z2-Cd_z1)/spacing;
    weight_z = weight_z + piecewise((z>z1)&(z<=z2), Cd_slope_z - mean_Cd_slope_z,0);
end

%save("../Buoyancy/weight_z.mat","weight_z")

%% WEIGHTING COEFF y
mean_Cd_slope_y = Cd_y(end,1)/1;

syms y
weight_y = 0;
for i = 1:19
    y1 = y_coord(i);  
    Cd_y1 = Cd_y(i);
    if i == 19
        y2 = 0;
        Cd_y2 = 0.6392;
    else
        y2 = y_coord(i+1);
        Cd_y2 = Cd_y(i+1);
    end
    spacing = (y2-y1)/delta_y;
    Cd_slope_y = (Cd_y2-Cd_y1)/spacing;
    weight_y = weight_y + piecewise((y>y1)&(y<=y2), Cd_slope_y - mean_Cd_slope_y,0);
end

%save("../Buoyancy/weight_y.mat","weight_y")

%% PLOT

figure
hold on
grid on
plot(Cd_z, z_coord, 'ko')
plot(Cd_z, z_coord, 'k')
title(['$$cumulative\;C_D$$'],'interpreter','latex','FontWeight','bold','Color','k');
xlabel('$$C_D$$','interpreter','latex');
ylabel('z','interpreter','latex');

figure
hold on
grid on
plot(y_coord,Cd_y, 'ko')
plot(y_coord,Cd_y, 'b')
title(['$$cumulative\;C_D$$'],'interpreter','latex','FontWeight','bold','Color','k');
xlabel('$$C_D$$','interpreter','latex');
ylabel('y','interpreter','latex');

figure
hold on
grid on
plot(x_coord,Cd_x, 'ko')
plot(x_coord,Cd_x, 'r')
title(['$$cumulative\;C_D$$'],'interpreter','latex','FontWeight','bold','Color','k');
xlabel('x','interpreter','latex');
ylabel('$$C_D$$','interpreter','latex');

 figure
 fplot(weight_z)

 figure
 fplot(weight_y)
