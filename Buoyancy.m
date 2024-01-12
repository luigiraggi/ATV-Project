clear all
clc
close all

%% LOADING OF WEIGHT FUNCTION
syms z
load("weight_z.mat")

%% DEFINITIONS
y_iter = floor(0.225/0.01);  %number in y spacing in plotoverline
z_iter = 1.6/0.05; %number in z spacing in plotoverline
x_spacing = 15/1001; %x spacing in plotoverline, [m]
corrective_coeff = 3; %coefficient needed to scale the weight function, set to 0 for avoid weight influence

cd_bike = 0.6392; %solo bike cd(CFD)
cd_wwb = 0.6043;  %cd coefficient of wagon with bikes at 3m
cd_w = 0.6150;  %cd coefficient of clean wagon at 3m
Cd = [0.5337, 0.5900, 0.6250, 0.6258, 0.6354, 0.6380, 0.5282, 0.5871, 0.6245, 0.6253, 0.6340, 0.6379]; %cd(CFD) values: 1-6 w, 7-12 wwb
distances = [1, 2, 4, 5, 6, 8]; %distances between bicycle and flagship, measured in meters
delta_Cd = cd_bike*ones(1,12) - Cd;


%% WEIGHTING ALGORITHM
for k = 1:2 
    counter = 1;
    if k == 1
       nome_vehicle = 'plot_p_w_';
    else 
       nome_vehicle = 'plot_p_wwb_';
    end
    y_coord = double(0.0);
    for i = 1:y_iter
        y_char = char(string(y_coord));
        if y_coord == 0
           y_char = [y_char '.' y_char];
        end
        z_coord = double(0.0);
        for j = 1:z_iter
            z = z_coord;
            w = double(subs(weight_z));
            if corrective_coeff*w < -1   %check that the correction dosen't leads to negative pressure contribution
               coeff_z = -1;
            else
            coeff_z = corrective_coeff*w;
            end
            if z_coord == 0 
               z_char = [char(string(z_coord)) '.' char(string(0))];
            elseif floor(z_coord) == 1 && mod(z_coord,1) < 0.05
               z_char = [char(string(z_coord)) '.' char(string(0))];
            else
               z_char = char(string(z_coord)); 
            end
            nome_file = [nome_vehicle y_char '_' z_char '.csv'];
            p(:, counter,k) = (1+coeff_z)*table2array(import_p(nome_file));  
            z_coord = z_coord + 0.05;
            counter = counter +1;
        end
        y_coord = y_coord - 0.01;
    end
end

%% P WEIGHTED TREND
p_w = mean(p(:,:,1), 2);   %clean wagon pressure distribution
p_wwb = mean(p(:,:,2), 2);  %wagon with bikes pressure distribution

%% LAMBDA COEFFICIENT EVALUATION
counter = 1;
for i = 1:2 
    if i == 1
       p_mean = p_w;
    else 
       p_mean = p_wwb;
    end
    figure 
    hold on
    grid on
    plot(linspace(0.5, 15.5, 1001), p_mean)
    for j = 1:6
        if j == 6   %x(1) coordinate definition, remove 0.25 for getting the real down stream coordinate
         w(j,i) = floor((distances(j)-0.5 + 0.25)/x_spacing);
        else 
         w(j,i) = floor((distances(j)-0.5 + 0.25)/x_spacing) +1;
        end
        z(j,i) = w(j,i) + floor(0.93/x_spacing);  %x(2) coordinate definition, from x(1) is 0.93m, change 0.93 with 1.4 for getting the real up stream coordinate
        dp(counter) = (p_mean(z(j,i))-p_mean(w(j,i)))/0.93;  %linear approximtion of the derivative
        counter = counter +1;
        plot(distances(j)+0.4, p_mean(w(j,i)),'x', distances(j) + 0.4+0.93, p_mean(z(j,i)),'x')
    end  
end

lambda = -dp'\delta_Cd';  %linear coefficient evaluation, through least squares method

%% WWB 3m ESTIMATION
x1 = floor((3-0.5+0.4)/x_spacing);
x2 = x1 + floor(0.93/x_spacing);
dp_wwb_3m = (p_wwb(x2)-p_wwb(x1))/0.93;

err_percentage_wwb = (cd_wwb - (cd_bike + lambda*dp_wwb_3m))/cd_wwb;
cd_wwb_3_estimated = cd_bike + lambda*dp_wwb_3m;

%% W 3m ESTIMATION
x1 = floor((3-0.5+0.4)/x_spacing);
x2 = x1 + floor(0.93/x_spacing);
dp_w_3m = (p_w(x2)-p_w(x1))/0.93;

err_percentage_w = (cd_w - (cd_bike + lambda*dp_w_3m))/cd_w;
cd_w_3_estimated = cd_bike + lambda*dp_w_3m;






