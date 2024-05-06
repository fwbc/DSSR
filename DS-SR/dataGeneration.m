clear; close all; clc;
format long
%--------------------------------------------------------------------------
%{ 
    生成数据集
  Runge-Kutta for data generation
  
%}

% create data directory
mkdir data

% Parameters
%35 3 28
%10 28 8/3
a=35;
b=3;
c=28;
%d=10;
% Define the time span and step size
tspan = [0 100];
h_step = 0.001; % Step size
t = tspan(1):h_step:tspan(2);

% Initial conditions
%{
x0 = -8;
y0 = -7;
z0 = 27;
%}
x0 = 0.5;
y0 = 0.5;
z0 = 0.5;

%% 混沌系统

%chen
f = @(t, x) [a * (x(2) - x(1)); 
            (c - a) * x(1) - x(1) .* x(3) + c * x(2);
            x(1) .* x(2) - b * x(3)]; 


%{
%洛伦兹
f = @(t, y) [a * (y(2) - y(1)); 
            y(1).*(b-y(3))-y(2);  
            y(1) .* y(2) - c * y(3)];
%}

%{
%Rossler
f = @(t, y) [-a*y(2)-y(3); 
            a*y(1)+b*y(2);
            c+y(3)*(y(1)-d)];            
 %}           

% Initialize arrays for the solution
x = zeros(size(t));
y = zeros(size(t));
z = zeros(size(t));

% Set initial values
x(1) = x0;
y(1) = y0;
z(1) = z0;


%% 四阶Runge-Kutta方法进行数值积分
for i = 1:length(t)-1

    k1 = h_step * f(t(i), [x(i), y(i), z(i)]);
    k2 = h_step * f(t(i) + h_step/2, [x(i) + k1(1)/2, y(i) + k1(2)/2, z(i) + k1(3)/2]);
    k3 = h_step * f(t(i) + h_step/2, [x(i) + k2(1)/2, y(i) + k2(2)/2, z(i) + k2(3)/2]);
    k4 = h_step * f(t(i) + h_step, [x(i) + k3(1), y(i) + k3(2), z(i) + k3(3)]);
    
    x(i+1) = x(i) + (1/6) * (k1(1) + 2*k2(1) + 2*k3(1) + k4(1));
    y(i+1) = y(i) + (1/6) * (k1(2) + 2*k2(2) + 2*k3(2) + k4(2));
    z(i+1) = z(i) + (1/6) * (k1(3) + 2*k2(3) + 2*k3(3) + k4(3));

end


%chen
rhs_1_anal = a * (y - x);
rhs_2_anal = (c - a) * x - x .* z + c * y;
rhs_3_anal = x .* y - b * z;


%{
%洛伦兹
rhs_1_anal = a * (y - x);
rhs_2_anal = x.*(b-z)-y;
rhs_3_anal = x .* y - c * z;
%}

%{
%Rossler
rhs_1_anal = -a*y-z;
rhs_2_anal = a*x+b*y;
rhs_3_anal = c+z.*(x-d);
%}
y1_anal = x;
y2_anal = y;
y3_anal = z;

%%保存文件
% save data
directory = 'data';

% Save the variables in the specified file
save(fullfile(directory, 'data_generated.mat'), 'y1_anal', 'y2_anal', 'y3_anal', 'rhs_1_anal', 'rhs_2_anal', 'rhs_3_anal', 't', 'h_step');
data_for_pysr = [t',  y1_anal', y2_anal',y3_anal', rhs_1_anal', rhs_2_anal', rhs_3_anal'];
writematrix(data_for_pysr, fullfile(directory, 'pysr_data_1.csv'));




%% Plotting the solutions

figure(1)
subplot(3,1,1)
set(gca,'Fontsize',12)
hold on
grid on 
plot(t,y1_anal,'LineWidth',2)
box on
title('y1', 'FontWeight', 'Normal')
xlabel('t');
ylabel('x(t)');

subplot(3,1,2)
set(gca,'Fontsize',12)
hold on
grid on 
plot(t,y2_anal,'LineWidth',2)
box on
title('y2', 'FontWeight', 'Normal')
xlabel('t');
ylabel('y(t)');

subplot(3,1,3)
set(gca,'Fontsize',12)
hold on
grid on 
plot(t,y3_anal,'LineWidth',2)
box on
title('y3', 'FontWeight', 'Normal')
xlabel('t');
ylabel('z(t)');

figure(2)
 plot3(y1_anal,y2_anal,y3_anal,'LineWidth',2, 'Color','#15607a')
 title('Chen混沌系统');
 xlabel('x');
 ylabel('y');
 zlabel('z');
 grid on; % 添加网格线
 axis equal; % 保持坐标轴比例一致
 view(3); % 设置三维视图

