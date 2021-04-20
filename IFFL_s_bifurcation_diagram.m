clear; clc;

% %radius of yeast cell in microns
% r_yc = 2.5;
% 
% m = 0;
% x_len = r_yc*pi;
% t_max = 10000000;
% x_num = 500;
% t_num = 5000;
% x = linspace(0,x_len,x_num);
% t = linspace(0,t_max,t_num);
% 
% k0 = 0.1;
% k1 = 0.1;
% k2 = 0.25;
% k3 = 0.7e-03;
% k4 = 0.2;
% k5 = 5e-04;
% 
% s1_vals = 0:0.0005:0.6;
% 
% drop = false;
% t_sig1 = t_max/3;
% t_sig2 = 2*t_max/3;
% 
% %diffusion coefficients in microns^2/s
% D_Ua = 0.0025;
% D_Ui = 0.25;
% D_X = 0.25;
% 
% X_init = 4.5900;
% 
% %initial conditions 1
% U_concentration = 1;
% U_all_1 = U_concentration*x_len;
% ex_U_start_1 = 0.9*x_len;
% ex_U_end_1 = x_len;
% Ua_init_max_1 = 0;
% Ua_init_min_1 = 0;
% Ui_init_1 = (U_all_1 - Ua_init_max_1*(ex_U_end_1-ex_U_start_1) - Ua_init_min_1*(x_len-ex_U_end_1+ex_U_start_1))/x_len;
% 
% %initial conditions 1
% U_all_2 = U_concentration*x_len;
% ex_U_start_2 = 0.9*x_len;
% ex_U_end_2 = x_len;
% Ua_init_max_2 = U_all_2/(ex_U_end_2-ex_U_start_2);
% Ua_init_min_2 = 0;
% Ui_init_2 = (U_all_2 - Ua_init_max_2*(ex_U_end_2-ex_U_start_2) - Ua_init_min_2*(x_len-ex_U_end_2+ex_U_start_2))/x_len;
% 
% fin_state_1 = zeros(length(s1_vals), x_num);
% fin_state_2 = zeros(length(s1_vals), x_num);
% 
% for i = 1:length(s1_vals)
%     disp([i, length(s1_vals)]);
% 
%     s1 = s1_vals(i);
%     s2 = s1;
%     s3 = s1;
%     
%     P1 = [x_len, t_max, k0, k1, k2, k3, k4, k5, s1, s2, s3, drop, t_sig1, t_sig2, D_Ua, D_Ui, D_X, ...
%           ex_U_start_1, ex_U_end_1, Ui_init_1, Ua_init_max_1, Ua_init_min_1, X_init];
%     
%     P2 = [x_len, t_max, k0, k1, k2, k3, k4, k5, s1, s2, s3, drop, t_sig1, t_sig2, D_Ua, D_Ui, D_X, ...
%           ex_U_start_2, ex_U_end_2, Ui_init_2, Ua_init_max_2, Ua_init_min_2, X_init];
%       
%     opts = ['MaxStep',0.01];
%     
%     sol1 = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,opts,P1);
%     sol2 = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,opts,P2);
%     
%     fin_state_1(i, :) = sol1(end,:,1);
%     fin_state_2(i, :) = sol2(end,:,1);
% end
% 
% save('IFFL_s_bifurcation_diagram.mat', 'r_yc', 'x_len', 't_max', 'x_num', 't_num', ...
%     'x', 't', 'k0', 'k1', 'k2', 'k3', 'k4', 'k5', 's1_vals', 'drop', ...
%     't_sig1', 't_sig2', 'D_Ua', 'D_Ui', 'D_X', 'X_init', 'U_concentration', ...
%     'U_all_1', 'ex_U_start_1', 'ex_U_end_1', 'Ua_init_max_1', 'Ua_init_min_1', ...
%     'Ui_init_1', 'U_all_2', 'ex_U_start_2', 'ex_U_end_2', 'Ua_init_max_2', ...
%     'Ua_init_min_2', 'Ui_init_2', 'fin_state_1', 'fin_state_2');

load('Z:\Siarhei Hladyshau\pheromone_hysteresis\1D_MCAS_type1\final_setup_2\IFFL\s_bifurc_diagram\s_bifurc_diagr_1.mat');

vals1 = zeros(1, size(fin_state_1,1));
vals2 = zeros(1, size(fin_state_1,1));

%scaling factor to convert concentrations from zeptomoles to molecules
SF = 10^(-21)*6.02*10^(23);

for i = 1:size(fin_state_1,1)
    vals1(i) = 2*sum(fin_state_1(i,:).*(x_len/x_num*ones(1,x_num)))*SF;
    vals2(i) = 2*sum(fin_state_2(i,:).*(x_len/x_num*ones(1,x_num)))*SF;
end

vals1(vals1 > (min(vals1)+max(vals1))/2) = NaN;
vals2(vals2 < 3) = NaN;

fig = figure('Position', [50 50 800 800]);
hold on;
grid off;
box on;
scatter(s1_vals/SF/SF, vals1, 'Marker', 'o', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [1, 0, 0]);
scatter(s1_vals/SF/SF, vals2, 'Marker', 'o', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [1, 0, 0]);
xlabel('Pheromone signal, $s \left( \frac{\mu m^2}{s \cdot mol.^2} \right)$', 'Interpreter','latex');
ylabel('Total ${C_A}$ (molecules)', 'Interpreter','latex');
set(gca, 'FontSize', 30);
xticks(0:0.0000007:1*max(s1_vals/SF/SF));
x_ticklabels = string((0:0.0000007:1*max(s1_vals/SF/SF))*1000000);
x_ticklabels(1) = "0^{ }";
x_ticklabels(2:end) = strcat({x_ticklabels{2:end}}, "\cdot10^{-6}");
xticklabels(x_ticklabels);
xlim([0, 1*max(s1_vals/SF/SF)]);
yticks(1000:1500:8500);
set(gca,'linewidth',4);

saveas(fig, 'IFFL_s_bifurcation_diagram_total_Ua.fig');
saveas(fig, 'IFFL_s_bifurcation_diagram_total_Ua.png');
saveas(fig, 'IFFL_s_bifurcation_diagram_total_Ua.pdf');


% --------------------------------------------------------------
function [c,f,s] = pdex4pde(x,t,u,DuDx,P)

t_max = P(2);
k0 = P(3);
k1 = P(4);
k2 = P(5);
k3 = P(6);
k4 = P(7);
k5 = P(8);
s1 = P(9);
s2 = P(10);
s3 = P(11);
drop = P(12);
t_sig1 = P(13);
t_sig2 = P(14);
D_Ua = P(15);
D_Ui = P(16);
D_X = P(17);

if drop
    if (t < t_sig1)
        s = s1 + (s2 - s1)*t/t_sig1 + 0.001*x;
    end
    if (t >= t_sig1) && (t < t_sig2)
        s = s2 + 0.001*x;
    end
    if (t >= t_sig2)
        s = s3 + 0.001*x;
    end
else
    if (t < t_sig1)
        s = s1 + (s2 - s1)*t/t_sig1 + 0.001*x;
    end
    if (t >= t_sig1) && (t < t_sig2)
        s = s2 + 0.001*x;
    end
    if (t >= t_sig2)
        s = s2 + (s3 - s2)*(t-t_sig2)/(t_max-t_sig2) + 0.001*x;
    end
end

Ua = u(1);
Ui = u(2);
X = u(3);

c = [1; 1; 1];
f = [D_Ua; D_Ui; D_X] .* DuDx; 
r = 0.001*x;

F1 = (k0 + r + (k1+s)*Ua^2)*Ui - (k2+k3*X)*Ua;
F2 = -F1;
F3 = k4*s - k5*X;
s = [F1; F2; F3];

end
% --------------------------------------------------------------
function u0 = pdex4ic(x,P)

ex_U_start = P(18);
ex_U_end = P(19);
Ui_init = P(20);
Ua_init_max = P(21);
Ua_init_min = P(22);
X_init = P(23);

if (x >= ex_U_start) && ( x <= ex_U_end) && (ex_U_start ~= ex_U_end)
    Ua_val = Ua_init_max;
else
    Ua_val = Ua_init_min;
end

r = 0;
u0 = [Ua_val+r; Ui_init-r; X_init];

end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t,P)
pl = [0;0;0];
ql = [1;1;1];
pr = [0;0;0];
qr = [1;1;1];
end