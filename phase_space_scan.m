clear; clc;

%radius of yeast cell in microns
r_yc = 2.5;

m = 0;
x_len = r_yc*pi;
t_max = 1000000;
x_num = 500;
t_num = 1000;
x = linspace(0,x_len,x_num);
t = linspace(0,t_max,t_num);

k0 = 0.1;

%diffusion coefficients in microns^2/s
D_Ua = 0.0025;
D_Ui = 0.25;

%initial conditions 1
U_concentration = 1;
U_all_1 = U_concentration*x_len;
ex_U_start_1 = 0;
ex_U_end_1 = 0.1;
Ua_init_max_1 = 0;
Ua_init_min_1 = 0;
Ui_init_1 = (U_all_1 - Ua_init_max_1*(ex_U_end_1-ex_U_start_1) - Ua_init_min_1*(x_len-ex_U_end_1+ex_U_start_1))/x_len;

%initial conditions 1
U_all_2 = U_concentration*x_len;
ex_U_start_2 = 0;
ex_U_end_2 = 0.1;
Ua_init_max_2 = U_all_2/(ex_U_end_2-ex_U_start_2);
Ua_init_min_2 = 0;
Ui_init_2 = (U_all_2 - Ua_init_max_2*(ex_U_end_2-ex_U_start_2) - Ua_init_min_2*(x_len-ex_U_end_2+ex_U_start_2))/x_len;

%k1_vals = 0:0.05:2;
%k2_vals = 0:0.02:1;

%k1_vals = 0:0.02:2;
%k2_vals = 0:0.01:1;

%k1_vals = 0:0.05:10;
%k2_vals = 0:0.01:2;

%k1_vals = 0:0.01:2;
%k2_vals = 0:0.01:2;

k1_vals = 0:0.05:10;
k2_vals = 0:0.02:4;

% k1_vals = 0:0.25:50;
% k2_vals = 0:0.02:4;

fin_state_1 = zeros(length(k1_vals), length(k2_vals), x_num);
fin_state_2 = zeros(length(k1_vals), length(k2_vals), x_num);

for i = 1:length(k1_vals)
    for j = 1:length(k2_vals)
        disp([i, length(k1_vals), j, length(k2_vals)]);
        
        k1 = k1_vals(i);
        k2 = k2_vals(j);
        
        P1 = [x_len, k0, k1, k2, ...
            D_Ua, D_Ui, ...
            ex_U_start_1, ex_U_end_1, Ui_init_1, Ua_init_max_1, Ua_init_min_1, ...
            t_max];
        
        P2 = [x_len, k0, k1, k2, ...
            D_Ua, D_Ui, ...
            ex_U_start_2, ex_U_end_2, Ui_init_2, Ua_init_max_2, Ua_init_min_2, ...
            t_max];

        opts = ['MaxStep',0.01];
        
        sol1 = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,opts,P1);
        sol2 = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,opts,P2);
        
        fin_state_1(i, j, :) = sol1(end,:,1);
        fin_state_2(i, j, :) = sol2(end,:,1);
    end
end

pars_file = strcat('phase_space_scan_k1_', num2str(min(k1_vals)), '_', num2str(max(k1_vals)), ...
    '_k2_', num2str(min(k2_vals)), '_', num2str(max(k2_vals)), '_gradient.mat');

save(pars_file, ...
    'm', 'x_len', 't_max', 'x_num', 't_num', 'x', 't', 'k0', 'k1_vals', 'k2_vals', ...
    'D_Ua', 'D_Ui', 'U_all_1', 'ex_U_start_1', 'ex_U_end_1', 'Ua_init_max_1', 'Ua_init_min_1', ...
    'Ui_init_1', 'U_all_2', 'ex_U_start_2', 'ex_U_end_2', 'Ua_init_max_2', 'Ua_init_min_2', ...
    'Ui_init_2', 'fin_state_1', 'fin_state_2');

max1 = zeros(length(k1_vals), length(k2_vals));
max2 = zeros(length(k1_vals), length(k2_vals));

for i = 1:length(k1_vals)
    for j = 1:length(k2_vals)
        max1(j,i) = max(fin_state_1(i,j,:));
        max2(j,i) = max(fin_state_2(i,j,:));
    end
end

thr1 = 0.99;
im1 = double(max1 < thr1);
im1_outline = cell_outline(im1);
[im1_r, im1_c] = find(im1_outline);

thr2 = 0.99;
im2 = double(max2 < thr2);
im2_outline = cell_outline(im2);
[im2_r, im2_c] = find(im2_outline);

fig2 = figure('Position', [6 1 899 816]);
hold on;
axis xy;
imagesc(max2);
plot(im1_c, im1_r, 'Color', [0 0 0], 'LineWidth', 10);
plot(im1_c, im1_r, 'Color', [0 1 0], 'LineWidth', 7);
plot(im2_c, im2_r, 'Color', [0 0 0], 'LineWidth', 10);
plot(im2_c, im2_r, 'Color', [0 1 0], 'LineWidth', 7);
colormap(hot);

xlabel('k1 (\mum^2/s)');
x_ticks = 1:20:length(k1_vals);
xticks(x_ticks);
xticklabels(k1_vals(x_ticks));
xlim([0, length(k1_vals)]);
xtickangle(90);

ylabel('k2 (1/s)');
y_ticks = 1:20:length(k2_vals);
yticks(y_ticks);
yticklabels(k2_vals(y_ticks));
ylim([0, length(k2_vals)]);

set(gca, 'FontSize', 25);

saveas(fig2, strrep(pars_file, '.mat', '_hot_outlines.png'));

% --------------------------------------------------------------
function [c,f,s] = pdex4pde(x,t,u,DuDx,P)

k0 = P(2);
k1 = P(3);
k2 = P(4);
D_Ua = P(5);
D_Ui = P(6);

Ua = u(1);
Ui = u(2);

r = 0.001*x;

c = [1; 1];
f = [D_Ua; D_Ui] .* DuDx; 
F1 = (k0 + r + k1*Ua^2)*Ui - k2*Ua; 
F2 = -F1;
s = [F1; F2];

end
% --------------------------------------------------------------
function u0 = pdex4ic(x,P)

ex_U_start = P(7);
ex_U_end = P(8);
Ui_init = P(9);
Ua_init_max = P(10);
Ua_init_min = P(11);

if (x >= ex_U_start) && ( x <= ex_U_end) && (ex_U_start ~= ex_U_end)
    Ua_val = Ua_init_max;
else
    Ua_val = Ua_init_min;
end

r = 0.01*abs(randn(1,1));

u0 = [Ua_val+r;Ui_init-r];

end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t,P)
pl = [0;0];
ql = [1;1];
pr = [0;0];
qr = [1;1];
end