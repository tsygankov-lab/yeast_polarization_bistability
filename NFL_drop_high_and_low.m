clear; clc;

phase_space_file = 'Z:\Siarhei Hladyshau\pheromone_hysteresis\1D_MCAS_type1\final_setup_2\k1_k2_phase_space\parameters_k1_0_2_k2_0_1_gradient.mat';

%radius of yeast cell in microns
r_yc = 2.5;

%scaling factor to convert concentrations from zeptomoles to molecules
SF = 10^(-21)*6.02*10^(23);

m = 0;
x_len = r_yc*pi;
t_max = 21600;
x_num = 500;
t_num = 10000;
x = linspace(0,x_len,x_num);
t = linspace(0,t_max,t_num);

k0 = 0.1;
k1 = 0.19/SF/SF;
k2 = 0.14;
k3 = 20/SF;
k4 = 0.0002;
k5 = 1.2e-03;
k6 = 1/SF/SF/SF;

k3 = 2000/SF;
k4 = 0.000002;

s1_high = 40*SF;
s2_high = s1_high;
s1_low = 2*SF;
s2_low = s1_low;
s3 = 0.5*SF;

drop = true;
t_sig1 = t_max/2;
t_sig2 = t_max/2;

%diffusion coefficients in microns^2/s
D_Ua = 0.0025;
D_Ui = 0.25;
D_X = 0.25;

%initial conditions 1
%concentrations are measured in zeptomoles
U_concentration = 1*SF;
U_all = U_concentration*x_len;
ex_U_start = 0.95*x_len;
ex_U_end = x_len;
Ua_init_max = 1*U_all/(ex_U_end-ex_U_start);
Ua_init_min = 0;
Ui_init = (U_all - Ua_init_max*(ex_U_end-ex_U_start) - Ua_init_min*(x_len-ex_U_end+ex_U_start))/x_len;
X_init = 0;

P_high = [x_len, t_max, k0, k1, k2, k3, k4, k5, k6, s1_high, s2_high, s3, drop, t_sig1, t_sig2, D_Ua, D_Ui, D_X, ...
     ex_U_start, ex_U_end, Ui_init, Ua_init_max, Ua_init_min, X_init];
P_low = [x_len, t_max, k0, k1, k2, k3, k4, k5, k6, s1_low, s2_low, s3, drop, t_sig1, t_sig2, D_Ua, D_Ui, D_X, ...
     ex_U_start, ex_U_end, Ui_init, Ua_init_max, Ua_init_min, X_init]; 

%RelTol, AbsTol, NormControl, InitialStep, MaxStep
opts = [];

tic;
sol_high = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,opts,P_high);
sol_low = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,opts,P_low);
toc;

% Ua = sol(end,:,1);
% Ui = sol(end,:,2);
% X = sol(end,:,3);

s_vals_high = zeros(1, t_num);
s_vals_low = zeros(1, t_num);

for i = 1:t_num
    ti = t(i);
    if drop
        if (ti < t_sig1)
            s_vals_high(i) = s1_high + (s2_high - s1_high)*ti/t_sig1;
            s_vals_low(i) = s1_low + (s2_low - s1_low)*ti/t_sig1;
        end
        if (ti >= t_sig1) && (ti < t_sig2)
            s_vals_high(i) = s2_high;
            s_vals_low(i) = s2_low;
        end
        if (ti >= t_sig2)
            s_vals_high(i) = s3;
            s_vals_low(i) = s3;
        end
    else
        if (ti < t_sig1)
            s_vals_high(i) = s1_high + (s2_high - s1_high)*ti/t_sig1;
            s_vals_low(i) = s1_low + (s2_low - s1_low)*ti/t_sig1;
        end
        if (ti >= t_sig1) && (ti < t_sig2)
            s_vals_high(i) = s2_high;
            s_vals_low(i) = s2_low;
        end
        if (ti >= t_sig2)
            s_vals_high(i) = s2_high + (s3 - s2_high)*(ti-t_sig2)/(t_max-t_sig2);
            s_vals_low(i) = s2_low + (s3 - s2_low)*(ti-t_sig2)/(t_max-t_sig2);
        end
    end
end

X_vals_high = mean(sol_high(:,:,3),2);
X_vals_low = mean(sol_low(:,:,3),2);

Ua_vals_high = sol_high(:,:,1)';
Ua_vals_high = [Ua_vals_high; flip(Ua_vals_high,1)];
Ua_vals_low = sol_low(:,:,1)';
Ua_vals_low = [Ua_vals_low; flip(Ua_vals_low,1)];

s = size(Ua_vals_high);

X_vals_high_all = zeros(1, size(sol_high(:,:,3),1));
for i = 1:size(sol_high(:,:,3),1)
    X_vals_high_all(i) = 2*sum(sol_high(i,:,3).*(x_len/x_num*ones(1,x_num)))*SF;
end

X_vals_low_all = zeros(1, size(sol_low(:,:,3),1));
for i = 1:size(sol_low(:,:,3),1)
    X_vals_low_all(i) = 2*sum(sol_low(i,:,3).*(x_len/x_num*ones(1,x_num)))*SF;
end

disp([num2str(max(X_vals_high_all)), '  ', num2str(max(X_vals_low_all))]);

ids = 1:10:length(t);

fig = figure('Position',[50 30 600 900]);

y_ticks = 0:0.00005:1.1*max(s_vals_high)*k6;
y_ticklabels = string(y_ticks*100000);
y_ticklabels(2) = strcat({y_ticklabels{2}}, "\cdot10^{-5}");
y_ticklabels(3) = "10^{-4}";

sp1 = subplot(4,1,1);
hold on;
grid off;
box on;
plot(t(ids), s_vals_high(ids)*k6, 'Color', [0 0 0], 'LineWidth', 4);
plot(t(ids), s_vals_low(ids)*k6, 'Color', [1 0 0], 'LineWidth', 4, 'LineStyle', ':');
xlim([0 t_max]);
xticks(0:3600:t_max);
xticklabels((0:3600:t_max)/60);
yticks(y_ticks);
yticklabels(y_ticklabels);
ylim([0, 1.1*max(s_vals_high)*k6]);
ylabel('$s \left( \frac{\mu m^2}{s \cdot mol.^2} \right)$', 'Interpreter','latex');
set(gca,'fontsize',15);
set(gca,'linewidth',2);
%set(sp1, 'InnerPosition', [0.1933 0.6061 0.65 0.2924]);

sp2 = subplot(4,1,2);
hold on;
axis ij;
imagesc(Ua_vals_high);
colormap(hot);
caxis([0 30*SF]);
set(gca, 'XLim', [0 s(2)], 'YLim', [0 s(1)]);
xticks((0:3600/(t_max/t_num):t_num));
xticklabels((0:3600:t_max)/60);
ylabel({'Distance (\mum)'});
yticks(linspace(0, 2*length(x),3));
yticklabels({'5\pi','2.5\pi','0'}); 
text(3350, 150, 'Active Cdc42', 'FontSize', 18, 'Color', [1 1 1], 'FontWeight', 'Bold');
set(gca,'fontsize',15);
%set(sp2, 'InnerPosition', [0.1933 0.1322 0.65 0.2924]);

sp3 = subplot(4,1,3);
hold on;
axis ij;
imagesc(Ua_vals_low);
colormap(hot);
caxis([0 30*SF]);
set(gca, 'XLim', [0 s(2)], 'YLim', [0 s(1)]);
xlabel('Time (min)');
xticks((0:3600/(t_max/t_num):t_num));
xticklabels((0:3600:t_max)/60);
ylabel({'Distance (\mum)'});
yticks(linspace(0, 2*length(x),3));
yticklabels({'5\pi','2.5\pi','0'}); 
%title('Active Cdc42','Color', [0,0,0]);
text(3350, 150, 'Active Cdc42', 'FontSize', 18, 'Color', [1 1 1], 'FontWeight', 'Bold');
set(gca,'fontsize',15);

c_ticks = 0:5000:30*SF;
c_ticklabels = string([0:5000:30*SF]/1000);
c_ticklabels(2:end) = strcat({c_ticklabels{2:end}}, "\cdot10^{3}");
cb = colorbar('Position', [0.28 0.25 0.5 0.015], ...
    'Ticks', c_ticks, ...
    'TickLabels', c_ticklabels, 'Location', 'SouthOutside', ...
    'LineWidth', 1);

saveas(fig, 'NFL_drop_from_high_and_low.png');
saveas(fig, 'NFL_drop_from_high_and_low.pdf');
saveas(fig, 'NFL_drop_from_high_and_low.fig');

save('NFL_drop_from_high_and_low.mat', ...
    'r_yc', 'SF', 'm', 'x_len', 't_max', 'x_num', 't_num', 'x', 't', ...
    'k0', 'k1', 'k2', 'k3', 'k4', 'k5', 'k6', 's1_high', 's2_high', ...
    's1_low', 's2_low', 's3', 'drop', 't_sig1', 't_sig2', 'D_Ua', ...
    'D_Ui', 'D_X', 'U_concentration', 'U_all', 'ex_U_start', 'ex_U_end', ...
    'Ua_init_max', 'Ua_init_min', 'Ui_init', 'X_init', 'P_high', 'P_low', ...
    'opts', 'sol_high', 'sol_low', 's_vals_high', 's_vals_low', 'X_vals_high', ...
    'X_vals_low', 'Ua_vals_high', 'Ua_vals_low', 's', 'X_vals_high_all', ...
    'X_vals_low_all');


% --------------------------------------------------------------
function [c,f,s] = pdex4pde(x,t,u,DuDx,P)

t_max = P(2);
k0 = P(3);
k1 = P(4);
k2 = P(5);
k3 = P(6);
k4 = P(7);
k5 = P(8);
k6 = P(9);
s1 = P(10);
s2 = P(11);
s3 = P(12);
drop = P(13);
t_sig1 = P(14);
t_sig2 = P(15);
D_Ua = P(16);
D_Ui = P(17);
D_X = P(18);

if drop
    if (t < t_sig1)
        s = s1 + (s2 - s1)*t/t_sig1;
    end
    if (t >= t_sig1) && (t < t_sig2)
        s = s2;
    end
    if (t >= t_sig2)
        s = s3;
    end
else
    if (t < t_sig1)
        s = s1 + (s2 - s1)*t/t_sig1;
    end
    if (t >= t_sig1) && (t < t_sig2)
        s = s2;
    end
    if (t >= t_sig2)
        s = s2 + (s3 - s2)*(t-t_sig2)/(t_max-t_sig2);
    end
end

Ua = u(1);
Ui = u(2);
X = u(3);

c = [1; 1; 1];
f = [D_Ua; D_Ui; D_X] .* DuDx; 
r = 0.001*x;

F1 = (k0 + r + (k1 + k6*s*(1+r))*Ua^2)*Ui - (k2+k3*X)*Ua;
F2 = -F1;
F3 = k4*Ua - k5*X;
s = [F1; F2; F3];

end
% --------------------------------------------------------------
function u0 = pdex4ic(x,P)

ex_U_start = P(19);
ex_U_end = P(20);
Ui_init = P(21);
Ua_init_max = P(22);
Ua_init_min = P(23);
X_init = P(24);

if (x >= ex_U_start) && ( x <= ex_U_end) && (ex_U_start ~= ex_U_end)
    Ua_val = Ua_init_max;
else
    Ua_val = Ua_init_min;
end

r = 0.03*abs(randn(1,1));
u0 = [Ua_val+r; Ui_init-r; X_init];

end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t,P)
pl = [0;0;0];
ql = [1;1;1];
pr = [0;0;0];
qr = [1;1;1];
end
