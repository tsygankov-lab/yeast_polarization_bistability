clear; clc;

phase_space_file = 'Z:\Siarhei Hladyshau\pheromone_hysteresis\1D_MCAS_type1\final_setup_2\k1_k2_phase_space\parameters_k1_0_50_k2_0_4_gradient.mat';

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
k1 = 0.22/SF/SF;
k2 = 0.5;
k3 = 7/SF;
k4 = 0.0001;
k5 = 1e-2;
k6 = 1/SF/SF/SF;

k4 = 0.0000001;
k3 = 7000/SF;

s1 = 48*SF;
s2 = s1;
%s3 = 0.5*SF;
s3 = 0*SF;

drop = false;
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

P = [x_len, t_max, k0, k1, k2, k3, k4, k5, k6, s1, s2, s3, drop, t_sig1, t_sig2, D_Ua, D_Ui, D_X, ...
     ex_U_start, ex_U_end, Ui_init, Ua_init_max, Ua_init_min, X_init];

%RelTol, AbsTol, NormControl, InitialStep, MaxStep
opts = [];

tic;
sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,opts,P);
toc;

Ua = sol(end,:,1);
Ui = sol(end,:,2);
X = sol(end,:,3);

s_vals = zeros(1, t_num);

for i = 1:t_num
    ti = t(i);
    if drop
        if (ti < t_sig1)
            s_vals(i) = s1 + (s2 - s1)*ti/t_sig1;
        end
        if (ti >= t_sig1) && (ti < t_sig2)
            s_vals(i) = s2;
        end
        if (ti >= t_sig2)
            s_vals(i) = s3;
        end
    else
        if (ti < t_sig1)
            s_vals(i) = s1 + (s2 - s1)*ti/t_sig1;
        end
        if (ti >= t_sig1) && (ti < t_sig2)
            s_vals(i) = s2;
        end
        if (ti >= t_sig2)
            s_vals(i) = s2 + (s3 - s2)*(ti-t_sig2)/(t_max-t_sig2);
        end
    end
end
X_vals = mean(sol(:,:,3),2);

Ua_vals = sol(:,:,1)';
Ua_vals = [Ua_vals; flip(Ua_vals,1)];
s = size(Ua_vals);

X_vals_all = zeros(1, size(sol(:,:,3),1));

for i = 1:size(sol(:,:,3),1)
    X_vals_all(i) = 2*sum(sol(i,:,3).*(x_len/x_num*ones(1,x_num)))*SF;
end

ids = 1:10:length(t);

fig1 = figure('Position',[50 30 1000 900]);

y_ticks = 0:0.00005:1.1*max(s_vals)*k6;
y_ticklabels = string(y_ticks*100000);
y_ticklabels(2) = strcat({y_ticklabels{2}}, "\cdot10^{-5}");
y_ticklabels(3) = "10^{-4}";

sp1 = subplot(3,1,1);
hold on;
grid off;
box on;
plot(t(ids), s_vals(ids)*k6, 'Color', [1 0 0], 'LineWidth', 4);
xticks(0:3600:t_max);
xticklabels((0:3600:t_max)/60);
ylabel('$s \left( \frac{\mu m^2}{s \cdot mol.^2} \right)$', 'Interpreter','latex');
set(gca,'fontsize',20);
ylim([0, 1.1*max(s_vals)*k6]);
yticks(y_ticks);
yticklabels(y_ticklabels);
set(gca,'linewidth',3);
xlim([0 t_max]);
set(sp1, 'InnerPosition', [0.1450 0.7093 0.74 0.1964]);

y_ticks = 0:1000:1.1*max(X_vals_all);
y_ticklabels = string(y_ticks/1000);
y_ticklabels(2) = "10^{3}";
y_ticklabels(3:end) = strcat({y_ticklabels{3:end}}, "\cdot10^{3}");
sp2 = subplot(3,1,2);
hold on;
grid off;
box on;
plot(t(ids), X_vals_all(ids), 'Color', [1 0 0], 'LineWidth', 4);
xticks(0:3600:t_max);
xticklabels((0:3600:t_max)/60);
ylabel('Y (molecules)');
title('Total Y','Color', [0,0,0]);
set(gca,'fontsize',20);
ylim([0, 1.1*max(X_vals_all)]);
yticks(y_ticks);
yticklabels(y_ticklabels);
set(gca,'linewidth',3);
xlim([0 t_max]);
set(sp2, 'InnerPosition', [0.1450 0.4096 0.74 0.1964]);

c_ticks = 0:5000:25*SF;
c_ticklabels = string([0:5000:25*SF]/1000);
c_ticklabels(2:end) = strcat({c_ticklabels{2:end}}, "\cdot10^{3}");
sp3 = subplot(3,1,3);
hold on;
axis ij;
imagesc(Ua_vals);
colormap(hot);
caxis([0 30*SF]);
set(gca, 'XLim', [0 s(2)], 'YLim', [0 s(1)]);
xlabel('Time (min)');
xticks((0:3600/(t_max/t_num):t_num));
xticklabels((0:3600:t_max)/60);
ylabel({'Distance (\mum)'});
yticks(linspace(0, 2*length(x),3));
%yticklabels(round(linspace(2*x_len,0,3), 2));
yticklabels({'5\pi','2.5\pi','0'}); 
title('Active Cdc42','Color', [0,0,0]);
set(gca,'fontsize',20);
cb = colorbar('Position', [0.89 0.1100 0.01 0.1967], ...
    'Ticks', c_ticks, ...
    'TickLabels', c_ticklabels);
set(sp3, 'InnerPosition', [0.1450 0.1100 0.74 0.1964]);

saveas(fig1, 'IFFL_decrease_from_high_dynamics.fig');
saveas(fig1, 'IFFL_decrease_from_high_dynamics.png');
set(fig1, 'PaperPosition', [0.25 0 8.25 8]);
saveas(fig1, 'IFFL_decrease_from_high_dynamics.pdf');

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
F3 = k4*s - k5*X;
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