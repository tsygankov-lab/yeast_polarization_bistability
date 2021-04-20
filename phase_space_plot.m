clear; clc;
 
phase_space_file = 'Z:\Siarhei Hladyshau\pheromone_hysteresis\1D_MCAS_type1\final_setup_2\k1_k2_phase_space\parameters_k1_0_10_k2_0_4_gradient.mat';

load(phase_space_file, 'k1_vals', 'k2_vals', 'fin_state_1', 'fin_state_2');

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

fig2 = figure('Position', [50 50 800 850]);
hold on;
axis xy;
imagesc(max2);
plot(im1_c, im1_r, 'Color', [0 0 0], 'LineWidth', 10);
plot(im1_c, im1_r, 'Color', [0 1 0], 'LineWidth', 7);
plot(im2_c, im2_r, 'Color', [0 0 0], 'LineWidth', 10);
plot(im2_c, im2_r, 'Color', [0 1 0], 'LineWidth', 7);
colormap(hot);

%scaling factor to convert concentrations from zeptomoles to molecules
SF = 10^(-21)*6.02*10^(23);

xlabel('$k_1 \left( \frac{\mu m^2}{s \cdot mol.^2} \right)$', 'Interpreter','latex');
x_ticks = 1:90:length(k1_vals);
x_ticklabels = string(round(k1_vals(x_ticks)/SF/SF*100000,1));
x_ticklabels(2:end) = strcat({x_ticklabels{2:end}}, "\cdot10^{-5}");
xticks(x_ticks);
xticklabels(x_ticklabels);
xlim([0, length(k1_vals)]);

ylabel('$k_2 \left( \frac{1}{s} \right)$', 'Interpreter','latex');
y_ticks = 1:50:length(k2_vals);
yticks(y_ticks);
yticklabels(k2_vals(y_ticks));
ylim([0, length(k2_vals)]);

set(gca, 'FontSize', 30);

saveas(fig2, 'phase_space_outlines.png');
saveas(fig2, 'phase_space_outlines.fig');
saveas(fig2, 'phase_space_outlines.pdf');