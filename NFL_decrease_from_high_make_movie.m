clear; clc;

inp_fold = 'Z:\Siarhei Hladyshau\pheromone_hysteresis\1D_MCAS_type1\final_setup_2\updated_code6';
inp_file = 'NFL_decrease_from_high.mat';

load(fullfile(inp_fold, inp_file));

save_fold = inp_fold;
pic_fold = fullfile(save_fold, strcat(strrep(inp_file, '.mat', ''), '_phase_space_pictures'));
mkdir(pic_fold);

%phase_space_file = 'Z:\Siarhei Hladyshau\pheromone_hysteresis\1D_MCAS_type1\final_setup_2\k1_k2_phase_space\parameters_k1_0_10_k2_0_4_gradient.mat';
phase_space_file = 'Z:\Siarhei Hladyshau\pheromone_hysteresis\1D_MCAS_type1\final_setup_2\k1_k2_phase_space\parameters_k1_0_50_k2_0_4_gradient.mat';
load(phase_space_file, 'fin_state_1', 'fin_state_2', 'k1_vals', 'k2_vals');

k1_vals = k1_vals/SF/SF;

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

Ua_all = sol(:,:,1);
T = 30;
fig = figure('Position', [6 162 1446 614]);
for i = 1:T:size(sol,1)
    disp(i);
    
    if i < size(sol,1)/2
        idx = [1:T:i, i];
    else
        idx = [1:T:size(sol,1)/2, size(sol,1)/2:T:i, i];
    end
    
    Ua = sol(i,:,1);
    X_t = mean(sol(idx,:,3),2);
    s_t = s_vals(idx);
    
    Ua_vals = [Ua'; flip(Ua',1)];
    
    k1_all_t = k1+k6*s_t;
    k2_all_t = k2+k3*X_t;
    k1_coord = (k1_all_t-min(k1_vals))/(max(k1_vals)-min(k1_vals))*length(k1_vals);
    k2_coord = (k2_all_t-min(k2_vals))/(max(k2_vals)-min(k2_vals))*length(k2_vals);
    
    clf;
    x_ticks = linspace(0, length(Ua_vals), 5);
    x_ticklabels = {'5\pi','3.75\pi','2.5\pi','1.25\pi','0'};
    y_ticks = 0:5000:max(Ua_all(:));
    y_ticklabels = string([0:5000:max(Ua_all(:))]/10000);
    y_ticklabels(2:end) = strcat({y_ticklabels{2:end}}, "\cdot10^{4}");
    
    subplot(1,2,1);
    hold on;
    %grid on;
    box on;
    ylim([0, max(Ua_all(:))]);
    plot(Ua_vals, 'LineWidth', 3, 'Color', [1 0 0]);
    txt = {strcat('k1_{all}=', num2str(k1_all_t(end))), ...
        strcat('k2_{all}=', num2str(k2_all_t(end)))};
    xticks(x_ticks);
    xticklabels(x_ticklabels);
    yticks(y_ticks);
    yticklabels(y_ticklabels);
    xlabel('Distance $(\mu m)$', 'Interpreter','latex');
    ylabel('Active Cdc42 $\left(\frac{mol.}{\mu m} \right)$', 'Interpreter','latex');
    set(gca, 'FontSize', 15);
    set(gca,'linewidth',3);
    
    subplot(1,2,2);
    hold on;
    axis xy;
    imagesc(max2);
    plot(im1_c, im1_r, 'Color', [0 0 0], 'LineWidth', 5);
    plot(im1_c, im1_r, 'Color', [0 1 0], 'LineWidth', 3);
    plot(im2_c, im2_r, 'Color', [0 0 0], 'LineWidth', 5);
    plot(im2_c, im2_r, 'Color', [0 1 0], 'LineWidth', 3);
    colormap(hot(200));
    plot(k1_coord, k2_coord, 'Color', [0 0 0], 'LineStyle', '-', 'LineWidth', 5);
    plot(k1_coord, k2_coord, 'Color', [1 1 1], 'LineStyle', '-', 'LineWidth', 3);
    scatter(k1_coord(end), k2_coord(end), 'Marker', 'o', 'SizeData', 250, 'MarkerEdgeColor', [0 0 0], 'LineWidth', 2);
    scatter(k1_coord(end), k2_coord(end), 'Marker', 'o', 'SizeData', 250, 'MarkerEdgeColor', [1 1 1], 'LineWidth', 1);
    
    xlabel('$k_1 \left( \frac{\mu m^2}{s} \right)$', 'Interpreter','latex');
    x_ticks = 1:50:length(k1_vals);
    x_ticklabels = string(round(k1_vals(x_ticks)*10000,2));
    x_ticklabels(2:end) = strcat({x_ticklabels{2:end}}, "\cdot10^{-4}");
    
    xticks(x_ticks);
    xticklabels(x_ticklabels);
    xlim([0, length(k1_vals)]);
    xtickangle(0);
    
    ylabel('$k_2 \left( \frac{1}{s} \right)$', 'Interpreter','latex');
    y_ticks = 1:50:length(k2_vals);
    yticks(y_ticks);
    yticklabels(k2_vals(y_ticks));
    ylim([0, length(k2_vals)]);
    set(gca,'fontsize',15);
    
    saveas(fig, fullfile(pic_fold, strcat(num2str(i), '.png')));
    
end
mov_file = fullfile(save_fold, strrep(inp_file, '.mat', '_dynamics'));

imageNames = dir(fullfile(pic_fold, '*.png'));
outputVideo = VideoWriter(mov_file);
outputVideo.FrameRate = 24;
outputVideo.Quality = 100;
open(outputVideo);

for j = 1:T:size(sol,1)
    disp(j);
    img = imread(fullfile(pic_fold, strcat(num2str(j),'.png')));
    writeVideo(outputVideo,img);
end
close(outputVideo);