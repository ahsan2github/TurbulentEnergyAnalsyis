parentDir = '/raid/home/mohammad/Ekman_Ugeo2_2048x2048x64_Lz750';
% parentDir = '/raid2/mohammad/ek02_768x768x96';
outputDir = [parentDir '/output/'];
readinputs(parentDir);
u_star_in = u_star; clear u_star;
BLH = 550.69;
ustar = 0.09;
refLevels = 3:3:21;  % total 21 levels
nscale = log(Nx)/log(2);
data = load('wav_pi_tmn_ek02_fixed.mat')
zvec = (0:Nz-1).*dz;
epsilon = ustar^3./dz/0.4;
pi_energy = data.frame_energy_pi;
pi_energy2 = data.frame_energy_pi2;
LineStyles = {'-', '--','-.',':','-+','-o','-s','-d','-^','-*'}; % 10 types total
figWidth = 3.5; figHeight = figWidth;  formatEng = '-depsc2';

%% extract pi^(n) at different heights for all scales
n_st = 3; n_pos_st = n_st-1; n_end = 10; n_pos_end = n_end-1;
heights = zeros([numel(pi_energy{1}) 1]);
k_indx = zeros([numel(pi_energy{1}) 1]);
for ii=1:length(heights)
    heights(ii) = pi_energy{1}{ii}.h;
    k_indx(ii) = int32(pi_energy{1}{ii}.h/dz);
end
pi_n_data_frame = zeros([numel(n_pos_st:n_pos_end) numel(refLevels) length(pi_energy)]);
pi_n_data_sigma_frame = zeros([numel(n_pos_st:n_pos_end) numel(refLevels) length(pi_energy)]);
pi_n_data_frame2 = zeros([numel(n_pos_st:n_pos_end) numel(refLevels) length(pi_energy)]);
pi_n_data_sigma_frame2 = zeros([numel(n_pos_st:n_pos_end) numel(refLevels) length(pi_energy)]);
for ff = 1:length(pi_energy2)
    for jj = 1:length(pi_energy2{1}) 
        count=0;
        for ii = n_pos_st:n_pos_end
            pi_n_data_frame2(count+1,jj,ff) = -pi_energy2{ff}{jj}.pi_energy(count+n_pos_st)...
                /epsilon/2^(2*nscale-2*pi_energy2{ff}{jj}.n(ii));
            pi_n_data_sigma_frame2(count+1,jj,ff) = pi_energy2{ff}{jj}.pi_energy_sigma(count+n_pos_st)/epsilon;
            count = count+1;
        end
    end
end
pi_n_data2 = mean(pi_n_data_frame2,3); % average over frames
pi_n_data_sigma2 = mean(pi_n_data_sigma_frame2,3); % average over frames
n_scale_in_m2 = zeros([size(pi_n_data2,1) 1]);
count=0;
for ii = n_pos_st:n_pos_end
    n_scale_in_m2(count+1)= pi_energy2{1}{1}.n_in_m(count+n_pos_st);
    count = count+1;
end

for ff = 1:length(pi_energy)
    for jj = 1:length(pi_energy{1}) 
        count=0;
        for ii = n_pos_st:n_pos_end
            pi_n_data_frame(count+1,jj,ff) = -pi_energy{ff}{jj}.pi_energy(count+n_pos_st)...
                /epsilon/2^(2*nscale-2*pi_energy{ff}{jj}.n(ii));
            pi_n_data_sigma_frame(count+1,jj,ff) = pi_energy{ff}{jj}.pi_energy_sigma(count+n_pos_st)/epsilon;
            count = count+1;
        end
    end
end
pi_n_data = mean(pi_n_data_frame,3); % average over frames
pi_n_data_sigma = mean(pi_n_data_sigma_frame,3); % average over frames
n_scale_in_m = zeros([size(pi_n_data,1) 1]);
count=0;
for ii = n_pos_st:n_pos_end
    n_scale_in_m(count+1)= pi_energy{1}{1}.n_in_m(count+n_pos_st);
    count = count+1;
end

% plot pi at different heights for all n

% %plot a fixed n vs all m at different heights
% figh = figure(); 
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figheight]);
% bounds = zeros([numel(pi_energy{1}.n) 2  numel(refLevels)]);
% for ii = 1:size(pi_n_data,1)
%     for mmm = 1:size(pi_n_data,2)
%         bounds(ii,:,mmm) = [pi_energy{ii}.pi_energy_sigma(mmm) pi_energy{ii}.pi_energy_sigma(mmm)];
%     end
% end
% [l,p]= boundedline(heights./BLH, -pi_n_data(:,:),...
%         bounds./ustar^3.*0.5,'cmap', parula(8), 'transparency', 0.3);

hf = figure(); hold all;
set(hf, 'Visible','off');
box on;
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
set(gcf,'Renderer','Painters');

for ii = 1:size(pi_n_data,1)
    plot(heights./BLH,pi_n_data(ii,:),LineStyles{ii},'MarkerSize',4);
end
xlabel('$z/\delta$','Interpreter','Latex'); 
ylabel('$\left < \pi^{n}_{sg}(z)\right >/(u_*^3\Delta z^{-1}\kappa^{-1})$','Interpreter','Latex');
lgndVal = cell([numel(n_pos_st:n_pos_end) 1]);
for ii = 1:size(pi_n_data,1)
    lgndVal{ii} = sprintf('r_n/\\delta = %2.2f',pi_energy{1}{1}.n_in_m(n_pos_st+ii-1)/BLH);
end
lgndh = legendflex(lgndVal, 'nrow', 4,'ncol',2,'box','off','ref', gca,...
        'anchor', [3 3], 'buffer', [-10 -10],'fontsize',7);
xlim([0 1.0]);
ylim([-1 4]);
set(gca,...
'XMinorTick',   'on',...
'LineWidth',    0.5);
set(gcf, 'Color', 'w');
SPR = sprintf('%s%s','','pi_ek02_diff_n_by_u3_dz.eps');  
print(formatEng,SPR);

hf = figure(); hold all;
set(hf, 'Visible','off');
box on;
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
set(gcf,'Renderer','Painters');
for ii = 1:size(pi_n_data,1)
    plot(heights./BLH,pi_n_data2(ii,:),LineStyles{ii},'MarkerSize',4);
end
xlabel('$z/\delta$','Interpreter','Latex'); 
ylabel('$\left < \pi^{\prime \ n}_{sg}(z)\right >/(u_*^3\Delta z^{-1}\kappa^{-1})$','Interpreter','Latex');
lgndVal = cell([numel(n_pos_st:n_pos_end) 1]);
for ii = 1:size(pi_n_data,1)
    lgndVal{ii} = sprintf('r_n/\\delta = %2.2f',pi_energy2{1}{1}.n_in_m(n_pos_st+ii-1)/BLH);
end
lgndh = legendflex(lgndVal, 'nrow', 4,'ncol',2,'box','off','ref', gca,...
        'anchor', [3 3], 'buffer', [-10 -10],'fontsize',7);
xlim([0 1.0]);
ylim([-1 4]);
set(gca,...
'XMinorTick',   'on',...
'LineWidth',    0.5);
set(gcf, 'Color', 'w');
SPR = sprintf('%s%s','','pi_ek02_diff_n_by_u3_dz_2.eps');  
print(formatEng,SPR);
%%
%
% text_x_pos = 0.4;
% text_y_pos = 0.15;
% xlim_max = 12;
% figure(); 
% set(gcf,'Renderer','Painters');
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0.1 figWidth*1.0 figHeight*1.0]);
% set(gcf,'PaperPositionMode','auto');
% subplot(231); hold all;
% select_h = 1;
% plot(n_scale_in_m./BLH, pi_n_data(:,select_h),'-om');
% % plot(n_scale_in_m./BLH, pi_n_data_sigma(:,select_h),'+c');
% % plot(n_scale_in_m./BLH, -pi_n_data_sigma(:,select_h),'sc');
% ann_str = [sprintf('z=%2.2f',heights(select_h)/BLH) sprintf('\\delta')];
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
% % ylabel('$\left < \pi^{n}_{sg}\right >(z)/\left < \epsilon \right>$','Interpreter','Latex');
% ylabel('$\left < \pi^{n}_{sg}(z)\right >/(u_*^3\Delta z^{-1}\kappa^{-1})$','Interpreter','Latex');
% % ylabel('$\left < \pi^{n}_{sg}\right >(z)/\left < \left < u_iu_j \right> \frac{\partial U_i}{\partial x_j}\right >$','Interpreter','Latex');
% xlim([-2 xlim_max]);
% box on;
% subplot(232)
% select_h = 3;
% plot(n_scale_in_m./BLH, -pi_n_data(:,select_h),'-om');
% ann_str = [sprintf('z=%2.2f',heights(select_h)/BLH) sprintf('\\delta')];
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
% xlim([-2 xlim_max]);
% 
% subplot(233)
% select_h = 5;
% plot(n_scale_in_m./BLH, pi_n_data(:,select_h),'-om');
% %xlabel('$2^{n}h/\delta$','Interpreter','Latex'); 
% %ylabel('$\left < \pi^{n}_{sg}\right >$','Interpreter','Latex');
% ann_str = [sprintf('z=%2.2f',heights(select_h)/BLH) sprintf('\\delta')];
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
% xlim([-2 xlim_max]);
% 
% subplot(234)
% select_h = 10;
% plot(n_scale_in_m./BLH, pi_n_data(:,select_h),'-om');
% xlabel('$r_n/\delta$','Interpreter','Latex'); 
% % ylabel('$\left < \pi^{n}_{sg}\right >(z)/\left < \epsilon \right >$','Interpreter','Latex');
% ylabel('$\left < \pi^{n}_{sg}(z)\right >/(u_*^3\Delta z^{-1}\kappa^{-1})$','Interpreter','Latex');
% % ylabel('$\left < \pi^{n}_{sg}\right >(z)/\left < \left < u_iu_j \right> \frac{\partial U_i}{\partial x_j}\right >$','Interpreter','Latex');
% ann_str = [sprintf('z=%2.2f',heights(select_h)/BLH) sprintf('\\delta')];
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
% xlim([-2 xlim_max]);
% 
% subplot(235)
% select_h = 18;
% plot(n_scale_in_m./BLH, -pi_n_data(:,select_h),'-om');
% xlabel('$r_n/\delta$','Interpreter','Latex'); 
% %ylabel('$\left < \pi^{n}_{sg}\right >$','Interpreter','Latex');
% ann_str = [sprintf('z=%2.2f',heights(select_h)/BLH) sprintf('\\delta')];
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
% xlim([-2 xlim_max]);
% 
% subplot(236)
% select_h = 20;
% plot(n_scale_in_m./BLH, pi_n_data(:,select_h),'-om');
% xlabel('$r_n/\delta$','Interpreter','Latex'); 
% %ylabel('$\left < \pi^{n}_{sg}\right >$','Interpreter','Latex');
% ann_str = [sprintf('z=%2.2f',heights(select_h)/BLH) sprintf('\\delta')];
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
% xlim([-2 xlim_max]);
% 
% SPR = sprintf('%s%s','','pi_ek02_fixed_h_by_u3_dz.eps');  
% set(gcf, 'Color', 'w');
% %print(formatEng,SPR);
% 
% %% what happens to each scale on a volume average sense
% pi_scale_vol_avg = zeros([1 numel(n_pos_st:n_pos_end)]);
% for kk = 1:length(pi_scale_vol_avg)
%     pi_scale_vol_avg(kk) = trapz(heights,pi_n_data(kk,:))/BLH;
% end
%figure();
%plot(n_scale_in_m./BLH, -pi_scale_vol_avg,'-om');
%xlabel('$2^{n}h/\delta$','Interpreter','Latex'); 

