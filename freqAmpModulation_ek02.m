
% This analysis is based on the paper Amplitude and Frequency modulation
% in wall turbulence, B. Ganapathisubramani, N. Hutchins, J. P. Monty,
% D. Chung & I. Marusic, JFM 712(2012)
clear; clc; close all;
%addpath('/raid/home/mohammad/rwt-master/bin');
frameVec= [38,39,40,41]; 
parentDir = '/media/ahsan/ds00/Ekman_Ugeo2_2048x2048x64_Lz750';
outputDir = [parentDir '/output/'];
readinputs(parentDir);
u_star_in = u_star; clear u_star;
BLH = 550.69;
u_star = 0.09;
heightVec = ((1:Nz-17)-1).*dz/BLH;
NN = 2; % NN\delta will be the filter length used to separate large and small scale fluctuation
u_l = zeros([Nx Ny length(heightVec) length(frameVec)]);
u_s = zeros([Nx Ny length(heightVec) length(frameVec)]);
uwsig = zeros([Nx Ny length(heightVec) length(frameVec)]);

for ff = 1:length(frameVec)
    % load frames
    frameStr = sprintf('%4.4i',frameVec(ff));
    fn = [parentDir '/output/u_frame/u_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    uu = fread(fh,'double');
    uu = reshape(uu, [Nx,Ny,Nz]).*u_star_in + Ugal; 
    uu = interpolate_on_w(uu);
    
    sfl = '/media/ahsan/ds01/stressFiles';
    fn = [sfl '/ek02_tau13_'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    sgs_tau13 = fread(fh,'double');
    sgs_tau13 = reshape(sgs_tau13, [Nx,Ny,Nz]);
    fclose(fh);     
    fn = [sfl '/ek02_uw_'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    sgs_uw = fread(fh,'double');
    sgs_uw = reshape(sgs_uw, [Nx,Ny,Nz]);
    fclose(fh);     
    total_uw = sgs_tau13 + sgs_uw;
    
    for kk =1:1:length(heightVec)
        tmp = uu(:,:,kk)-mean(mean(uu(:,:,kk)));
        uufilt = cuttoff_2D_multiF(tmp,2*pi*z_i/(NN*BLH));
        uufluc = tmp-uufilt(:,:); 
        uw_plane = total_uw(:,:,kk);
        for jj=1:Ny
            u_l(:,jj,kk,ff) = uufilt(:,jj); %mean(uufilt,2);
            u_s(:,jj,kk,ff) = uufluc(:,jj); %mean(uufluc,2);
            uwsig(:,jj,kk,ff) = uw_plane(:,jj);
        end
    end
end

%% sample u_l and get varaince of u_s at each 2\delta segment
horzBinWidth = round(BLH*NN/dx);
horzEdges = (1:horzBinWidth:Nx);
bin_midVal_ul = zeros([length(horzEdges)-1 Ny length(heightVec) length(frameVec)]);
bin_variance_us = zeros([length(horzEdges)-1  Ny length(heightVec) length(frameVec)]);
bin_std_uw = zeros([length(horzEdges)-1  Ny length(heightVec) length(frameVec)]);
bin_mean_uw = zeros([length(horzEdges)-1  Ny length(heightVec) length(frameVec)]);
bin_waveLen_us = zeros([length(horzEdges)-1  Ny length(heightVec) length(frameVec)]);
bin_zeroCrossing_us = zeros([length(horzEdges)-1  Ny length(heightVec) length(frameVec)]);
bin_zeroCrossing_uw = zeros([length(horzEdges)-1  Ny length(heightVec) length(frameVec)]);
waveLen = 2*BLH./(1:round(Nx*NN*BLH/(2*pi*z_i))/2);
%wave no = 2*pi*z_i/waveLen
for ff=1:length(frameVec)
    for kk = 1:length(heightVec)
        for jj =1:Ny
            for ii = 1:length(horzEdges)-1
                bin_midVal_ul(ii,jj,kk,ff) = u_l(horzEdges(ii)+round(horzBinWidth/2),jj,kk,ff);
                bin_variance_us(ii,jj,kk,ff) = var(u_s(horzEdges(ii)+1:horzEdges(ii+1)-1,jj,kk,ff));
                bin_std_uw(ii,jj,kk,ff) = std(uwsig(horzEdges(ii)+1:horzEdges(ii+1)-1,jj,kk,ff));
                bin_mean_uw(ii,jj,kk,ff) = mean(uwsig(horzEdges(ii)+1:horzEdges(ii+1)-1,jj,kk,ff));
                tmp_sig = u_s(horzEdges(ii)+1:horzEdges(ii+1)-1,jj,kk,ff);
                tmp = abs(fft(tmp_sig))/ length(tmp_sig);
                [~,indx] = max(tmp(2:floor(numel(tmp)/2)+1));
                bin_waveLen_us(ii,jj,kk,ff) = waveLen(indx); 
                bin_zeroCrossing_us(ii,jj,kk,ff) = hzerocross(u_s(horzEdges(ii)+1:horzEdges(ii+1)-1,jj,kk,ff))/2; 
                bin_zeroCrossing_uw(ii,jj,kk,ff) = hzerocross(uwsig(horzEdges(ii)+1:horzEdges(ii+1)-1,jj,kk,ff))/2;
            end
        end
    end
end
horz_bin_sample_size = numel(u_s(horzEdges(1)+1:horzEdges(1+1)-1,1,1,1));
disp(['Horizontal Binning Done !']);
% now do data binning for sampled u_l values and keep track where each sampled values go
ul_binwidth = 0.01;
midval_ul = -0.65:ul_binwidth:0.65;
ul_edges=[midval_ul(1)-ul_binwidth/2:ul_binwidth:midval_ul(end)+ul_binwidth/2];
bin_midVal_ul_binIndx = zeros([length(horzEdges)-1 Ny length(heightVec) length(frameVec)]);
bin_midVal_waveLen_binIndx = zeros([length(horzEdges)-1 Ny length(heightVec) length(frameVec)]);
bin_midVal_zeroCrossing_binIndx = zeros([length(horzEdges)-1 Ny length(heightVec) length(frameVec)]);
bin_midVal_ul_binFreq = zeros([length(ul_edges)-1 Ny length(heightVec) length(frameVec)]);
for ff=1:length(frameVec)
    for kk = 1:length(heightVec)
        for jj =1:Ny
            [bin_midVal_ul_binFreq(:,jj,kk,ff),ul_edges,bin_midVal_ul_binIndx(:,jj,kk,ff)] = ...
                histcounts(bin_midVal_ul(:,jj,kk,ff),ul_edges);
        end
    end
end
disp(['u_l  Binning Done !']);

%% calculate <u_s^2+(u_l^+,z)>
u_s_fluc_cond_u_l_sum = zeros([length(midval_ul) length(heightVec)]);
uw_variance_fluc_cond_u_l_sum = zeros([length(midval_ul) length(heightVec)]);
uw_mean_cond_u_l_sum = zeros([length(midval_ul) length(heightVec)]);
u_s_waveLen_cond_u_l_sum = zeros([length(midval_ul) length(heightVec)]);
u_s_zeroCrossings_u_l_sum = zeros([length(midval_ul) length(heightVec)]);
N_u_l = zeros([length(midval_ul) length(heightVec)]);
bin_midVal_ul_binIndx_permuted = permute(bin_midVal_ul_binIndx,[1 2 4 3]);
bin_variance_us_permuted = permute(bin_variance_us, [1 2 4 3]); 
bin_std_uw_permuted = permute(bin_std_uw, [1 2 4 3]); 
bin_mean_uw_permuted = permute(bin_mean_uw, [1 2 4 3]);
bin_waveLen_us_permuted = permute(bin_waveLen_us, [1 2 4 3]);
bin_zeroCrossing_us_permuted = permute(bin_zeroCrossing_us,[1 2 4 3]);
combined_mean_uw = mean(bin_mean_uw(:));
for kk = 1:length(heightVec)
  for ff = 1:length(frameVec)
    for jj =1:Ny
      for ii = 1:length(horzEdges)-1
        indx = bin_midVal_ul_binIndx_permuted(ii,jj,ff,kk);
         if(indx)
            u_s_fluc_cond_u_l_sum(indx,kk) = u_s_fluc_cond_u_l_sum(indx,kk) + bin_variance_us_permuted(ii,jj,ff,kk);
            u_s_waveLen_cond_u_l_sum(indx,kk) = u_s_waveLen_cond_u_l_sum(indx,kk) + bin_waveLen_us_permuted(ii,jj,ff,kk);
            u_s_zeroCrossings_u_l_sum(indx,kk) = u_s_zeroCrossings_u_l_sum(indx,kk) + bin_zeroCrossing_us_permuted(ii,jj,ff,kk);
            uw_variance_fluc_cond_u_l_sum(indx,kk) = uw_variance_fluc_cond_u_l_sum(indx,kk)+ (bin_std_uw_permuted(ii,jj,ff,kk)^2+...
                        (bin_mean_uw_permuted(ii,jj,ff,kk)-combined_mean_uw)^2)*horz_bin_sample_size;
            uw_mean_cond_u_l_sum(indx,kk) = uw_mean_cond_u_l_sum(indx,kk) + bin_mean_uw_permuted(ii,jj,ff,kk)*...
                horz_bin_sample_size;
            N_u_l(indx,kk) = N_u_l(indx,kk)+1;
         end
      end
    end
  end
end
uw_std_cond_u_l = zeros([length(midval_ul) length(heightVec)]);
uw_mean_cond_u_l = zeros([length(midval_ul) length(heightVec)]);
% calculate combined standard deviationof uw
for kk = 1:length(heightVec)
    uw_std_cond_u_l(:,kk) = sqrt(uw_variance_fluc_cond_u_l_sum(:,kk)./(numel(bin_mean_uw_permuted(:,:,:,kk))*...
                horz_bin_sample_size));
    uw_mean_cond_u_l(:,kk) = uw_mean_cond_u_l_sum(:,kk)./(numel(bin_mean_uw_permuted(:,:,:,kk))*...
                horz_bin_sample_size);
end
% calculate pdf of u_l & average u_s^2, zerocrossings, wavelength
pdf_N_u_l = zeros([length(midval_ul) length(heightVec)]);
us_conditional_u_l = zeros([length(midval_ul) length(heightVec)]);
us_zeroCrossing_conditional_u_l = zeros([length(midval_ul) length(heightVec)]);
us_waveLen_conditional_u_l =  zeros([length(midval_ul) length(heightVec)]);
for kk =1:size(pdf_N_u_l,2)
  for ii = 1:size(pdf_N_u_l,1)
    pdf_N_u_l(ii,kk) = N_u_l(ii,kk)/sum(N_u_l(:,kk));
    us_conditional_u_l(ii,kk) = u_s_fluc_cond_u_l_sum(ii,kk)/N_u_l(ii,kk); % this is equation 4.1
    us_zeroCrossing_conditional_u_l(ii,kk) = u_s_zeroCrossings_u_l_sum(ii,kk)/N_u_l(ii,kk); % avg wavelen in u_L bin
    us_waveLen_conditional_u_l(ii,kk) = u_s_waveLen_cond_u_l_sum(ii,kk)/N_u_l(ii,kk); % average no of zero crossings
  end
end

figWidth = 3.0; figHeight = figWidth;  formatEng = '-depsc2';
xvec = (0:Nx-1).*dx./BLH;
% figure(); 
% set(gcf,'Renderer','OpenGL');
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
% set(gcf,'PaperPositionMode','Auto');
% set(gcf,'Renderer','painters');
% hold all;
% plot(xvec, u_s(:,Ny/2,2,1),'-xc','LineWidth',1);
% plot(xvec, u_l(:,Ny/2,2,1),'-xm','LineWidth',2);
% xlabel('x/\delta');
% SPR = sprintf('%s%s','','ek02_u_l_u_s.eps');  
% print(formatEng,SPR);

% figure(); 
% set(gcf,'Renderer','OpenGL');
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
% set(gcf,'PaperPositionMode','Auto');
% set(gcf,'Renderer','painters');
% hold all;
% plot(xvec, uwsig(:,Ny/2,2,1),'-xc','LineWidth',1);
% plot(xvec, u_l(:,Ny/2,2,1),'-xm','LineWidth',2);
% xlabel('x/\delta');
% SPR = sprintf('%s','ek02_u_l_uw.eps');  
% print(formatEng,SPR);
save('freqAmpMod_ek02.mat','us_conditional_u_l','midval_ul','pdf_N_u_l',...
    'uw_std_cond_u_l', 'uw_mean_cond_u_l','us_zeroCrossing_conditional_u_l',...
    'us_waveLen_conditional_u_l','frameVec','u_star','BLH','heightVec','-v7.3');
%%
parentDir = '/media/ahsan/ds00/Ekman_Ugeo2_2048x2048x64_Lz750';
readinputs(parentDir);
figWidth = 3.0; figHeight = figWidth;  formatEng = '-depsc2';
LineStyles = {'-o', '--d','-.',':^','-+','-','-s','-d','-^','-*'}; % 10 types total
LineStyles1 = {'-co', '--md','-.g','-b^','-+r','-*m','-sg','-bd','-c^','-r*'}; % 10 types total
LineStyles2 = {'co', 'md','.g','b^','+r','*m','sg','bd','c^','r*'}; % 10 types total
data=load('freqAmpMod_ek02.mat');
midval_ul=data.midval_ul;
us_conditional_u_l = data.us_conditional_u_l;
uw_std_cond_u_l = data.uw_std_cond_u_l;
uw_mean_cond_u_l = data.uw_mean_cond_u_l;
us_zeroCrossing_conditional_u_l = data.us_zeroCrossing_conditional_u_l;
us_waveLen_conditional_u_l = data.us_waveLen_conditional_u_l;
u_star = data.u_star;
heightVec = data.heightVec;
BLH = data.BLH;
pdf_N_u_l = data.pdf_N_u_l;
% plot average u_s^2 for 7 different u_L bin
% figure();
% % set(gcf,'Renderer','OpenGL');
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
% % set(gcf,'PaperPositionMode','Auto');
% set(gcf,'Renderer','painters');
% hold all;
% ul_val_vec = [57 60 63 66 69 72 75];% 66 is u_l = 0;
% for i = 1:length(ul_val_vec)
%   plot(heightVec(2:end),us_conditional_u_l(ul_val_vec(i),2:end),LineStyles{i});  
% end
% lgndVal = cell([numel(ul_val_vec) 1]);
% for ii = 1:numel(ul_val_vec)
%     lgndVal{ii} = sprintf('u_L = %2.3f',midval_ul(ul_val_vec(ii)));
% end
% lgndh = legendflex(lgndVal, 'nrow', 4,'ncol',2,'box','off','ref', gca, 'anchor', [3 3],...
%     'buffer', [-10 -10],'FontSize',8);
% set(gca, 'Xscale','log');
% xlabel('$z/\delta$','Interpreter','Latex');
% ylabel('$\left < u^{2+}_s\right >$','Interpreter','Latex');
% xlim([0.01 1.5]);
% SPR = sprintf('%s','ek02_u_s2_z.eps');  
% print(formatEng,SPR);

% plot average std(uw) for 7 different u_L bin
figure(); 
% set(gcf,'Renderer','OpenGL');
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
% set(gcf,'PaperPositionMode','Auto');
set(gcf,'Renderer','painters');
hold all;
ul_val_vec = [51 56 61 66 71 76 81];% 66 is u_l = 0;
for i = 1:length(ul_val_vec)
%     plot(heightVec(2:end), uw_mean_cond_u_l(ul_val_vec(i),2:end)./u_star,LineStyles1{i});
    plot(heightVec(2:end), 3.*uw_std_cond_u_l(ul_val_vec(i),2:end)./u_star^2, LineStyles2{i}); 
%     plot(heightVec(2:end),-3.*uw_std_cond_u_l(ul_val_vec(i),2:end)./u_star, LineStyles1{i});
end
lgndVal = cell([numel(ul_val_vec) 1]);
for ii = 1:numel(ul_val_vec)
    lgndVal{ii} = sprintf('u_L = %2.3f',midval_ul(ul_val_vec(ii)));
end
lgndh = legendflex(lgndVal, 'nrow', 4,'ncol',2,'box','off','ref', gca, 'anchor', [3 3],...
    'buffer', [-10 -10],'FontSize',8);
set(gca, 'Xscale','log');
for i = 1:length(ul_val_vec)
    plot(heightVec(2:end), uw_mean_cond_u_l(ul_val_vec(i),2:end)./u_star^2,LineStyles1{i},'LineWidth',0.5,'MarkerSize',4);
%     plot(heightVec(2:end), 3.*uw_std_cond_u_l(ul_val_vec(i),2:end)./u_star, LineStyles1{i}); 
    plot(heightVec(2:end),-3.*uw_std_cond_u_l(ul_val_vec(i),2:end)./u_star^2, LineStyles2{i},'LineWidth',0.5,'MarkerSize',4);
end
xlabel('$z/\delta$','Interpreter','Latex');
ylabel('$\left < {u^{\prime} w^{\prime}}\right >/u_*^2$','Interpreter','Latex');
xlim([0.007 1.1]); ylim([-2.0 4]);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [1 1 1], ...
  'YColor'      , [1 1 1], ...
  'XTick'       , [.001 0.01 0.1 1], ...
  'LineWidth'   , 0.5         );
SPR = sprintf('%s','ek02_u_l_uw_std.eps');  
print(formatEng,SPR);

% equation 4.3, 5.2 <u_s^2+(u_l^+ == 0,z)>
ul_zero_indx = find(~midval_ul);
% plot Delta u_s^2 vs y
tmp_vec = [60 63 69 72]; % 66 is u_l = 0;
figure(); 
% set(gcf,'Renderer','OpenGL');
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
% set(gcf,'PaperPositionMode','Auto');
set(gcf,'Renderer','painters');
hold all;
for i = 1:length(tmp_vec)
  dat= (us_conditional_u_l(tmp_vec (i),:)-us_conditional_u_l(ul_zero_indx,:))./...
      us_conditional_u_l(ul_zero_indx,:).*100; % eqn 4.3
  plot(heightVec,dat,LineStyles{i},'LineWidth',0.5,'MarkerSize',4);  
end
lgndVal = cell([numel(tmp_vec) 1]);
for ii = 1:numel(tmp_vec)
    lgndVal{ii} = sprintf('u_L = %2.3f',midval_ul(tmp_vec(ii)));
end
lgndh = legendflex(lgndVal, 'nrow', 4,'ncol',2,'box','off','ref', gca, 'anchor', [1 1],...
    'buffer', [10 -10],'FontSize',8);
set(gca, 'Xscale','log');
xlabel('$z/\delta$','Interpreter','Latex');
ylabel('$\Delta \left < u^{2+}_s\right >$','Interpreter','Latex');
xlim([0.01 1.1]);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XTick'       , [.001 0.01 0.1 1], ...
  'LineWidth'   , 0.5         );
SPR = sprintf('%s','ek02_u_s2_z_multiple_u_l_0.eps');  
print(formatEng,SPR);

% % plot Delta u_s^2 vs u_L
% tmp_vec = [2 4 8 12];
% figure(); 
% set(gcf,'Renderer','OpenGL');
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
% set(gcf,'PaperPositionMode','Auto');
% set(gcf,'Renderer','painters');
% hold all;
% for i = 1:length(tmp_vec)
%   dat= (us_conditional_u_l(:,tmp_vec (i))-us_conditional_u_l(ul_zero_indx,tmp_vec(i)))./...
%       us_conditional_u_l(ul_zero_indx,tmp_vec (i)).*100; % eqn 4.3
%   plot(midval_ul,dat,LineStyles{i});  
% end
% lgndVal = cell([numel(tmp_vec) 1]);
% for ii = 1:numel(tmp_vec)
%     lgndVal{ii} = sprintf('z=%2.3f\\delta',heightVec(tmp_vec(ii)));
% end
% lgndh = legendflex(lgndVal, 'nrow', 2,'ncol',2,'box','off','ref', gca, 'anchor', [1 1],...
%     'buffer', [10 -10],'FontSize',8);
% xlabel('$u_L$','Interpreter','Latex');
% ylabel('$\Delta \left < u^{2+}_s\right >$','Interpreter','Latex');
% xlim([-0.3 0.3]);
% ylim([-60 120]);
% SPR = sprintf('%s','ek02_u_s2_uL_multiple_z_ul_0.eps');  
% print(formatEng,SPR);

% convert zeroCrossings to wave no
ul_val_vec = [60 66 72];% 66 is u_l = 0;
us_waveNo_conditional_u_l = (2*pi*z_i)./(2*BLH./us_zeroCrossing_conditional_u_l);
% plot waveNo at different heights for 7 different u_L
figure(); 
% set(gcf,'Renderer','OpenGL');
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
% set(gcf,'PaperPositionMode','Auto');
set(gcf,'Renderer','painters');
hold all;
for i = 1:length(ul_val_vec)
  plot(heightVec(2:end),us_waveNo_conditional_u_l(ul_val_vec(i),2:end),LineStyles{i},'LineWidth',0.5,'MarkerSize',4);  
end
lgndVal = cell([numel(ul_val_vec) 1]);
for ii = 1:numel(ul_val_vec)
    lgndVal{ii} = sprintf('u_L = %2.3f',midval_ul(ul_val_vec(ii)));
end
lgndh = legendflex(lgndVal, 'nrow', 4,'ncol',2,'box','off','ref', gca, 'anchor', [1 1],...
    'buffer', [10 -10],'FontSize',8);
set(gca, 'Xscale','log');
xlabel('$z/\delta$','Interpreter','Latex');
ylabel('$\left < k_{x}|u_L\right >$','Interpreter','Latex');
xlim([0.01 1.1]);
% ylim([-60 120]);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XTick'       , [.001 0.01 0.1 1], ...
  'LineWidth'   , 0.5         );
SPR = sprintf('%s','ek02_kx_ul_z.eps');  
print(formatEng,SPR);

% % plot waveNo at different heights for 7 different u_L from wavelength determination
% figure();
% set(gcf,'Renderer','OpenGL');
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
% set(gcf,'PaperPositionMode','Auto');
% set(gcf,'Renderer','painters');
% hold all;
% for i = 1:length(ul_val_vec)
%   plot(heightVec(2:end),(2*pi*z_i)./us_waveLen_conditional_u_l(ul_val_vec(i),2:end),LineStyles{i});  
% end
% lgndVal = cell([numel(ul_val_vec) 1]);
% for ii = 1:numel(ul_val_vec)
%     lgndVal{ii} = sprintf('u_L = %2.3f',midval_ul(ul_val_vec(ii)));
% end
% lgndh = legendflex(lgndVal, 'nrow', 4,'ncol',2,'box','off','ref', gca, 'anchor', [3 3],...
%     'buffer', [-10 -10],'FontSize',8);
% set(gca, 'Xscale','log');
% xlabel('$z/\delta$','Interpreter','Latex');
% ylabel('$\left < k_{x}|u_L\right >$','Interpreter','Latex');
% xlim([0.01 1.5]);

% % plot pdf of u_l
% figure();
% set(gcf,'Renderer','OpenGL');
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
% set(gcf,'PaperPositionMode','Auto');
% % set(gcf,'Renderer','painters');
% hold all;
% [XX, YY] = meshgrid(midval_ul, heightVec);
% pcolor(XX,YY, pdf_N_u_l');
% xlabel('u_{L}');
% ylabel('z/\delta');
% colormap parula;
% shading interp;
% hcbar = colorbar('northoutside');
% SPR = sprintf('%s','ek02_pcolor_pdf_uL_against_uL_z.eps');  
% print(formatEng,SPR);
% 
% % plot of conditional u^2_s
% figure();
% % set(gcf,'Renderer','OpenGL');
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
% set(gcf,'PaperPositionMode','Auto');
% set(gcf,'Renderer','painters');
% pcolor(XX,YY, us_conditional_u_l');
% xlabel('u_{L}');
% ylabel('z/\delta');
% colormap parula;
% shading interp;
% hcbar = colorbar('northoutside');
% ylabel(hcbar, 'u_s^{2+}');
% SPR = sprintf('%s','ek02_pcolor_us2_against_uL_z.eps');  
% print(formatEng,SPR);
