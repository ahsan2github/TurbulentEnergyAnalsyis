parentDir = '/raid/home/mohammad/Ekman_Ugeo2_2048x2048x64_Lz750';
outputDir = [parentDir '/output/'];
readinputs(parentDir);
u_star_in = u_star; clear u_star;
BLH = 550.0;
ustar = 0.09;
refLevels = 3:3:62; %64; % total 21 levels
nscale = log(2048)/log(2);
data = load('wav_t_m_n_ek02_1.mat');
tenergy = data.tenergy;
LineMarkers = {'co-', 'm--', 'b.-', 'r:', 'c-+', 'm-o', 'b-s', 'r-d', 'c-^','m-*'}; % 10 types total
LineColors = {rgb('Salmon'), rgb('Pink'), rgb('Orange'), rgb('Cyan'), rgb('SkyBlue'),...
              rgb('Violet'), rgb('Silver'), rgb('Tomato'), rgb('Olive'), rgb('Blue')};  %10 types total         
dummy = mClass(nscale);
nframe = 4;
%%
%
figWidth = 6.5; figHeight = figWidth/1.2;  formatEng = '-depsc';
opengl('save','software');
m_st = 5; 
%plot a fixed n vs all m at different heightsi
figh = figure(); set(gcf,'Renderer','painters');
set(figh, 'Visible','on');
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
set(gcf, 'PaperPositionMode','Auto');
text_x_pos = 0.55;
text_y_pos = 0.85;
%----------------------------------------------------------------------------
subplot(5,3,[1,4]);
select_n = 10;
possible_m = m_st:select_n;
m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(refLevels) length(possible_m)]);
for ff = 1:nframe 
    for hh = 1:length(refLevels)
        for mmm = 1:size(m_data_frame,2)
            m = possible_m(mmm); n_pos = select_n-m+1;
            m_data_frame(hh, mmm, ff) = data.frame_energy{ff}{hh}.n_data(m).tmn(n_pos);
            sigma_m_data_frame(hh,mmm,ff)= data.frame_energy{ff}{hh}.n_data(m).sigma_tmn(n_pos);
            m_in_m(hh,mmm) = tenergy{hh}.m_in_m(m); 
        end
    end
end
m_data = mean(m_data_frame,3);
sigma_m_data = mean(sigma_m_data_frame,3);
heightVec = zeros([1 length(refLevels)]);
for hh=1:length(refLevels)
    heightVec(hh) = (data.tenergy{hh}.h);
end
hold all;
for mmm = 1:size(m_data,2)
    plot(heightVec./BLH, m_data(:,mmm)./ustar^3,LineMarkers{mmm},'MarkerSize',3);
end
plot([0 1], [0 0],'--k');
ylim([-20 10]);
lgndVal = cell([numel(possible_m) 1]);
for ll = 1:numel(possible_m)
    lgndVal{ll} = sprintf('r_m/\\delta=%2.2f',tenergy{1}.m_in_m(possible_m(ll))/BLH);
end
ann_str = sprintf('r_n\\equiv%2.2f\\delta',n_scale/BLH);
text(text_x_pos,text_y_pos*0.17, ann_str,'Units','Normalized','FontSize',8);
set(gca,'FontSize',8);
ylabel('T^{(m,n)}(z,k_m|k_n)/u_*^3','FontSize',8);
set(gca,'XMinorTick','on','YMinorTick','on');
box on;
%----------------------------------------------------------------------------
subplot(5,3,[2,5]);
select_n = select_n-1;
possible_m = m_st:select_n;
m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(refLevels) length(possible_m)]);
for ff = 1:nframe
    for hh = 1:length(refLevels)
        for mmm = 1:size(m_data_frame,2)
            m = possible_m(mmm); n_pos = select_n-m+1;
            m_data_frame(hh, mmm, ff) = data.frame_energy{ff}{hh}.n_data(m).tmn(n_pos);
            sigma_m_data_frame(hh,mmm,ff)= data.frame_energy{ff}{hh}.n_data(m).sigma_tmn(n_pos);
            m_in_m(hh,mmm) = tenergy{hh}.m_in_m(m); 
        end
    end
end
m_data = mean(m_data_frame,3);
sigma_m_data = mean(sigma_m_data_frame,3);
heightVec = zeros([1 length(refLevels)]);
for hh=1:length(refLevels)
    heightVec(hh) = (data.tenergy{hh}.h);
end
hold all;
% bounds = zeros([length(heightVec) 2 length(possible_m)]);
% for mmm = 1:size(m_data,2)
%     bounds(:,:,mmm) = repmat(ones([length(heightVec) 1]).*sigma_m_data(:,mmm),[1 2]);
% end
% [l,p]= boundedline(heightVec./BLH, m_data(:,:)./ustar^3,...
%         bounds./ustar^3.*0.5,'cmap', parula(10), 'transparency', 0.3);
for mmm = 1:size(m_data,2)
    plot(heightVec./BLH, m_data(:,mmm)./ustar^3,LineMarkers{mmm});
end
plot([0 1], [0 0],'--k');
ylim([-20 10]);
ann_str = sprintf('r_n\\equiv%2.2f\\delta',n_scale/BLH);
text(text_x_pos,text_y_pos*0.17, ann_str,'Units','Normalized','FontSize',8);
set(gca,'FontSize',8);
set(gca,'XMinorTick','on','YMinorTick','on');
box on;
%----------------------------------------------------------------------------
subplot(5,3,[3,6]);
select_n = select_n-1;
possible_m = m_st:select_n;
m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(refLevels) length(possible_m)]);
for ff = 1:nframe
    for hh = 1:length(refLevels)
        for mmm = 1:size(m_data_frame,2)
            m = possible_m(mmm); n_pos = select_n-m+1;
            m_data_frame(hh, mmm, ff) = data.frame_energy{ff}{hh}.n_data(m).tmn(n_pos);
            sigma_m_data_frame(hh,mmm,ff)= data.frame_energy{ff}{hh}.n_data(m).sigma_tmn(n_pos);
            m_in_m(hh,mmm) = tenergy{hh}.m_in_m(m); 
        end
    end
end
m_data = mean(m_data_frame,3);
sigma_m_data = mean(sigma_m_data_frame,3);
heightVec = zeros([1 length(refLevels)]);
for hh=1:length(refLevels)
    heightVec(hh) = (data.tenergy{hh}.h);
end
hold all;
% bounds = zeros([length(heightVec) 2 length(possible_m)]);
% for mmm = 1:size(m_data,2)
%     bounds(:,:,mmm) = repmat(ones([length(heightVec) 1]).*sigma_m_data(:,mmm),[1 2]);
% end
% [l,p]= boundedline(heightVec./BLH, m_data(:,:)./ustar^3,...
% bounds./ustar^3.*0.5,'cmap', parula(10), 'transparency', 0.3);
for mmm = 1:size(m_data,2)
    plot(heightVec./BLH, m_data(:,mmm)./ustar^3,LineMarkers{mmm});
end
plot([0 1], [0 0],'--k');
ylim([-20 10]);
ann_str = sprintf('r_n\\equiv%2.2f\\delta',n_scale/BLH);
text(text_x_pos,text_y_pos*0.17, ann_str,'Units','Normalized','FontSize',8);
xlabel('z/\delta','FontSize',8);
set(gca,'FontSize',8);
set(gca,'XMinorTick','on','YMinorTick','on');
box on;
%----------------------------------------------------------------------------
subplot(5,3,[7,10]);
select_n = select_n-1;
possible_m = m_st:select_n;
m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(refLevels) length(possible_m)]);
for ff = 1:nframe
    for hh = 1:length(refLevels)
        for mmm = 1:size(m_data_frame,2)
            m = possible_m(mmm); n_pos = select_n-m+1;
            m_data_frame(hh, mmm, ff) = data.frame_energy{ff}{hh}.n_data(m).tmn(n_pos);
            sigma_m_data_frame(hh,mmm,ff)= data.frame_energy{ff}{hh}.n_data(m).sigma_tmn(n_pos);
            m_in_m(hh,mmm) = tenergy{hh}.m_in_m(m); 
        end
    end
end
m_data = mean(m_data_frame,3);
sigma_m_data = mean(sigma_m_data_frame,3);
heightVec = zeros([1 length(refLevels)]);
for hh=1:length(refLevels)
    heightVec(hh) = (data.tenergy{hh}.h);
end
hold all;
% bounds = zeros([length(heightVec) 2 length(possible_m)]);
% for mmm = 1:size(m_data,2)
%     bounds(:,:,mmm) = repmat(ones([length(heightVec) 1]).*sigma_m_data(:,mmm),[1 2]);
% end
% [l,p]= boundedline(heightVec./BLH, m_data(:,:)./ustar^3,...
%         bounds./ustar^3.*0.5,'cmap', parula(10), 'transparency', 0.3);
for mmm = 1:size(m_data,2)
    plot(heightVec./BLH, m_data(:,mmm)./ustar^3,LineMarkers{mmm});
end
plot([0 1], [0 0],'--k');
ylim([-20 10]);
ann_str = sprintf('r_n\\equiv%2.2f\\delta',n_scale/BLH);
text(text_x_pos,text_y_pos*.17, ann_str,'Units','Normalized','FontSize',8);
xlabel('z/\delta','FontSize',8);
ylabel('T^{(m,n)}(z,k_m|k_n)/u_*^3','FontSize',8);
set(gca,'FontSize',8);
%ylim([-0.6 0.6]);
set(gca,'XMinorTick','on','YMinorTick','on');
box on;
%----------------------------------------------------------------------------
subplot(5,3,[8,11]);
select_n = select_n-1;
possible_m = m_st:select_n;
m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(refLevels) length(possible_m)]);
for ff = 1:nframe
    for hh = 1:length(refLevels)
        for mmm = 1:size(m_data_frame,2)
            m = possible_m(mmm); n_pos = select_n-m+1;
            m_data_frame(hh, mmm, ff) = data.frame_energy{ff}{hh}.n_data(m).tmn(n_pos);
            sigma_m_data_frame(hh,mmm,ff)= data.frame_energy{ff}{hh}.n_data(m).sigma_tmn(n_pos);
            m_in_m(hh,mmm) = tenergy{hh}.m_in_m(m); 
        end
    end
end
m_data = mean(m_data_frame,3);
sigma_m_data = mean(sigma_m_data_frame,3);
heightVec = zeros([1 length(refLevels)]);
for hh=1:length(refLevels)
    heightVec(hh) = (data.tenergy{hh}.h);
end
hold all;
% bounds = zeros([length(heightVec) 2 length(possible_m)]);
% for mmm = 1:size(m_data,2)
%     bounds(:,:,mmm) = repmat(ones([length(heightVec) 1]).*sigma_m_data(:,mmm),[1 2]);
% end
% [l,p]= boundedline(heightVec./BLH, m_data(:,:)./ustar^3,...
%         bounds./ustar^3.*0.5,'cmap', parula(10), 'transparency', 0.3);
for mmm = 1:size(m_data,2)
    plot(heightVec./BLH, m_data(:,mmm)./ustar^3,LineMarkers{mmm});
end    
plot([0 1], [0 0],'--k');
ylim([-20 10]);
ann_str = sprintf('r_n\\equiv%2.2f\\delta',n_scale/BLH);
text(text_x_pos,text_y_pos*0.17, ann_str,'Units','Normalized','FontSize',8);
xlabel('z/\delta','FontSize',8);
set(gca,'FontSize',8);
%ylim([-0.6 0.6]);
set(gca,'XMinorTick','on','YMinorTick','on');
box on;
%----------------------------------------------------------------------------
subplot(5,3,[9,12]);
select_n = select_n-1;
possible_m = m_st:select_n;
m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(refLevels) length(possible_m)]);
for ff = 1:nframe
    for hh = 1:length(refLevels)
        for mmm = 1:size(m_data_frame,2)
            m = possible_m(mmm); n_pos = select_n-m+1;
            m_data_frame(hh, mmm, ff) = data.frame_energy{ff}{hh}.n_data(m).tmn(n_pos);
            sigma_m_data_frame(hh,mmm,ff)= data.frame_energy{ff}{hh}.n_data(m).sigma_tmn(n_pos);
            m_in_m(hh,mmm) = tenergy{hh}.m_in_m(m); 
        end
    end
end
m_data = mean(m_data_frame,3);
sigma_m_data = mean(sigma_m_data_frame,3);
heightVec = zeros([1 length(refLevels)]);
for hh=1:length(refLevels)
    heightVec(hh) = (data.tenergy{hh}.h);
end
hold all;
% bounds = zeros([length(heightVec) 2 length(possible_m)]);
% for mmm = 1:size(m_data,2)
%     bounds(:,:,mmm) = repmat(ones([length(heightVec) 1]).*sigma_m_data(:,mmm),[1 2]);
% end
% [l,p]= boundedline(heightVec./BLH, m_data(:,:)./ustar^3,...
%         bounds./ustar^3.*0.5,'cmap', parula(10), 'transparency', 0.3);
for mmm = 1:size(m_data,2)
    plot(heightVec./BLH, m_data(:,mmm)./ustar^3,LineMarkers{mmm});
end
plot([0 1], [0 0],'--k');
ylim([-20 10]);
ann_str = sprintf('r_n\\equiv%2.2f\\delta',n_scale/BLH);
text(text_x_pos,text_y_pos*0.17, ann_str,'Units','Normalized','FontSize',8);
xlabel('z/\delta','FontSize',8);
set(gca,'FontSize',8);
%ylim([-0.6 0.6]);
set(gca,'XMinorTick','on','YMinorTick','on');
box on;
%----------------------------------------------------------------------------
hax = subplot(5,3,[13,14,15], 'Position',[0.15 0.07 0.72 0.12]);
box on;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
hold all;
plot([0,1],[2.2,2.2],LineMarkers{1}); text(1.2,2.3,lgndVal{1});
plot([3,4],[2.2,2.2],LineMarkers{2}); text(4.2,2.3,lgndVal{2});
plot([6,7],[2.2,2.2],LineMarkers{3}); text(7.2,2.3,lgndVal{3});
plot([0,1],[1.3,1.35],LineMarkers{4}); text(1.2,1.35,lgndVal{4});
plot([3,4],[1.3,1.35],LineMarkers{5}); text(4.2,1.35,lgndVal{5});
plot([6,7],[1.3,1.35],LineMarkers{6}); text(7.2,1.35,lgndVal{6});

xlim([-0.5 9]);
%ylim([0.5 3]);

SPR = sprintf('%s%','tmn_ek02_fixed_n-m_n_equal.eps');  
set(gcf, 'Color', 'w');
print(formatEng,SPR);


%% plot m vs smallest n at different heights with one standard deviation
%===================================================================================================
text_x_pos = 0.15;
text_y_pos = 0.10;
fig = figure(); set(gcf,'Renderer','painters');
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6 3.0]);
set(gcf,'PaperPositionMode','auto');
m_st = 1;
% == how all scales larger than 1.33\delta interact with smaller scales ===
select_n = 6; select_h = 1; 
possible_m = m_st:select_n-1;
m_data_frame = zeros([length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for ff = 1:nframe
    for mmm = 1:size(m_data_frame,1)
        m = possible_m(mmm); n_pos = select_n-m+1;
        m_data_frame(mmm,ff) = data.frame_energy{ff}{select_h}.n_data(m).tmn(n_pos);
        sigma_m_data_frame(mmm,ff)= data.frame_energy{ff}{select_h}.n_data(m).sigma_tmn(n_pos);
        m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
    end
end
m_data = mean(m_data_frame,2);
sigma_m_data = mean(sigma_m_data_frame,2);

subplot(241);
hold all;
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3,'-or', 'transparency', 0.3);
    
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH, m_data./ustar^3 +sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH, m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.5,text_y_pos, ann_str,'Units','Normalized');
ylabel('T^{(m,n)}(z,k_m|k_n)/u_*^3','FontSize',8);
set(gca,'XMinorTick','on','YMinorTick','on');

subplot(242);
hold all;
select_n = 6; select_h = 2;
possible_m = m_st:select_n-1;
m_data_frame = zeros([length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for ff = 1:nframe
    for mmm = 1:size(m_data_frame,1)
        m = possible_m(mmm); n_pos = select_n-m+1;
        m_data_frame(mmm,ff) = data.frame_energy{ff}{select_h}.n_data(m).tmn(n_pos);
        sigma_m_data_frame(mmm,ff)= data.frame_energy{ff}{select_h}.n_data(m).sigma_tmn(n_pos);
        m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
    end
end
m_data = mean(m_data_frame,2);
sigma_m_data = mean(sigma_m_data_frame,2);

bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH, m_data./ustar^3 + sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH, m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.5,text_y_pos, ann_str,'Units','Normalized');
set(gca,'XMinorTick','on','YMinorTick','on');

subplot(243); hold all;
select_n = 6; select_h = 6;
possible_m = m_st:select_n-1;
m_data_frame = zeros([length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for ff = 1:nframe
    for mmm = 1:size(m_data_frame,1)
        m = possible_m(mmm); n_pos = select_n-m+1;
        m_data_frame(mmm,ff) = data.frame_energy{ff}{select_h}.n_data(m).tmn(n_pos);
        sigma_m_data_frame(mmm,ff)= data.frame_energy{ff}{select_h}.n_data(m).sigma_tmn(n_pos);
        m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
    end
end
m_data = mean(m_data_frame,2);
sigma_m_data = mean(sigma_m_data_frame,2);
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH, m_data./ustar^3 + sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH, m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.5,text_y_pos, ann_str,'Units','Normalized');
% xlabel('2^{m}h/\delta');
set(gca,'XMinorTick','on','YMinorTick','on');

subplot(244);
hold all;
select_n = 6; select_h = 15;
possible_m = m_st:select_n-1;
m_data_frame = zeros([length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for ff = 1:nframe
    for mmm = 1:size(m_data_frame,1)
        m = possible_m(mmm); n_pos = select_n-m+1;
        m_data_frame(mmm,ff) = data.frame_energy{ff}{select_h}.n_data(m).tmn(n_pos);
        sigma_m_data_frame(mmm,ff)= data.frame_energy{ff}{select_h}.n_data(m).sigma_tmn(n_pos);
        m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
    end
end
m_data = mean(m_data_frame,2);
sigma_m_data = mean(sigma_m_data_frame,2);
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH, m_data./ustar^3 + sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH, m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.5,text_y_pos, ann_str,'Units','Normalized');
% xlabel('2^{m}h/\delta');
set(gca,'XMinorTick','on','YMinorTick','on');

% == how all scales larger than 5.33\delta interact with smaller scales ===
subplot(245);
hold all;
select_n = 4; select_h = 1;
possible_m = m_st:select_n-1;
m_data_frame = zeros([length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for ff = 1:nframe
    for mmm = 1:size(m_data_frame,1)
        m = possible_m(mmm); n_pos = select_n-m+1;
        m_data_frame(mmm,ff) = data.frame_energy{ff}{select_h}.n_data(m).tmn(n_pos);
        sigma_m_data_frame(mmm,ff)= data.frame_energy{ff}{select_h}.n_data(m).sigma_tmn(n_pos);
        m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
    end
end
m_data = mean(m_data_frame,2);
sigma_m_data = mean(sigma_m_data_frame,2);
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH, m_data./ustar^3 + sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH, m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.5,text_y_pos, ann_str,'Units','Normalized');
xlabel('r_{m}/\delta');
ylabel('T^{(m,n)}(z,k_m|k_n)/u_*^3','FontSize',8);
set(gca,'XMinorTick','on','YMinorTick','on');

subplot(246);
hold all;
select_n = 4; select_h = 2;
m_data_feame = zeros([length(possible_m) nframe]);
m_data_frame = zeros([length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for ff = 1:nframe
    for mmm = 1:size(m_data_frame,1)
        m = possible_m(mmm); n_pos = select_n-m+1;
        m_data_frame(mmm,ff) = data.frame_energy{ff}{select_h}.n_data(m).tmn(n_pos);
        sigma_m_data_frame(mmm,ff)= data.frame_energy{ff}{select_h}.n_data(m).sigma_tmn(n_pos);
        m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
    end
end
m_data = mean(m_data_frame,2);
sigma_m_data = mean(sigma_m_data_frame,2);
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH, m_data./ustar^3 +sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH, m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.5,text_y_pos, ann_str,'Units','Normalized');
xlabel('r_{m}/\delta');
set(gca,'XMinorTick','on','YMinorTick','on');

subplot(247);
hold all;
select_n = 4; select_h = 6;
possible_m = m_st:select_n-1;
m_data_frame = zeros([length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for ff = 1:nframe
    for mmm = 1:size(m_data_frame,1)
        m = possible_m(mmm); n_pos = select_n-m+1;
        m_data_frame(mmm,ff) = data.frame_energy{ff}{select_h}.n_data(m).tmn(n_pos);
        sigma_m_data_frame(mmm,ff)= data.frame_energy{ff}{select_h}.n_data(m).sigma_tmn(n_pos);
        m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
    end
end
m_data = mean(m_data_frame,2);
sigma_m_data = mean(sigma_m_data_frame,2);
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH, m_data./ustar^3 +sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH, m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.5,text_y_pos, ann_str,'Units','Normalized');
xlabel('r_{m}/\delta');
set(gca,'XMinorTick','on','YMinorTick','on');

subplot(248);
hold all;
select_n = 4; select_h = 15;
possible_m = m_st:select_n-1;
m_data_frame = zeros([length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for ff = 1:nframe
    for mmm = 1:size(m_data_frame,1)
        m = possible_m(mmm); n_pos = select_n-m+1;
        m_data_frame(mmm,ff) = data.frame_energy{ff}{select_h}.n_data(m).tmn(n_pos);
        sigma_m_data_frame(mmm,ff)= data.frame_energy{ff}{select_h}.n_data(m).sigma_tmn(n_pos);
        m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
    end
end
m_data = mean(m_data_frame,2);
sigma_m_data = mean(sigma_m_data_frame,2);
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH, m_data./ustar^3 +sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH, m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.5,text_y_pos, ann_str,'Units','Normalized');
xlabel('r_{m}/\delta');
set(gca,'XMinorTick','on','YMinorTick','on');

SPR = sprintf('%s','tmn_chnl_fixed_n_fixed_h.eps');  
set(gcf, 'Color', 'w');
print(formatEng,SPR);
