parentDir = '/raid/home/mohammad/Ekman_Ugeo2_2048x2048x64_Lz750';
outputDir = [parentDir '/output/'];
readinputs(parentDir);
u_star_in = u_star; clear u_star;
BLH = 550.69;
ustar = 0.09;
refLevels = 3:3:62; %64; % total 21 levels
nscale = log(2048)/log(2);
data = load('wav_t_m_n_ek02.mat');
tenergy = data.tenergy;
LineMarkers = {'co-', 'm--', 'b.-', 'r:', 'c-+', 'm-o', 'b-s', 'r-d', 'c-^','m-*'}; % 10 types total
LineColors = {rgb('Salmon'), rgb('Pink'), rgb('Orange'), rgb('Cyan'), rgb('SkyBlue'),...
              rgb('Violet'), rgb('Silver'), rgb('Tomato'), rgb('Olive'), rgb('Blue')};  %10 types total         
dummy = mClass(nscale);
nframe = 3;
%%
figWidth = 10; figHeight = figWidth/1.6;  formatEng = '-depsc';
opengl('save','software');
m_st = 4; 
%plot a fixed n vs all m at different heights
figh = figure(); set(gcf,'Renderer','painters');
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figHeight]);
set(gcf,'PaperPositionMode','auto');
text_x_pos = 0.60;
text_y_pos = 0.85;

subplot(5,3,[1,4]);
select_n = 10;
possible_m = m_st:select_n-1;
m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(refLevels) length(possible_m)]);
for ff = 1:nframe 
    for hh = 1:length(refLevels)
        for mmm = 1:size(m_data_frame,2)
            m = possible_m(mmm); n_pos = select_n-m;
            m_data_frame(hh, mmm, ff) = -data.frame_energy{ff}{hh}.n_data(m).tmn(n_pos);
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
    plot(heightVec./BLH, m_data(:,mmm)./ustar^3,LineMarkers{mmm});
end
plot([0 1], [0 0],'--k');

lgndVal = cell([numel(possible_m) 1]);
for ll = 1:numel(possible_m)
    lgndVal{ll} = sprintf('r_n/\\delta=%2.2f',tenergy{1}.m_in_m(possible_m(ll))/BLH);
end
% lgndh = legendflex(lgndVal, 'nrow', 1,'ncol',6,'box','off','ref', ...
%     gca, 'anchor', [6 2], 'buffer', [0 -10],'FontSize',8);
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
text(text_x_pos,text_y_pos*0.9, ann_str,'Units','Normalized','FontSize',8);
set(gca,'FontSize',8);
ylabel('t^{(m,n)}(z,k_m|k_n)/u_*^3','FontSize',8);
xlim([0 1.45]);
box on;
%----------------------------------------------------------------------------
subplot(5,3,[2,5]);
select_n = 9;
possible_m = m_st:select_n-1;
m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(refLevels) length(possible_m)]);
for ff = 1:nframe
    for hh = 1:length(refLevels)
        for mmm = 1:size(m_data_frame,2)
            m = possible_m(mmm); n_pos = select_n-m;
            m_data_frame(hh, mmm, ff) = -data.frame_energy{ff}{hh}.n_data(m).tmn(n_pos);
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
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
text(text_x_pos,text_y_pos*0.9, ann_str,'Units','Normalized','FontSize',8);
set(gca,'FontSize',8);
xlim([0 1.45]);
box on;
%----------------------------------------------------------------------------
subplot(5,3,[3,6]);
select_n = 8;
possible_m = m_st:select_n-1;
m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(refLevels) length(possible_m)]);
for ff = 1:nframe
    for hh = 1:length(refLevels)
        for mmm = 1:size(m_data_frame,2)
            m = possible_m(mmm); n_pos = select_n-m;
            m_data_frame(hh, mmm, ff) = -data.frame_energy{ff}{hh}.n_data(m).tmn(n_pos);
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
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
text(text_x_pos,text_y_pos*0.9, ann_str,'Units','Normalized','FontSize',8);
xlabel('z/\delta','FontSize',8);
set(gca,'FontSize',8);
xlim([0 1.45]);
box on;
%----------------------------------------------------------------------------
subplot(5,3,[7,10]);
select_n = 7;
possible_m = m_st:select_n-1;
m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(refLevels) length(possible_m)]);
for ff = 1:nframe
    for hh = 1:length(refLevels)
        for mmm = 1:size(m_data_frame,2)
            m = possible_m(mmm); n_pos = select_n-m;
            m_data_frame(hh, mmm, ff) = -data.frame_energy{ff}{hh}.n_data(m).tmn(n_pos);
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
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
text(text_x_pos,text_y_pos*1, ann_str,'Units','Normalized','FontSize',8);
xlabel('z/\delta','FontSize',8);
ylabel('t^{(m,n)}(z,k_m|k_n)/u_*^3','FontSize',8);
set(gca,'FontSize',8);
ylim([-0.6 0.6]);
xlim([0 1.45]);
box on;
%----------------------------------------------------------------------------
subplot(5,3,[8,11]);
select_n = 6;
possible_m = m_st:select_n-1;
m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(refLevels) length(possible_m)]);
for ff = 1:nframe
    for hh = 1:length(refLevels)
        for mmm = 1:size(m_data_frame,2)
            m = possible_m(mmm); n_pos = select_n-m;
            m_data_frame(hh, mmm, ff) = -data.frame_energy{ff}{hh}.n_data(m).tmn(n_pos);
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
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
text(text_x_pos,text_y_pos*0.9, ann_str,'Units','Normalized','FontSize',8);
xlabel('z/\delta','FontSize',8);
set(gca,'FontSize',8);
ylim([-0.6 0.6]);
xlim([0 1.45]);
box on;
%----------------------------------------------------------------------------
subplot(5,3,[9,12]);
select_n = 5;
possible_m = m_st:select_n-1;
m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
sigma_m_data_frame = zeros([length(refLevels) length(possible_m) nframe]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(refLevels) length(possible_m)]);
for ff = 1:nframe
    for hh = 1:length(refLevels)
        for mmm = 1:size(m_data_frame,2)
            m = possible_m(mmm); n_pos = select_n-m;
            m_data_frame(hh, mmm, ff) = -data.frame_energy{ff}{hh}.n_data(m).tmn(n_pos);
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
ann_str = sprintf('n\\equiv%2.2f\\delta',n_scale/BLH);
text(text_x_pos,text_y_pos*0.9, ann_str,'Units','Normalized','FontSize',8);
xlabel('z/\delta','FontSize',8);
%legendflex(lgndVal, 'nrow', 2,'ncol',3,'box','off','ref', gca, 'anchor', [2 6], 'buffer', [-2 -2]);
set(gca,'FontSize',8);
ylim([-0.6 0.6]);
xlim([0 1.45]);
box on;
%----------------------------------------------------------------------------
hax = subplot(5,3,[13,14,15], 'Position',[0.15 0.08 0.72 0.11]);
box on;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
hold all;
plot([0,1],[2.2,2.2],LineMarkers{1}); text(1.2,2.3,lgndVal{1});
plot([3,4],[2.2,2.2],LineMarkers{2}); text(4.2,2.3,lgndVal{2});
plot([6,7],[2.2,2.2],LineMarkers{3}); text(7.2,2.3,lgndVal{3});
plot([0,1],[1,1],LineMarkers{4}); text(1.2,1,lgndVal{4});
plot([3,4],[1,1],LineMarkers{5}); text(4.2,1,lgndVal{5});
plot([6,7],[1,1],LineMarkers{6}); text(7.2,1,lgndVal{6});

xlim([-0.5 9]);
ylim([0.5 3]);

SPR = sprintf('%s%s','','tmn_ek02_fixed_n.eps');  
set(gcf, 'Color', 'w');
print(formatEng,SPR);


%% plot m vs smallest n at different heights with one standard deviation
%===================================================================================================
text_x_pos = 0.15;
text_y_pos = 0.10;
fig = figure(); set(gcf,'Renderer','OpenGL');
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6 3.0]);
set(gcf,'PaperPositionMode','auto');
m_st = 1;
% == how all scales larger than 1.82\delta interact with smaller scales ===
select_n = 7; select_h = 1;
possible_m = m_st:select_n-1;
m_data = zeros([length(possible_m) 1]);
sigma_m_data = zeros([length(possible_m) 1]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for mmm = 1:size(m_data,1)
    m = possible_m(mmm); n_pos = select_n-m;
    m_data(mmm) = -tenergy{select_h}.n_data(m).tmn(n_pos);
    sigma_m_data(mmm)= tenergy{select_h}.n_data(m).sigma_tmn(n_pos);
    m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
end
subplot(241);
hold all;
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3,'-or', 'transparency', 0.3);
    
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH,  m_data./ustar^3 + sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH,  m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
xlim([0 120]);
ann_str = sprintf('n=%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.6,text_y_pos, ann_str,'Units','Normalized');
ylabel('T^{(m,n)}(z,k_m|k_n)/u_*^3','FontSize',8);

subplot(242);
select_n = 7; select_h = 2;
possible_m = m_st:select_n-1;
m_data = zeros([length(possible_m) 1]);
sigma_m_data = zeros([length(possible_m) 1]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for mmm = 1:size(m_data,1)
    m = possible_m(mmm); n_pos = select_n-m;
    m_data(mmm) = -tenergy{select_h}.n_data(m).tmn(n_pos);
    sigma_m_data(mmm)= tenergy{select_h}.n_data(m).sigma_tmn(n_pos);
    m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
end
hold all;
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH,  m_data./ustar^3 + sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH,  m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
xlim([0 120]);
ann_str = sprintf('n=%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.6,text_y_pos, ann_str,'Units','Normalized');

subplot(243);
select_n = 7; select_h = 5;
possible_m = m_st:select_n-1;
m_data = zeros([length(possible_m) 1]);
sigma_m_data = zeros([length(possible_m) 1]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for mmm = 1:size(m_data,1)
    m = possible_m(mmm); n_pos = select_n-m;
    m_data(mmm) = -tenergy{select_h}.n_data(m).tmn(n_pos);
    sigma_m_data(mmm)= tenergy{select_h}.n_data(m).sigma_tmn(n_pos);
    m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
end
hold all;
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH,  m_data./ustar^3,'mo-');  
plot(m_in_m./BLH,  m_data./ustar^3 + sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH,  m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
xlim([0 120]);
ann_str = sprintf('n=%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.6,text_y_pos, ann_str,'Units','Normalized');
% xlabel('2^{m}h/\delta');

subplot(244);
select_n = 7; select_h = 12;
possible_m = m_st:select_n-1;
m_data = zeros([length(possible_m) 1]);
sigma_m_data = zeros([length(possible_m) 1]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for mmm = 1:size(m_data,1)
    m = possible_m(mmm); n_pos = select_n-m;
    m_data(mmm) = -tenergy{select_h}.n_data(m).tmn(n_pos);
    sigma_m_data(mmm)= tenergy{select_h}.n_data(m).sigma_tmn(n_pos);
    m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
end
hold all;
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH,  m_data./ustar^3 + sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH,  m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
ann_str = sprintf('n=%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.6,text_y_pos, ann_str,'Units','Normalized');
xlim([0 120]);

% == how all scales larger than 7.26\delta interact with smaller scales ===
subplot(245);
select_n = 5; select_h = 1;
possible_m = m_st:select_n-1;
m_data = zeros([length(possible_m) 1]);
sigma_m_data = zeros([length(possible_m) 1]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for mmm = 1:size(m_data,1)
    m = possible_m(mmm); n_pos = select_n-m;
    m_data(mmm) = -tenergy{select_h}.n_data(m).tmn(n_pos);
    sigma_m_data(mmm)= tenergy{select_h}.n_data(m).sigma_tmn(n_pos);
    m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
end
hold all;
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH,  m_data./ustar^3 + sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH,  m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
xlim([0 120]);
ann_str = sprintf('n=%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.6,text_y_pos, ann_str,'Units','Normalized');
xlabel('r_{m}/\delta');
ylabel('T^{(m,n)}(z,k_m|k_n)/u_*^3','FontSize',8);

subplot(246);
select_n = 5; select_h = 2;
possible_m = m_st:select_n-1;
m_data = zeros([length(possible_m) 1]);
sigma_m_data = zeros([length(possible_m) 1]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for mmm = 1:size(m_data,1)
    m = possible_m(mmm); n_pos = select_n-m;
    m_data(mmm) = -tenergy{select_h}.n_data(m).tmn(n_pos);
    sigma_m_data(mmm)= tenergy{select_h}.n_data(m).sigma_tmn(n_pos);
    m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
end
hold all;
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH,  m_data./ustar^3 + sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH,  m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
xlim([0 120]);
ann_str = sprintf('n=%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.6,text_y_pos, ann_str,'Units','Normalized');
xlabel('r_{m}/\delta');

subplot(247);
select_n = 5; select_h = 5;
possible_m = m_st:select_n-1;
m_data = zeros([length(possible_m) 1]);
sigma_m_data = zeros([length(possible_m) 1]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for mmm = 1:size(m_data,1)
    m = possible_m(mmm); n_pos = select_n-m;
    m_data(mmm) = -tenergy{select_h}.n_data(m).tmn(n_pos);
    sigma_m_data(mmm)= tenergy{select_h}.n_data(m).sigma_tmn(n_pos);
    m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
end
hold all;
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH,  m_data./ustar^3 + sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH,  m_data./ustar^3 - sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
xlim([0 120]);
ann_str = sprintf('n=%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.6,text_y_pos, ann_str,'Units','Normalized');
xlabel('r_{m}/\delta');

subplot(248);
select_n = 5; select_h = 12;
possible_m = m_st:select_n-1;
m_data = zeros([length(possible_m) 1]);
sigma_m_data = zeros([length(possible_m) 1]);
n_scale = 2.^(nscale-select_n).*(dx*dy)^(0.5);
m_in_m = zeros([length(possible_m) 1]);
for mmm = 1:size(m_data,1)
    m = possible_m(mmm); n_pos = select_n-m;
    m_data(mmm) = -tenergy{select_h}.n_data(m).tmn(n_pos);
    sigma_m_data(mmm)= tenergy{select_h}.n_data(m).sigma_tmn(n_pos);
    m_in_m(mmm) = tenergy{select_h}.m_in_m(m); 
end
hold all;
bounds =  repmat(ones([length(m_data) 1]).*sigma_m_data,[1 2]);
% [l,p]= boundedline(m_in_m./BLH, m_data./ustar^3,...
%         bounds./ustar^3.*0.5,'-or', 'transparency', 0.3);
plot(m_in_m./BLH, m_data./ustar^3,'mo-');  
plot(m_in_m./BLH,  m_data./ustar^3 + sigma_m_data./ustar^3,'c+'); 
plot(m_in_m./BLH,  m_data./ustar^3 -sigma_m_data./ustar^3,'c+'); 
plot([min(m_in_m(:)) max(m_in_m(:))]./BLH, [0 0],'--k');
xlim([0 120]);
ann_str = sprintf('n=%2.2f\\delta',n_scale/BLH);
% text(text_x_pos,text_y_pos, ann_str,'Units','Normalized');
ann_str = sprintf('h=%2.2f\\delta',tenergy{select_h}.h/BLH);
text(0.6,text_y_pos, ann_str,'Units','Normalized');
xlabel('r_{m}/\delta');

SPR = sprintf('%s%s','','tmn_ek02_fixed_n_fixed_h.eps');  
set(gcf, 'Color', 'w');
print(formatEng,SPR);
