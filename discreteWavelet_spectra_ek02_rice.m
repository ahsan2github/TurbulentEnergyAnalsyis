clc; clear;
addpath('/raid/home/mohammad/rwt-master/bin');
lpassFilter = daubcqf(10,'min');
%chnl 
frameVec = [38,39,40,41];
heightVec = 2:4:40;
parentDir = '/raid/home/mohammad/Ekman_Ugeo2_2048x2048x64_Lz750';
outputDir = [parentDir '/output/'];
readinputs(parentDir);
u_star_in = u_star; clear u_star;
ustar = 0.09;
BLH = 550.69;
M = log2(Nx);
scale = 1:M-1;
rm = 2.^(M-scale).*dx;
kwave = 2*pi*z_i./rm;
kwave_delta = kwave.*(BLH/z_i);
kx_fourier = (0:1023)./z_i;
kx_del_fourier = kx_fourier(2)-kx_fourier(1);
kx_fourier = kx_fourier(2:end-1);
kx_delta_fourier = kx_fourier.*(BLH);
wavelet_spectra_u = zeros([length(scale) 1]);
fourier_spectra_u = zeros([Nx/2-2 1]);
wavelet_sigma_u = zeros([length(scale) 1]);
disp( [ 'kx_delta_fourier(1:10): ' num2str(kx_delta_fourier(1:10)) ] );
disp( [ 'kwave_delta(1:10): ' num2str(kwave_delta(1:10)) ] );
wavelet_spectra_u_frame = zeros([length(heightVec) length(scale)]);
fourier_spectra_u_frame = zeros([length(heightVec) Nx/2-2]);
wavelet_sigma_u_frame = zeros([length(heightVec) length(scale)]);
for ff = 1:length(frameVec)
    % load frames
    frameStr = sprintf('%4.4i',frameVec(ff));
    fn = [parentDir '/output/u_frame/u_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    uu = fread(fh,'double');
    uu = reshape(uu, [Nx,Ny,Nz]).*u_star_in + Ugal; 
    uu = interpolate_on_w(uu);
    fclose(fh);
    [dumm_decomp_wv, bk] = mywavedec1(uu(:,1,1),lpassFilter,scale(end));
    for h = 1:length(heightVec)
        tmp_fft = zeros([Nx Ny]);
        for jj = 1:Ny
            tmp_fft(:,jj) = fft(uu(:,jj,heightVec(h))).*conj(fft(uu(:,jj,heightVec(h))));
        end
        E_w_fft = 0.5 .*mean(tmp_fft./Nx^2,2);
        fourier_spectra_u = E_w_fft(2:Nx/2-1); 
        tmp_u = uu(:,:,heightVec(h));% - mean(mean(uu(:,:,heightVec(h))));
        tmp_trans = zeros([length(dumm_decomp_wv) Ny]);
        for jj = 1:Ny
            [tmp_trans(:,jj), bk] = mywavedec1(tmp_u(:,jj),lpassFilter,scale(end));
        end
        tmp_trans_streamwise = 0.5.* mean(tmp_trans.^2,2);
        for ss = 1:length(scale)
            tmp = wv_get_coeff_n_level_1d(tmp_trans_streamwise, bk, scale(ss));
            E_w = (dx/2/pi/z_i/log(2))* mean(tmp(:));
            wavelet_spectra_u(ss) = E_w; 
            wavelet_sigma_u(ss) = (dx/2/z_i/pi/log(2))*...
            sqrt(mean(tmp(:).^2)- mean(tmp(:))^2);
        end        
        wavelet_spectra_u_frame(h, :) = wavelet_spectra_u_frame(h, :) + wavelet_spectra_u';
        fourier_spectra_u_frame(h, :) = fourier_spectra_u_frame(h, :) + fourier_spectra_u';
        wavelet_sigma_u_frame(h, :) = wavelet_sigma_u_frame(h, :) + wavelet_sigma_u';        
    end

end
wavelet_spectra_u_frame = wavelet_spectra_u_frame./length(frameVec);
fourier_spectra_u_frame = fourier_spectra_u_frame./length(frameVec);
wavelet_sigma_u_frame = wavelet_sigma_u_frame./length(frameVec);

disp( [ 'wavelet_spectra_u_frame(1:10): ' num2str(wavelet_spectra_u_frame(1:10)) ] );
disp( [ 'fourier_spectra_u_frame(1:10): ' num2str(fourier_spectra_u_frame(1:10)) ] );

% fourier_spectra_u_intrp = interp1(km_delta_fourier,fourier_spectra_u_frame(2,:),km_delta);
%%
LineStyles = {'-', '--','-.',':','-+','-o','-s','-d','-^','-*'}; % 10 types total
figWidth = 2; figheight = figWidth/1;  formatEng = '-depsc2';
figure(); set(gcf,'Renderer','painters');
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figheight]);
set(gcf,'PaperPositionMode','auto');
hold all;
counter = 0;
for kk = 1:2:length(heightVec)
    counter = counter + 1;
    plot(kx_delta_fourier, fourier_spectra_u_frame(kk,:)./ustar^2./BLH, LineStyles{counter}); hold on;
    plot(kwave_delta, wavelet_spectra_u_frame(kk,:)./ustar^2./BLH,'o'); hold on;
%     plot(kwave_delta, wavelet_spectra_u_frame(kk,:)+ wavelet_sigma_u_frame(kk,:),'+'); 
    % xticks([0.1 1 10 100]);
    % yticks([0.00001 0.0001 0.001 0.01]);
    xlabel('$k_x \delta$','Interpreter','Latex','FontSize',13);
    ylabel('$E_{11}(k_x)u_*^{-2}\delta^{-1}$','Interpreter','Latex','FontSize',13);
end
xlim([0.05 50]);
ylim([8e-9 0.00004]);
xx = [10 40];
loglog(xx,xx.^(-5/3).*0.0001,'LineWidth',2); % -5/3 slope
text(18, 0.0000012, '-5/3','FontSize',10);
xx = [1.5 9];
loglog(xx,xx.^(-1).*0.00002,'LineWidth',2); % -1 slope
text(3, 0.000009, '-1','FontSize',10);
text(0.1,0.1,'(c)','Units','Normalized','FontSize',11);
set(gca,'XScale','log'); set(gca,'YScale','log');
SPR = sprintf('%s%s','','wavelet_n_fourier_spectra_ek02.eps');  
set(gcf, 'Color', 'w');
print(formatEng,SPR);