clc; clear;
addpath('/raid/home/mohammad/rwt-master/bin');
lpassFilter = daubcqf(10,'min');
%chnl 
frameVec = [9:19];
heightVec = 2:4:40;
parentDir = '/raid/home/mohammad/NeutralChannel_2048x2048x64_free';
outputDir = [parentDir '/output/'];
readinputs(parentDir);
z_vec = ((1:Nz)-1).*dz;
u_star_in = u_star; clear u_star;
ustar = 0.42;
BLH = 1500.0;
M = log2(Nx);
scale = 1:M-1;
rm = 2.^(M-scale).*dx;
kwave = 2*pi*z_i./rm;
kx_fourier = (0:1023);
kx_fourier = kx_fourier(2:end-1);
kx_delta_fourier = kx_fourier.*(BLH);
wavelet_spectra_u = zeros([length(scale) 1]);
fourier_spectra_u = zeros([Nx/2-2 1]);
wavelet_sigma_u = zeros([length(scale) 1]);
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
    
    fn = [parentDir '/output/v_frame/v_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    vv = fread(fh,'double');
    vv = reshape(vv, [Nx,Ny,Nz]).*u_star_in + Vgal; 
    vv = interpolate_on_w(vv);
    
    fclose(fh);
    [dumm_decomp_wv, bk] = mywavedec1(uu(:,1,1),lpassFilter,scale(end));
    for h = 1:length(heightVec)      
        tmp_fft = zeros([Nx Ny]);
        for jj = 1:Ny
            tmp_fft(:,jj) = fft(uu(:,jj,heightVec(h))).*conj(fft(uu(:,jj,heightVec(h))));
        end
        E_w_fft = mean(tmp_fft./Nx^2,2);
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
LineStyles = {'-c', '--m','-.r',':g','-+b','-oc','-sm','-d','-^','-*'}; % 10 types total
LineStyles1 = {'oc', 'om','or','og','ob','oc','-s','-d','-^','-*'}; % 10 types total
figWidth = 3.2; figheight = figWidth/1;  formatEng = '-depsc2';
hf = figure(); 
set(gcf,'Renderer','painters');
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figWidth figheight]);
set(hf,'Visible','off');
hold all;
counter = 0;
for kk = 1:2:length(heightVec)
    counter = counter + 1;
    %plot(kx_fourier.*(BLH/z_i), fourier_spectra_u_frame(kk,:).*kx_fourier./ustar^2, LineStyles{counter}); 
    plot(kx_fourier.*z_vec(heightVec(kk))/z_i, fourier_spectra_u_frame(kk,:).*kx_fourier./ustar^2, LineStyles{counter}); 
    %plot(kwave.*(BLH/z_i), wavelet_spectra_u_frame(kk,:).*kwave./ustar^2,LineStyles1{counter}); 
    plot(kwave.*z_vec(heightVec(kk))/z_i, wavelet_spectra_u_frame(kk,:).*kwave./ustar^2,LineStyles1{counter}); 
    % xticks([0.1 1 10 100]);
    % yticks([0.00001 0.0001 0.001 0.01]);
end
xlabel('$k_x z$','Interpreter','Latex','FontSize',10);
ylabel('$k_x E_{11}(k_x)u_*^{-2}$','Interpreter','Latex','FontSize',10);
colormap parula;
%xlim([0.05 100]);
% ylim([2e-7 0.00003]);
xx = [3 60];
loglog(xx,xx.^(-2/3).*0.79,'LineWidth',2); % -2/3 slope
text(11.5, 0.175, '-2/3','FontSize',8);
set(gca,'XScale','log'); set(gca,'YScale','log');
box on;
SPR = sprintf('%s%s','','wavelet_n_fourier_spectra_chnl_pre.eps');  
set(gcf, 'Color', 'w');
print(formatEng,SPR);
