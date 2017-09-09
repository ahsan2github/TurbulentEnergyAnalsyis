clc; clear;
addpath('/raid/home/mohammad/rwt-master/bin');
lpassFilter = daubcqf(10,'min');
%chnl 
frameVec = [32,33,34,35,36,37,38,39,40,41];
heightVec = 2:4:40;
parentDir = '/raid/home/mohammad/Ekman_Ugeo2_2048x2048x64_Lz750';
outputDir = [parentDir '/output/'];
readinputs(parentDir);
z_vec = ((1:Nz)-1).*dz;
u_star_in = u_star; clear u_star;
ustar = 0.09;
BLH = 550.69;
M = log2(Nx);
scale = 1:M-1;
rm = 2.^(M-scale).*dx;
kwave = 2*pi*z_i./rm;
kwave_delta = kwave.*(BLH/z_i);
kx_fourier = (0:1023);
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
    fn = [parentDir '/output/v_frame/v_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    vv = fread(fh,'double');
    vv = reshape(vv, [Nx,Ny,Nz]).*u_star_in + Vgal; 
    vv = interpolate_on_w(vv);
    fclose(fh);
    [dumm_decomp_wv, bk] = mywavedec1(uu(:,1,1),lpassFilter,scale(end));
    for h = 1:length(heightVec)
        avgAngle = atan(mean(mean(vv(:,:,heightVec(h))))./mean(mean(uu(:,:,heightVec(h)))));
          % clock wise rotation matrix
        rotMat = [cos(avgAngle)   sin(avgAngle);...
                -sin(avgAngle)  cos(avgAngle)]; 
        uup = zeros([Nx Ny]);
        vvp = zeros([Nx Ny]);       
        for kk = 1:size(uup,1)
            for mm = 1:size(uup,2)
                res = rotMat*[uu(kk,mm,heightVec(h)); vv(kk,mm,heightVec(h))];
                uup(kk,mm) = res(1);
                vvp(kk,mm) = res(2); 
            end
        end                            
        tmp_fft = zeros([Nx Ny]);
        for jj = 1:Ny
            tmp_fft(:,jj) = fft(uup(:,jj)).*conj(fft(uup(:,jj)));
        end
        E_w_fft = mean(tmp_fft./Nx^2,2);
        fourier_spectra_u = E_w_fft(2:Nx/2-1); 
        tmp_trans = zeros([length(dumm_decomp_wv) Ny]);
        for jj = 1:Ny
            [tmp_trans(:,jj), bk] = mywavedec1(uup(:,jj),lpassFilter,scale(end));
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
LineStyles = {'-c', '--m','-.r',':g','-+b','-oc','-s','-d','-^','-*'}; % 10 types total
LineStyles1 = {'oc', 'om','or','og','ob','oc','om','-d','-^','-*'}; % 10 types total
figWidth = 3.2; figheight = figWidth/1;  formatEng = '-depsc2';
set(groot,'defaultFigurePaperPositionMode','manual');
hf = figure(); 
ax = gca;
set(gcf,'Renderer','painters'); box on;
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, figWidth, figheight]);
set(hf, 'Visible', 'off');
hold all;
counter = 0;
for kk = 1:2:length(heightVec)
    counter = counter + 1;
    %plot(kx_fourier.*(BLH/z_i), fourier_spectra_u_frame(kk,:).*kx_fourier./ustar^2, LineStyles{counter}); 
    plot(kx_fourier.*z_vec(heightVec(kk))/z_i, fourier_spectra_u_frame(kk,:).*kx_fourier./ustar^2, LineStyles{counter}); 
    %plot(kwave.*(BLH/z_i), wavelet_spectra_u_frame(kk,:).*kwave./ustar^2,LineStyles1{counter}); 
    plot(kwave.*z_vec(heightVec(kk))/z_i, wavelet_spectra_u_frame(kk,:).*kwave./ustar^2,LineStyles1{counter}); 
end
xlabel('$k_x z$','Interpreter','Latex','FontSize',10);
ylabel('$k_x E_{11}(k_x)u_*^{-2}$','Interpreter','Latex','FontSize',10);
xlim([0.0003 100]);
xx = [0.3 19];
loglog(xx,xx.^(-2/3).*1,'LineWidth',2); % -2/3 slope
text(2, 1, '-2/3','FontSize',8);
set(gca,'XScale','log'); set(gca,'YScale','log');
SPR = sprintf('%s','wavelet_n_fourier_spectra_ek02_pre.eps');  
set(gcf, 'Color', 'w');
print(SPR, formatEng);
