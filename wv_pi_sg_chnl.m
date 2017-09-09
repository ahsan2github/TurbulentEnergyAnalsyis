clear; clc; close all;
addpath('/home/ahsan/rwt-master/bin');
frameVec=[9,10,11,12]; %,13,14,15,16,17,18,19];
parentDir = '/media/ahsan/ds00/NeutralChannel_2048x2048x64_free';
outputDir = [parentDir '/output/'];
readinputs(parentDir);
u_star_in = u_star; clear u_star;
BLH = 1500.0;
refLevels = 3:3:62; % total 21 levels
nscale = log(2048)/log(2);
lpassFilter = daubcqf(10,'min');
tic
% create data structure to save the results 
pi_energy = cell([1 length(refLevels)]);
for i = 1:length(pi_energy)
	pi_energy{i} = piClass(nscale, dx, dy, dz);
	pi_energy{i}.h = refLevels(i)*dz;
end
for i = 1:length(pi_energy)
	pi_energy2{i} = piClass_(nscale, dx, dy, dz);
	pi_energy2{i}.h = refLevels(i)*dz;
end

frame_energy_pi = cell([1 length(frameVec)]);
frame_energy_pi2 = cell([1 length(frameVec)]);
%%
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

    fn = [parentDir '/output/w_frame/w_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    ww = fread(fh,'double');
    ww = reshape(ww, [Nx,Ny,Nz]).*u_star_in; 
    fclose(fh);
    
    fn = [parentDir '/output/p_frame/p_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    pp = fread(fh,'double');
    pp = reshape(pp, [Nx,Ny,Nz]);
    pp = interpolate_on_w(pp);
    fclose(fh);

    sfl = '/media/ahsan/ds01/stressFiles';
    fn = [sfl '/chnl_tau12_'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    sgs_uv = fread(fh,'double');
    sgs_uv = reshape(sgs_uv, [Nx,Ny,Nz]);
    fclose(fh);   

    fn = [sfl '/chnl_tau13_'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    sgs_uw = fread(fh,'double');
    sgs_uw = reshape(sgs_uw, [Nx,Ny,Nz]);
    fclose(fh); 

    fn = [sfl '/chnl_tau23_'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    sgs_vw = fread(fh,'double');
    sgs_vw = reshape(sgs_vw, [Nx,Ny,Nz]);
    fclose(fh);  
    [~, S] = mywavedec2(uu(:,:,refLevels(1)), lpassFilter, nscale);
    for ii = 1:length(pi_energy)
        if(refLevels(ii) < 3) 
            error(['Calculation must be done above :' num2str(2) ' verical node !']); 
        end        
        for nn = 1:length(pi_energy{ii}.n)
            disp(['calculating n-scale : ' num2str(pi_energy{ii}.n(nn))]);
            tmpPI = zeros([S(pi_energy{ii}.n(nn)+1,1), S(pi_energy{ii}.n(nn)+1)]);
            tmpPI2 = zeros([S(pi_energy{ii}.n(nn)+1,1), S(pi_energy{ii}.n(nn)+1)]);
            % reconstruct the velocity field from smoothing coefficients and subtract 
            % from the original velocity field
            [wtu, S] = mywavedec2(uu(:,:,refLevels(ii)), lpassFilter, nscale);
            [wtv, ~] = mywavedec2(vv(:,:,refLevels(ii)), lpassFilter, nscale);
            [wtw, ~] = mywavedec2(ww(:,:,refLevels(ii)), lpassFilter, nscale);        
            [wtp, ~] = mywavedec2(pp(:,:,refLevels(ii)), lpassFilter, nscale);
            sgs_uv_lvl = sgs_uv(:,:,refLevels(ii));
            sgs_uw_lvl = sgs_uw(:,:,refLevels(ii));
            sgs_vw_lvl = sgs_vw(:,:,refLevels(ii));

            [wtu_above, ~] = mywavedec2(uu(:,:,refLevels(ii)+1), lpassFilter, nscale);
            [wtv_above, ~] = mywavedec2(vv(:,:,refLevels(ii)+1), lpassFilter, nscale);
            [wtw_above, ~] = mywavedec2(ww(:,:,refLevels(ii)+1), lpassFilter, nscale);        
            [wtp_above, ~] = mywavedec2(pp(:,:,refLevels(ii)+1), lpassFilter, nscale);     
            sgs_uv_above = sgs_uv(:,:,refLevels(ii)+1);
            sgs_uw_above = sgs_uw(:,:,refLevels(ii)+1);
            sgs_vw_above = sgs_vw(:,:,refLevels(ii)+1);                

            uu_gt_n = uu(:,:,refLevels(ii)).*uu(:,:,refLevels(ii)) - ...
                wv_CutGreater_than_n(wtu, S, lpassFilter, pi_energy{ii}.n(nn)).* ...
                wv_CutGreater_than_n(wtu, S, lpassFilter, pi_energy{ii}.n(nn));

            uv_gt_n = uu(:,:,refLevels(ii)).*vv(:,:,refLevels(ii)) - ...
                wv_CutGreater_than_n(wtu, S, lpassFilter, pi_energy{ii}.n(nn)).* ...
                wv_CutGreater_than_n(wtv, S, lpassFilter,pi_energy{ii}.n(nn));

            uw_gt_n = uu(:,:,refLevels(ii)).*ww(:,:,refLevels(ii)) - ...
                wv_CutGreater_than_n(wtu, S, lpassFilter,pi_energy{ii}.n(nn)).* ...
                wv_CutGreater_than_n(wtw, S, lpassFilter,pi_energy{ii}.n(nn));               

            uw_gt_n_above = uu(:,:,refLevels(ii)+1).*ww(:,:,refLevels(ii)+1) - ...
                wv_CutGreater_than_n(wtu_above, S, lpassFilter,pi_energy{ii}.n(nn)).* ...
                wv_CutGreater_than_n(wtw_above, S, lpassFilter,pi_energy{ii}.n(nn));                  

            vv_gt_n = vv(:,:,refLevels(ii)).*vv(:,:,refLevels(ii)) - ...
                wv_CutGreater_than_n(wtv, S, lpassFilter,pi_energy{ii}.n(nn)).* ...
                wv_CutGreater_than_n(wtv, S, lpassFilter,pi_energy{ii}.n(nn));                

            vw_gt_n = vv(:,:,refLevels(ii)).*ww(:,:,refLevels(ii)) - ...
                wv_CutGreater_than_n(wtv, S, lpassFilter,pi_energy{ii}.n(nn)).* ...
                wv_CutGreater_than_n(wtw, S, lpassFilter,pi_energy{ii}.n(nn)); 

            vw_gt_n_above = vv(:,:,refLevels(ii)+1).*ww(:,:,refLevels(ii)+1) - ...
                wv_CutGreater_than_n(wtv_above, S, lpassFilter,pi_energy{ii}.n(nn)).* ...
                wv_CutGreater_than_n(wtw_above, S, lpassFilter,pi_energy{ii}.n(nn));         

            ww_gt_n = ww(:,:,refLevels(ii)).*ww(:,:,refLevels(ii)) - ...
                wv_CutGreater_than_n(wtw, S, lpassFilter,pi_energy{ii}.n(nn)).* ...
                wv_CutGreater_than_n(wtw, S, lpassFilter,pi_energy{ii}.n(nn)); 

            ww_gt_n_above = ww(:,:,refLevels(ii)+1).*ww(:,:,refLevels(ii)+1) - ...
                wv_CutGreater_than_n(wtw_above, S, lpassFilter, pi_energy{ii}.n(nn)).* ...
                wv_CutGreater_than_n(wtw_above, S, lpassFilter, pi_energy{ii}.n(nn)); 

            p_gt_n = pp(:,:,refLevels(ii)) - ...
                wv_CutGreater_than_n(wtp, S, lpassFilter, pi_energy{ii}.n(nn));
            p_gt_n_above = pp(:,:,refLevels(ii)+1) - ...
                wv_CutGreater_than_n(wtp_above, S, lpassFilter, pi_energy{ii}.n(nn));                               

            [ddx_p_gt_n, ddy_p_gt_n] = ddx(p_gt_n, z_i,l_r,Nx);          
            ddz_p_gt_n = ddz_w_wn(p_gt_n, dz, p_gt_n_above);

            [ddx_uu_gt_n, ~] = ddx(uu_gt_n, z_i,l_r,Nx);
            [ddx_uv_gt_n, ddy_uv_gt_n] = ddx(uv_gt_n, z_i,l_r, Nx);
            [~, ddy_vv_gt_n] = ddx(vv_gt_n, z_i,l_r, Nx);
            [ddx_uw_gt_n, ~] = ddx(uw_gt_n, z_i, l_r, Nx);
            [~, ddy_vw_gt_n] = ddx(vw_gt_n, z_i, l_r, Nx);

            [ddx_sgs_uv, ddy_sgs_uv] = ddx(sgs_uv_lvl, z_i,l_r,Nx);
            [ddx_sgs_uw, ~] = ddx(sgs_uw_lvl, z_i,l_r,Nx);
            [~, ddy_sgs_vw] = ddx(sgs_vw_lvl, z_i, l_r, Nx);

            ddz_uw_gt_n = ddz_w_wn(uw_gt_n, dz, uw_gt_n_above);
            ddz_vw_gt_n = ddz_w_wn(vw_gt_n, dz, vw_gt_n_above);
            ddz_ww_gt_n = ddz_w_wn(ww_gt_n, dz, ww_gt_n_above);

            ddz_sgs_uw = ddz_w_wn(sgs_uw_lvl, dz,sgs_uw_above);
            ddz_sgs_vw = ddz_w_wn(sgs_vw_lvl, dz,sgs_vw_above);

            [wv_adv_u_2nd_factor, Sadv] = mywavedec2(ddx_uu_gt_n + ...
                                  ddy_uv_gt_n + ddz_uw_gt_n + ddy_sgs_uv + ...
                                  ddz_sgs_uw + ddx_p_gt_n, lpassFilter, nscale);
            [wv_adv_v_2nd_factor, ~] = mywavedec2(ddx_uv_gt_n + ...
                                  ddy_vv_gt_n + ddz_vw_gt_n + ddx_sgs_uv + ...
                                  ddz_sgs_vw + ddy_p_gt_n, lpassFilter, nscale);
            [wv_adv_w_2nd_factor, ~] = mywavedec2(ddx_uw_gt_n + ...
                                  ddy_vw_gt_n + ddz_ww_gt_n + ddx_sgs_uw + ...
                                  ddy_sgs_vw + ddz_p_gt_n, lpassFilter, nscale);            
            avgval = 0;
            sq_avgval = 0;
            for mm = 1:length(pi_energy{ii}.m_data(nn).m)         
                tmpRes = (wv_get_coeff_n_level(wtu, S, pi_energy{ii}.m_data(nn).m(mm)).*...
                    wv_get_coeff_n_level(wv_adv_u_2nd_factor, S, pi_energy{ii}.m_data(nn).m(mm)) + ...
                    wv_get_coeff_n_level(wtv, S, pi_energy{ii}.m_data(nn).m(mm)).* ...
                    wv_get_coeff_n_level(wv_adv_v_2nd_factor, S, pi_energy{ii}.m_data(nn).m(mm)) +...
                    wv_get_coeff_n_level(wtw, S, pi_energy{ii}.m_data(nn).m(mm)).*...
                    wv_get_coeff_n_level(wv_adv_w_2nd_factor,S, pi_energy{ii}.m_data(nn).m(mm)));                 
                
                disp(['m: ' num2str(pi_energy{ii}.m_data(nn).m(mm)) ' ::  Size(tmpRes) :' num2str(size(tmpRes))]);
                sum3coeff =  sum(tmpRes,3);
                normScale = 2^((nscale-pi_energy{ii}.m_data(nn).m(mm)))*(dx*dy)^(0.5);
                pi_energy{ii}.m_data(nn).tmn(mm) = mean(sum3coeff(:))/2/pi/log(2)/normScale;                
								if(pi_energy{ii}.n(nn) == pi_energy{ii}.m_data(nn).m(mm))
									sign = 1;
								else
									sign = 1;
								end
                tmpPI = tmpPI + sign.*stretch_by_zero(sum3coeff,2, ...
                                            pi_energy{ii}.n(nn)-pi_energy{ii}.m_data(nn).m(mm));   
                % excludes the case n=m
                if(mm < length(pi_energy{ii}.m_data(nn).m))
                  pi_energy2{ii}.m_data(nn).tmn(mm) = mean(sum3coeff(:))/2/pi/log(2)/normScale;
                  tmpPI2 = tmpPI2 + stretch_by_zero(sum3coeff,2, ...
                                            pi_energy{ii}.n(nn)-pi_energy{ii}.m_data(nn).m(mm)); 
                  disp(['n: ' num2str(pi_energy{ii}.n(nn)) '  m: ' num2str(pi_energy{ii}.m_data(nn).m(mm))]);
                end      
            end 
            % calculate pi_sg
            normScale = 2^((nscale-pi_energy{ii}.n(nn)))*(dx*dy)^(0.5);
            pi_energy{ii}.pi_energy(nn) = mean(tmpPI(:));
            pi_energy{ii}.pi_energy_sigma(nn) = ((mean(mean(tmpPI.^2))-(mean(tmpPI(:)))^2)^0.5);		
            pi_energy2{ii}.pi_energy(nn) = mean(tmpPI2(:));
            pi_energy2{ii}.pi_energy_sigma(nn) = ((mean(mean(tmpPI2.^2))-(mean(tmpPI2(:)))^2)^0.5);																									
            disp(['--------------------------------------------------------------------------------- ']);
            disp([' ']);
        end
        disp(['Calculated at Height : ' num2str(refLevels(ii))]);
        disp(['===================================================================================']);
        disp([' ']);
    end
    frame_energy_pi{ff} = pi_energy;
    frame_energy_pi2{ff} = pi_energy2;
    disp(['Calculated frame : ' num2str(frameVec(ff))]);
		save('wav_pi_tmn_chnl.mat', 'pi_energy', 'pi_energy2','frame_energy_pi','frame_energy_pi2','-v7.3');
    disp(['v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v']);
    disp([' ']);
end
save('wav_pi_tmn_chnl.mat', 'pi_energy','pi_energy2','frame_energy_pi','frame_energy_pi2','-v7.3');
toc
