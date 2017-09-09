function [horizontal_shear_stress, uw,uv, vw, u_star_mean,...
 tau13, tau12, tau23]= ...
   getStress3d_(parentDir, Nx, Ny, Nz, dx, dy, dz, u_star, Ugal, Vgal,zo, ...
   l_r, z_i,frameVec)

uw = zeros([Nx, Ny, Nz]);
uv = zeros([Nx, Ny, Nz]);
vw = zeros([Nx, Ny, Nz]);


tau12 = zeros([Nx, Ny, Nz]); 
tau13 = zeros([Nx, Ny, Nz]);
tau23 = zeros([Nx, Ny, Nz]);     
u_mean_nz = zeros([Nz, 1]);
v_mean_nz = zeros([Nz, 1]);
w_mean_nz = zeros([Nz, 1]);
tau_reynolds_s = zeros([Nx Ny]);
u_star_mean = 0.0;
horizontal_shear_stress = zeros([Nx, Ny, Nz]);
delta = (dx*dy*dz)^(1.0/3.0);
framecnt = 0;

for frameNo = frameVec
    frameStr = sprintf('%4.4i',frameNo);
    fn = [parentDir '/output/u_frame/u_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    uu = fread(fh,'double');
    uu = reshape(uu, [Nx,Ny,Nz]).*u_star + Ugal;     
    fclose(fh);    

    fn = [parentDir '/output/v_frame/v_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    vv = fread(fh,'double');
    vv = reshape(vv, [Nx,Ny,Nz]).*u_star + Vgal; 
    fclose(fh);     

    fn = [parentDir '/output/w_frame/w_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    ww = fread(fh,'double');
    ww = reshape(ww, [Nx,Ny,Nz]).*u_star; 
    fclose(fh); 

    fn = [parentDir '/output/Cs2_frame/Cs2_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    Cs2 = fread(fh,'double');
    Cs2 = reshape(Cs2, [Nx,Ny,Nz]); 
    fclose(fh); 
    % calculate mean at u,v & w nodes    
    for jj = 1:Nz                
        u_mean_nz(jj) = mean(mean(uu(:,:,jj)));
        v_mean_nz(jj) = mean(mean(vv(:,:,jj)));
        % at w-node
        w_mean_nz(jj) = mean(mean(ww(:,:,jj)));
    end
    uw_nxnynz = zeros([Nx, Ny, Nz]);
    uv_nxnynz = zeros([Nx, Ny, Nz]);
    vw_nxnynz = zeros([Nx, Ny, Nz]);
    % calculate reynolds stress at uvp nodes, dimension [Nx, Ny, Nz]      
    for i = 1:Nz-1
        uw_nxnynz(:,:,i) = uw_nxnynz(:,:,i) + (uu(:,:,i)-u_mean_nz(i)).*...
          ( 0.5.*(ww(:,:,i+1)+ww(:,:,i))- 0.5.*(w_mean_nz(i+1)+w_mean_nz(i)) );
        uv_nxnynz(:,:,i) = uv_nxnynz(:,:,i) + (uu(:,:,i)-u_mean_nz(i)).*(vv(:,:,i)-v_mean_nz(i));
        vw_nxnynz(:,:,i) = vw_nxnynz(:,:,i) + (vv(:,:,i)-v_mean_nz(i)).*...           
          ( 0.5.*(ww(:,:,i+1)+ww(:,:,i))- 0.5.*(w_mean_nz(i+1)+w_mean_nz(i)) );
    end
    % calculate the surface shear stress from actual velocity field,
    % similarity theory
    u_r = sqrt((uu(:,:,1)).^2+vv(:,:,1).^2);
    tau_xzs = -(u_r(:,:)./log(dz/2/zo)).^2 .*(uu(:,:,1)./...
        u_r(:,:));
    tau_yzs= -(u_r(:,:)./log(dz/2/zo)).^2 .*(vv(:,:,1)./...
        u_r(:,:));

    tau_reynolds_s = tau_reynolds_s + (tau_xzs.^2 + tau_yzs.^2).^0.5;
    u_star_mean = u_star_mean + mean(mean(((tau_xzs.^2 + tau_yzs.^2).^0.25)));
    % set the stress at ground zero
    uw(:,:,1) = uw(:,:,1) + zeros([Nx Ny]);
    uv(:,:,1) = uv(:,:,1) + zeros([Nx Ny]);
    vw(:,:,1) = vw(:,:,1) + zeros([Nx Ny]);
    % interpolate uw, uv, vw at w-nodes, 2:Nz w nodes
    uw(:,:,2:Nz) = uw(:,:,2:Nz) + 0.5.*(uw_nxnynz(:,:,1:Nz-1) + uw_nxnynz(:,:,2:Nz));
    uv(:,:,2:Nz) = uv(:,:,2:Nz) + 0.5.*(uv_nxnynz(:,:,1:Nz-1) + uv_nxnynz(:,:,2:Nz));
    vw(:,:,2:Nz) = vw(:,:,2:Nz) + 0.5.*(vw_nxnynz(:,:,1:Nz-1) + vw_nxnynz(:,:,2:Nz));
    clearvars uw_nxnynz uv_nxnynz vw_nxnynz;
    % calculate gradients ddx, ddy of  u,v at u-v-p nodes
    [dudx,dudy]=ddx(uu,z_i,l_r,Nx);
    [dvdx,dvdy]=ddx(vv,z_i,l_r,Nx);
    % interpolate ddx, ddy of u,v at w nodes (2:Nz)
    dudx(:,:,2:Nz) = 0.50 .* (dudx(:,:,1:Nz-1) + dudx(:,:,2:Nz)); dudx(:,:,Nz) = zeros([Nx Ny]);
    dudy(:,:,2:Nz) = 0.50 .* (dudy(:,:,1:Nz-1) + dudy(:,:,2:Nz)); dudy(:,:,Nz) = zeros([Nx Ny]);
    dvdx(:,:,2:Nz) = 0.50 .* (dvdx(:,:,1:Nz-1) + dvdx(:,:,2:Nz)); dvdx(:,:,Nz) = zeros([Nx Ny]);
    dvdy(:,:,2:Nz) = 0.50 .* (dvdy(:,:,1:Nz-1) + dvdy(:,:,2:Nz)); dvdy(:,:,Nz) = zeros([Nx Ny]);
    % set x, y derivatives at surface to zero
    dudx(:,:,1) = zeros([Nx Ny]);
    dudy(:,:,1) = zeros([Nx Ny]);
    dvdx(:,:,1) = zeros([Nx Ny]);
    dvdy(:,:,1) = zeros([Nx Ny]);

    % calculate ddx of w at w-nodes (2:Nz)
    [dwdx,dwdy] = ddx(ww,z_i,l_r,Nx); 
    % set dwdx at surface to zero
    dwdx(:,:,1) = zeros([Nx Ny]);
    % calculate ddz of u, v at w nodes (2:Nz)
    dudz = ddz_uvp(uu, dz); dvdz = ddz_uvp(vv, dz);
    % calculate ddz of u,v at surface  aka 1st w node using forward difference of 1st order accuracy
    dudz(:,:,1) = (uu(:,:,1))/dz/2;	
    dvdz(:,:,1) = (vv(:,:,1))/dz/2;	
    % calculate ddz of w at uvp nodes, will have a vertical dimension Nz-1
    dwdz = ddz_w(ww, dz); 
    % interpolate dwdz at w-nodes, this will leave dwdz with a vertical
    % dimension of Nz-1, 2:Nz
    dwdz(:,:,2:Nz) = 0.5.*(dwdz(:,:,1:Nz-1)+dwdz(:,:,2:Nz)); 
    % set w derivative at ground, 1st order accurate !!
    dwdz(:,:,1) = ww(:,:,2)/dz;          
    % calculate strain rates at w nodes, will have vertical dimensiom
    % of Nz, on w-nodes
    %disp(['Size of dudx, dudy, dudz : ' num2str(size(dudx)) '   ' num2str(size(dudy)) '     ' num2str(size(dudz))]);
    %disp(['Size of dvdx, dvdy, dvdz : ' num2str(size(dvdx)) '   ' num2str(size(dvdy)) '     ' num2str(size(dvdz))]);
    %disp(['Size of dwdx, dwdy, dwdz : ' num2str(size(dwdx)) '   ' num2str(size(dwdy)) '     ' num2str(size(dwdz))]);
    S11 = 0.5.*(dudx(:,:,1:end) + dudx(:,:,1:end)); 
    S12 = 0.5.*(dudy(:,:,1:end) + dvdx(:,:,1:end));    
    S13 = 0.5.*(dudz(:,:,1:end) + dwdx(:,:,1:end)); 
    S22 = 0.5.*(dvdy(:,:,1:end) + dvdy(:,:,1:end)); 
    S23 = 0.5.*(dvdz(:,:,1:end) + dwdy(:,:,1:end)); 
    S33 = 0.5.*(dwdz(:,:,1:end) + dwdz(:,:,1:end));
    abs_S = sqrt(2.0*(S11.*S11 + S12.*S12 + S13.*S13 + S12.*S12 + S22.*S22 + ...
            S23.*S23 + S13.*S13 + S23.*S23 + S33.*S33));
	
	  clearvars dudx dudy dudz dvdx dvdy dvdz dwdx dwdy dwdz;
    % at this point all quantities are available at w-nodes 1:Nz   
    % if mom_nodes == 0 in LESinputs, the first CS2 value is calculated
    % at dz/2 location, then subsequent Cs values are calculated at
    % regular w-nodes starting from 2nd w-node
    tau12(:,:,1) = tau12(:,:,1) + zeros([Nx, Ny]); 
		tau13(:,:,1) = tau13(:,:,1) + tau_xzs; 
    tau23(:,:,1) = tau23(:,:,1) + tau_yzs;
    tau12(:,:,2:Nz) = tau12(:,:,2:Nz) + ( -2.0 .* (delta .* sqrt(Cs2(:,:,2:Nz))).^2 .* abs_S(:,:,2:Nz) .* S12(:,:,2:Nz) );
    tau13(:,:,2:Nz) = tau13(:,:,2:Nz) + ( -2.0 .* (delta .* sqrt(Cs2(:,:,2:Nz))).^2 .* abs_S(:,:,2:Nz) .* S13(:,:,2:Nz) );
    tau23(:,:,2:Nz) = tau23(:,:,2:Nz) + ( -2.0 .* (delta .* sqrt(Cs2(:,:,2:Nz))).^2 .* abs_S(:,:,2:Nz) .* S23(:,:,2:Nz) );
%     disp(['tau_13 vertical profile: ' num2str((squeeze(mean(mean(tau13))))')]);
%     figure();
%     plot(squeeze(mean(mean(tau13))), (0:63).*dz );
%     xlabel('tau_{13}(m^2/s^2)'); ylabel('z(m)');	
    clearvars S11 S12 S13 S22 S23 S33;
    %disp(['Size of dudx, dudy, dudz : ' num2str(size(dudx)) '   ' num2str(size(dudy)) '     ' num2str(size(dudz))]);
    %disp(['Size of dvdx, dvdy, dvdz : ' num2str(size(dvdx)) '   ' num2str(size(dvdy)) '     ' num2str(size(dvdz))]);
    %disp(['Size of dwdx, dwdy, dwdz : ' num2str(size(dwdx)) '   ' num2str(size(dwdy)) '     ' num2str(size(dwdz))]);
    % add turb stress to resolved stress, stresses are at 2:Nz-1 w-nodes
    horizontal_shear_stress(:,:,1) = horizontal_shear_stress(:,:,1)+...
        (tau13(:,:,1).^2 + tau23(:,:,1).^2).^(0.5);
    horizontal_shear_stress(:,:,:) = horizontal_shear_stress(:,:,:) + ...
        ((uw(:,:,:) + tau13(:,:,:)).^2 + ...
        (vw(:,:,:) + tau23(:,:,:)).^2).^(0.5);
framecnt = framecnt + 1;    
end
u_star_mean = u_star_mean/framecnt;
uw = uw./framecnt;
uv = uv./framecnt;
vw = vw./framecnt;
tau12 = tau12./framecnt;
tau13 = tau13./framecnt;
tau23 = tau23./framecnt;
horizontal_shear_stress = horizontal_shear_stress/framecnt;
end

