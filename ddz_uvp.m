
%  First deriv in z direction for boundary layer (2nd order numerics)
%  F is on uvp nodes and dFdz is on w nodes

function [dfdz] = ddz_uvp(f,dz)
    dfdz = zeros([size(f,1), size(f,2), size(f,3)]);
    for k=2:size(f,3)
        dfdz(:,:,k)=(f(:,:,k)-f(:,:,k-1))./dz;
    end 
end
