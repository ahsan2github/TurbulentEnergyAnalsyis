
%...F is on w nodes and dFdz is on uvp nodes
function [dfdz] = ddz_w(f,dz)
    dfdz = zeros([size(f,1), size(f,2), size(f,3)]); % 2:Nz  uvp nodes    
    for k=1:size(f,3)-1
        dfdz(:,:,k)=(f(:,:,k+1)-f(:,:,k))./dz;
    end
    dfdz(:,:, size(f,3)) = zeros([size(f,1) size(f,2)]);
end 
