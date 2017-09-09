function [out] = stretch_by_zero(arr, scale, exponent)
    flag = false;
    if(mod(size(arr,2),scale)), flag = true; end
    if(mod(size(arr,1),scale)), flag = true; end
    if (size(arr,1)==1 && size(arr,2)==1), flag = false; end       
    if(flag)
      error(['Target Grid resolution in row direction is not a multiple of: ' num2str(scale)]); 
    end
    if(size(arr,1)==1 && size(arr,2)==1), scale = 2; end
    out = zeros([size(arr,1)*scale^exponent size(arr,2)*scale^exponent]);
    out(scale^exponent:scale^exponent:size(out,1), scale^exponent:scale^exponent:size(out,2)) = arr;    
end
