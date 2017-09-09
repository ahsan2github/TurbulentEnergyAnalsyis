%%% this function gives the wavelet coefficients of the n-th level of
%%% decomposition. Assuming, n-th level is calculated begining from larger scales
%%% and progressing towards the smaller scales. Level 1 will have the lowest resolution
%%% Nx/2^(n-1) x Ny/2^(n-1). level 2 will have resolution Nx/2^(n-2) x Ny/2^(n-2). level 3 will
%%% have resolution Nx/2^(n-3) x Ny/2^(n-3) and last level will have res. Nx/2 x Ny/2.
%%% In other words, it gives coefficients corresponding to
%%% the book keeping matrix row length(S)-n
function [out, varargout] = wv_get_coeff_n_level(C, S, n,varargin)
    % determine the level of decomposition
    lvlDecomp = length(S)-2;
    if(n>lvlDecomp), error(['Requested coefficients do not exist, increase the decomposition level']); end;
    if(n>1)
      st = S(1, 1) * S(1, 2) + sum(S(2:n, 1) .* S(2:n, 2) * 3) + 1;
      ed = S(1, 1) * S(1, 2) + sum(S(2:n + 1, 1) .* S(2:n + 1, 2) * 3);
    else
      st = S(1, 1) * S(1, 2) + 1;
      ed = S(1, 1) * S(1, 2) + sum(S(2, 1) .* S(2, 2) * 3);
    end
      out = reshape(C(st:ed),[S(n+1,1), S(n+1,2),3]);
      varargout{1}= st;
      varargout{2}= ed;  
%       C(1:st-1) = 0.0; C(ed+1:end) = 0.0;
%       varargout{3}= waverec2(C, S, varargin{1});   
end
