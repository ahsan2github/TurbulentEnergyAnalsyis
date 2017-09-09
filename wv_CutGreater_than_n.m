%%% this function gives large scale field after wavelet filtering
%%% it keeps the large scale 1:n-1 levels and discards all other smaller
%%% scales ie. 'length(S)-2-n' levels. In other words it keeps the coeffiecients corresponding to
%%% the book keeping matrix S's first 1+n rows

function [reconField] = wv_CutGreater_than_n(C, S, lpf, n)
    % determine the level of decomposition
    C = reshape(C,[1 numel(C)]);
    lvlDecomp = length(S)-2;
    if(n>lvlDecomp), error(['Requested filtering can''t be done, increase the decomposition level']); end;
    if(n>1)	
    	keepcno = S(1,1) * S(1,2) + sum(S(2:n+1,1) .* S(2:n+1,2) * 3);
    else
        keepcno = S(1,1) * S(1,2) + sum( S(2,1).* S(2,2) * 3);
    end

    C(keepcno+1:end) = 0.0;
%    reconField = waverec2(C, S, waveletName);
    reconField = mywaverec2(C,S,lpf);

end
