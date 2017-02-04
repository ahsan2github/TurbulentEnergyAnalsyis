function [varargout] = mywavedec2(dat,lpf,levels)
C = java.util.Vector;
S = zeros([levels+2 2]);
S(end,:) = size(dat,1);
counter = levels+1;
datatmp = dat;
CC = zeros([1 numel(datatmp)]);
for ii = 1:levels  
    [y, ~] = mdwt(datatmp,lpf,1);  
    S(counter, :) = [size(datatmp,1)/2, size(datatmp,1)/2];
    LL = y(1:S(counter,1),1:S(counter,2));
    LH = y(S(counter,1)+1:end,1:S(counter,2));
    HL = y(1:S(counter,1), S(counter,2)+1:end);
    HH = y(S(counter,1)+1:end,S(counter,2)+1:end);
    CC(1:S(counter,1)*S(counter,2)) = LL(:);
    CC(S(counter,1)*S(counter,2)+1:S(counter,1)*S(counter,2)*2) = LH(:);
    CC(S(counter,1)*S(counter,2)*2+1:S(counter,1)*S(counter,2)*3) = HL(:);
    CC(S(counter,1)*S(counter,2)*3+1:S(counter,1)*S(counter,2)*4) = HH(:);
    if(ii==levels)
        S(counter-1, :) = size(LL);
    end     
    counter = counter-1;
    datatmp = LL;
end
varargout{1}=CC;
varargout{2}=S;
