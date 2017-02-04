function [LL] = mywaverec2(datatmp,S,lpf)
    LL = datatmp(1:S(1,1)*S(1,2));
    for ii = 2:size(S,1)-1 
        LL = reshape(LL, [S(ii,1) S(ii,2)]);
        LH = datatmp(S(ii,1)*S(ii,2)+1:S(ii,1)*S(ii,2)*2);
        LH = reshape(LH, [S(ii,1) S(ii,2)]);
        HL = datatmp(S(ii,1)*S(ii,2)*2+1:S(ii,1)*S(ii,2)*3);
        HL = reshape(HL, [S(ii,1) S(ii,2)]);
        HH = datatmp(S(ii,1)*S(ii,2)*3+1:S(ii,1)*S(ii,2)*4);
        HH = reshape(HH, [S(ii,1) S(ii,2)]);
        LL = midwt([LL,HL;LH,HH],lpf,1);
    end
end

