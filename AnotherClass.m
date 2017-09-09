classdef AnotherClass
   properties
     n;
     n_in_m;
     tmn;
     sigma_tmn;
   end
   methods
       function obj = AnotherClass(m, nscale)
          if nargin == 2
             if isnumeric(m)
                obj(1:length(m)) = AnotherClass; 
                for i = 1:length(m)
                    obj(i).n = m(i):nscale;
                    dx = 62.5001; dy = dx;
                    obj(i).n_in_m = 2.^(nscale-obj(i).n).*(dx*dy)^(0.5);                    
                    obj(i).tmn = zeros([1 length(obj(i).n)]);
                    obj(i).sigma_tmn = zeros([1 length(obj(i).n)]);
                end
             else
                error('Two Numeric Paramters Required to Initialize Class')
             end
          end
       end
   end
end
