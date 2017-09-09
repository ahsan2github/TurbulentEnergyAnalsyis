classdef M4NClass
   properties
     m;
     m_in_m;
     tmn
   end
   methods
       function obj = M4NClass(n, nscale, dx, dy)
            if nargin >0
       			obj(length(n)) = M4NClass; 
                for i = 1:length(n)
                    obj(i).m = 1:n(i);
                    obj(i).tmn = zeros([1 length(1:n(i))]);
                    obj(i).m_in_m = 2.^(nscale-obj(i).m).*(dx*dy)^(0.5);                    
                end
            end
       end
   end
end
