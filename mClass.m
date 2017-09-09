classdef mClass
   properties
      m;
      m_in_m;
      h;
      n_data;
      dx = 62.5001;
      dz = 23.8095;
      strucFound;
   end
   methods
       function obj = mClass(nscale)
          if nargin ~= 1
                error('"nscale" Paramters Required to Initialize Class')
          else            
                obj.m = 1:nscale;
                dy = obj.dx;
                obj.m_in_m = 2.^(nscale-obj.m).*(obj.dx*dy)^(0.5);
                obj.n_data = AnotherClass(obj.m, nscale);
                obj.strucFound = zeros([1 nscale]);                
          end
       end
   end
end
