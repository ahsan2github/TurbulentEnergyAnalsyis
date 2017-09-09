classdef piClass_
  properties
		n;
		n_in_m;
		h;
		m_data;
		dx; dy; dz;
		pi_energy; pi_energy_sigma;
        strucFound;
	end
	methods
		function obj = piClass_(nscale, dx,dy,dz)
			if nargin ~=4
				error(['Four parameters required']);
			else
				obj.n = 2:nscale;
				obj.dx = dx;
				obj.dy = dy;
				obj.dz = dz;
				pi_energy = zeros([1 length(obj.n)]);
				pi_energy_sigma = zeros([1 length(obj.n)]);
				obj.n_in_m = 2.^(nscale-obj.n).*(dx*dy)^(0.5);
				obj.m_data = M4NClass_(obj.n, nscale, dx, dy);
                obj.strucFound = zeros([1 length(obj.n)]);
			end
		end
	end
end
