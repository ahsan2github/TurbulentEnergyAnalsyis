function [varargout] = getBBox(cx, cy, a, b, ort, boundx, boundy)
    jvector_x = java.util.Vector;
    jvector_y = java.util.Vector;
    phi = [0,pi/2, pi, 2*pi*3/4.0];
    cosphi = cos(phi);
    sinphi = sin(phi);
    % convrt to radian
    theta = pi* ort/180;   
    R = [ cos(theta)   sin(theta)
         -sin(theta)   cos(theta)];
    % calculate the four end points of major and minor axes 
    xy = [a*cosphi; b*sinphi];
    % generate points on major axis
    xx = ceil(min(xy(1,:))):floor(max(xy(1,:)));
    % for each point on major axis calculate all possible
    % y locations that fall within ellipse

    for i = xx
        proot = [-sqrt((1-i*i/a/a)*b*b), sqrt((1-i*i/a/a)*b*b)];
        yy = [real((min(proot))) real((max(proot)))];
        xcopy = ones(1,length(yy)).*i;
        % rotate the ellipse that was generated around (0,0) point 
        % to match the oreintation given and translate the ellipse  to the
        % centroid location
        rtxy = R*[xcopy;yy];
        xcor = (rtxy(1,:)+cx);
        ycor = (rtxy(2,:)+cy);
        % establish bound check
        xout1 = find(xcor < 1);
        xout2 = find(xcor > boundx);
        xcor([xout1, xout2]) = [];
        ycor([xout1, xout2]) = [];
        yout1 = find(ycor < 1);
        yout2 = find(ycor > boundy);
        xcor([yout1, yout2]) = [];
        ycor([yout1, yout2]) = [];
        xcor = round(xcor);
        ycor = round(ycor);
        if ~isempty(xcor)            
            for ii = 1:length(xcor)
                jvector_x.addElement(xcor(ii));
                jvector_y.addElement(ycor(ii));
            end

        end        
        
    end
    vector_x = zeros([jvector_x.size 1]);
    vector_y = zeros([jvector_x.size 1]);
    for tt = 0:jvector_x.size()-1
        vector_x(tt+1) =  jvector_x.elementAt(tt);
        vector_y(tt+1) =  jvector_y.elementAt(tt);
    end  
    max_x = max(vector_x(:)); min_x = min(abs(vector_x(:)));
    max_y = max(vector_y(:)); min_y = min(abs(vector_y(:)));

    horz_dim = (max_x - min_x); 
    vert_dim = (max_y - min_y); 
    box_cntr_horzc = 0.5*(max_x + min_x);
    box_cntr_vertc = 0.5*(max_y + min_y);
    switch(nargout)
        case 1
            varargout{1}= min_x;
        case 2            
            varargout{1} = min_x;
            varargout{2} = min_y;
        case 3
            varargout{1} = min_x;
            varargout{2} = min_y;
            varargout{3} = horz_dim;
        case 4
            varargout{1} = min_x;
            varargout{2} = min_y;
            varargout{3} = horz_dim;  
            varargout{4} = vert_dim;
        case 5
            varargout{1} = min_x;
            varargout{2} = min_y;
            varargout{3} = horz_dim;  
            varargout{4} = vert_dim;            
            varargout{5} = box_cntr_horzc;              
        case 6
            varargout{1} = min_x;
            varargout{2} = min_y;
            varargout{3} = horz_dim;
            varargout{4} = vert_dim;           
            varargout{5} = box_cntr_horzc;
            varargout{6} = box_cntr_vertc;     
        case 7
            varargout{1} = min_x;
            varargout{2} = min_y;
            varargout{3} = horz_dim;  
            varargout{4} = vert_dim;              
            varargout{5} = box_cntr_horzc;     
            varargout{6} = box_cntr_vertc;     
            varargout{7} = vector_x;  
        case 8
            varargout{1} = min_x;
            varargout{2} = min_y;
            varargout{3} = horz_dim;  
            varargout{4} = vert_dim;              
            varargout{5} = box_cntr_horzc;     
            varargout{6} = box_cntr_vertc;     
            varargout{7} = vector_x;   
            varargout{8} = vector_y;             
    end
end
