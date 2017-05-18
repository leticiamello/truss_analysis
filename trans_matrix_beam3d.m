function T = trans_matrix_beam3d(varargin)
%trans_matrix_beam3d   Transformation matrix of a beam in space
%
%   Function T = trans_matrix_beam3d(x0,y0,z0,x1,y1,z1,alfa)
% 
%   Function gives the Global to Local coordinate transformation
%   Matrix T. x0,y0,z0,x1,y1,z1 are the coordinate in Global system
%   of beam joints 0 and 1. Angles alfa [rad], defined if the beam 
%   is rotate around its axis. If alfa is ommited, function sets
%   alfa to 0.
%
%   Example
%       T = trans_matrix_beam3d(0,0,0,1,0,0) or
%       T = trans_matrix_beam3d(0,0,0,1,0,0,0)
%   give the identity matrix.
% 

%   version 1.0. Developed by Paulo J. Paupitz Goncalves.



if length(varargin) < 6
    error('Not enough input arguments')
elseif length(varargin) > 7    
    error('Too many input arguments')
else
    x0 = varargin{1}; x1 = varargin{4};
    y0 = varargin{2}; y1 = varargin{5};
    z0 = varargin{3}; z1 = varargin{6};
    
    if length(varargin) == 6
        alfa = 0;
    else
        alfa = varargin{7};
    end
    L = sqrt((x1-x0)^2 + (y1-y0)^2 + (z1-z0)^2);    % Beam length
    
    Cx = (x1-x0)/L;     %%% direction cosine 
    Cy = (y1-y0)/L;     %%% direction cosine 
    Cz = (z1-z0)/L;     %%% direction cosine 
    ca = cos(alfa);
    sa = sin(alfa);
    
    if (Cx == 0) && (Cz == 0)
        % When the bar is vertical
        T = [ 0              Cy    0
             -Cy*ca   0     sa
              Cy*sa   0     ca];
    else
        % all other cases
        Cxz = sqrt(Cx^2 + Cz^2);
        T = [  Cx                                    Cy                Cz
             -(Cx*Cy*ca+Cz*sa)/Cxz     Cxz*ca   -(Cy*Cz*ca+Cx*sa)/Cxz
              (Cx*Cy*sa-Cz*ca)/Cxz    -Cxz*sa    (Cy*Cz*sa+Cx*ca)/Cxz];
    end
end
 