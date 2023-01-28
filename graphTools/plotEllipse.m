%% Plot ellipse
%  PLOTELLIPSE plots an ellipse given its parameters
%
%   M. A. Isa UoN, 2021
function [u2,v2,varargout]=plotEllipse(center,rad,theta,col,varargin)
    x0=center(1); y0=center(2); a=rad(1); b=rad(2);
    n=130; show_plot=true;
    if nargin>=5
        n=varargin{1};
    end
    if nargin>=6
        show_plot=varargin{2};
    end
    th=(0:2*pi/(n-1):2*pi)'; 
    u1=a*cos(th); v1=b*sin(th);
    u2=x0+u1*cos(theta)-v1*sin(theta);
    v2=y0+u1*sin(theta)+v1*cos(theta);
    
    if nargout>=3
        varargout{1}=th;
    end
    
    if show_plot
        plot(u2,v2,col, 'LineWidth',1.2);
        axis equal;

    end
    

end
