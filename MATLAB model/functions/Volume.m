function [V,dVdtheta]=Volume(Vdead,Vdisp,theta)
%% This function returns the value of volume and volume dervative at a specific crank angle
V = Vdead+Vdisp/2*(1-cos(theta));
dVdtheta = Vdisp/2*sin(theta);

end