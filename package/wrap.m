function y=wrap(x)
%        y=wrap(x)
%
% wrap operator y=mod(x+pi,2*pi)-pi
%
%  Authors: V. Katkovnik. Januauy,  2007.
y= mod(x+pi,2*pi)-pi;

% Z2= cos(x); Z1=sin(x);
% y=angle(Z2+i*Z1);
% for s=1:length(x)
% if y(s)==pi
%     y(s)=-pi;
% end
% end


