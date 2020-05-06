function [dr,badg] = drsnbrck(x)
%function dr = drsnbrck(x)
dr=zeros(2,1)
dr(1,1) = 2*(x(1)-1) - 8*105*x(1)*(x(2)-x(1)^2)^3;
dr(2,1) = 4*105*(x(2)-x(1)^2)^3;
badg=0