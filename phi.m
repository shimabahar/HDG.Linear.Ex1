% Base function of degree i on element e
function out = phi(i, e, x)
global h

c = (e-1/2)*h; % center of the elemnt
out = ((2/h)*(x-c))^(i-1);
