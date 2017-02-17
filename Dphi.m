% Derivative of base function of degree i on element e
function out = Dphi(i, e, x)
global h

if i == 1
    out = 0;
else
    c = (e-1/2)*h; % center of the elemnt
    out = (2/h)^(i-1) * (i-1)*(x-c)^(i-2);
end
