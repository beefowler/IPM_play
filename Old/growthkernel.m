function [K] = growthkernel(x, y, t); 

beeta = sin(pi*t/24); 

K = heaviside((x - y)) .* 1/beeta .* exp(-abs(x - y) ./ beeta); 


end
