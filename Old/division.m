%integrand of division portion of IPM 
%ouput is contribution of division to n(x,t+1) given x, y, ny = n(y,t), and parameters
%x should be a single value 
function [nx] = division(ydata, x, pars, edges)
%pars: [a b] s.t. a>= 0, b>0
    a = pars(1);
    b = pars(2);

realsize = 2.^(ydata(:,1)); 
    s = (realsize.^a)./(realsize.^a + b^a); %function for probability of dividing
    
%preset nx 
    nx = zeros(length(ydata),1); 
%ydata = [y ny] at time t, for any nuber of y values 
    ybin = discretize(x+1, edges); 
    if ~isnan(ybin) % sometimes NaN appears when y is too big
        nx(ybin) = 2 * s(ybin) * ydata(ybin, 2); % contribution of division to n(x,t+1)
    end
    
end %division function



