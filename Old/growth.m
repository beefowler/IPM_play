
%Integrand of growth portion of IPM
%output is contribution of growth to n(x,t+1) given x, y, ny = n(y,t), and parameters
function [nx] = growth(ydata, x, E, pars, edges)
%pars = [ a b phi, psi, ro] 
    a = pars(1);
    b = pars(2);
    phi = pars(3); 
    psi = pars(4); 
    ro = pars(5); 

    realsize = 2.^(ydata(:,1)); 
    s = (realsize.^a)./(realsize.^a + b^a); %function for probability of dividing
    
    beta = 2; %This is an arbitrary, but reasonable choice. Could change later. 
    
    alpha = beta*x - (beta*ro*E^phi)/(E^phi + psi^phi) + 1; %helpful term for gamma dist  
    if alpha <= 0 %alpha needs to be positive
        disp('err: alpha < 0') 
        %keyboard
    end 
    
    xbin = discretize(x,edges); 

    k = gampdf((x - ydata(:,1)), beta, alpha); %growth kernel is a gamma distribution
    
    nx = (1-s).* k .* ydata(:,2); %new x size cells from growth of y size cells
    
end %growth function 
