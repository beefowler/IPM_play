%integrand of division portion of IPM, with dirichlet distribution 
%ouput is contribution of division to n(x,t+1) given y, ny = n(y,t), and parameters
%let x be equal to all meshpts
function [nx] = division2(ydata, x, pars)
%pars: [a b] s.t. a>= 0, b>0
    alpha1 = pars(1)
    phi = pars(2)
    
realsize = 2.^(ydata(:,1)); 
    s = (realsize.^a)./(realsize.^a + b^a); %function for probability of dividing
    
%preset Fecundity matrix 
    Div = zeros(length(realsize)); 
    
    for j = 1:length(ydata)
        f1 = (1/(beta(alpha1, ((alpha1)*(1-phi)/phi)))) .*((x/y(j)).^(alpha1 - 1)).* ((1-(x/y(j))).^((alpha1*(1-phi)/(phi))-1));
        f2 = (1/(beta(alpha1, ((alpha1)*(1-phi)/phi)))) .*((1-x/y(j)).^(alpha1 - 1)).* ((x/y(j)).^((alpha1*(1-phi)/(phi))-1));
        Div(:,j) = 2 .* s .* (1/y(j) * [f1 + f2]) .* ny(j); 
    end
    Div = real(Div); 
    Div(isnan(Div)) = 0;

    pt1 = trapz(meshpts, [Div']); 
    
    
end %division function
