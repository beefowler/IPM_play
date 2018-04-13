function [F] = divisionkernel(x,y, pars)

alpha = pars(1); 
phi = pars(2); 

Alph1 = alpha; 
Alph2 = alpha*round((1-phi)/phi); 


f1 = (1/beta(Alph1, Alph2))*((x/y)^(Alph1 - 1))*((1-(x/y))^(Alph2 - 1)); %first f function
f2 = (1/beta(Alph1, Alph2))*((1-(x./y))^(Alph1-1))*((x/y)^(Alph2 - 1)); %second f function
F = (f1 + f2) / y; 

if F < 0; F = 0; 
    

end
 