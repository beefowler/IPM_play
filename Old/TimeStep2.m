function [nx] = TimeStep2(ny, rs, pars, t)

Alph1 = pars(1); 
phi = pars(2); 
Alph2 = Alph1*(1-phi)/phi; 

nx = []; 
for x = rs 
      f1 = (1/beta(Alph1, Alph2)).*(x./rs).^(Alph1 - 1).*((1-(x./rs)).^(Alph2 - 1)); %first f function
      f2 = (1/beta(Alph1, Alph2)).*(1-(x./rs)).^(Alph1-1).*(x./rs).^(Alph2 - 1); %second f function
      F = (f1 + f2)./ rs; 
      F(find(F<0)) = 0; 
      Div = F .* ny * (1/2); 
    
      beeta = sin(pi*t/24); 
      Hkern = heaviside(x - rs) .* 1/beeta .* exp(-(x - rs) ./  beeta); 
      Hkern(find(isnan(Hkern))) = 0; 
      Grow = (1/2)* ny .* Hkern; 
      
      nx1 = sum(trapz(rs, [Div' Grow'])); 
      nx = [nx nx1]; 
end 


end 
