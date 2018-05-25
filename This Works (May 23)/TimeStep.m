%rs = realsize of middle cell in each bin 
%ny = counts of cells in each bin 
%pars = 2 parameters defined below 
%t = time of day after dawn in hours 

function [nx] = TimeStep(ny, rs, pars, t) 

alpha = pars(1); 
phi = pars(2); 

%s = (rs+2).^10 ./ ((rs+2).^10 + 1); %just choosing these parameters for now 
%s = (rs).^4 ./ ((rs).^4 + 35); %just choosing these parameters for now
s = (rs).^6 ./ ((rs).^6 + 600); %just choosing these parameters for now



%Division 
Alph1 = alpha; 
Alph2 = alpha*(1-phi)/phi; 

xovery = rs' * (1./rs); %get values of x/y, each collumn a y, each row an x 
xovery(find(xovery >= 1)) = 0;  %dirichlet funciton only applies to values between 0 and 1, ie when x is less than y 
f1 = (1/beta(Alph1, Alph2)).*(xovery.^(Alph1 - 1)).*((1-xovery).^(Alph2 - 1)); %first f function
f2 = (1/beta(Alph1, Alph2)).*((1-xovery).^(Alph1-1)).*((xovery).^(Alph2 - 1)); %second f function
F = (f1 + f2) ./ rs; 

Div = F .* ny .* s; 

%Growth (AND STASIS)

xminusy = ((rs)' * ones(1, length(rs))- (ones(length(rs),1)*(rs))); %Matrix with values of x - y for collumns y and rows x
beeta = .05*sin(pi*t/24); %trying to mimic sunlight patterns 

Hkern = heaviside(xminusy) .* 1./beeta .* exp(-abs(xminusy) ./ (beeta)); 
Grow = (1-s) .* ny .* Hkern; 

keyboard
nx = trapz(rs, Div') + trapz(rs, Grow') ; 


end
