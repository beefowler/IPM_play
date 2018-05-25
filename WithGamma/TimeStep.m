%rs = realsize of middle cell in each bin 
%ny = counts of cells in each bin 
%pars = 2 parameters defined below 
%t = time of day after dawn in hours 

function [nx] = TimeStep(ny, rs, pars, t) 

alpha = pars(1); 
phi = pars(2); 
G = pars(3); 
b = pars(4); 

%s = (rs+2).^10 ./ ((rs+2).^10 + 1); %just choosing these parameters for now 
%s = (rs).^4 ./ ((rs).^4 + 35); %just choosing these parameters for now
s = (rs).^6 ./ ((rs).^6 + 6000); %just choosing these parameters for now



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
xminusy(find(xminusy<0))=0; 
xminusyscaled = xminusy ./ max(xminusy); %convert to domain (0,1)
xminusyscaled(find(isnan(xminusyscaled))) = 0; %remove NaNs

a = G*sin(pi*t/24)+1; %trying to mimic sunlight patterns 
%b = 5; 

%Kumaraswamy distribution, scaled to integrate to 1 over rs. 
Hkern = heaviside(xminusyscaled) .*a.*b.*xminusyscaled.^(a - 1).*(1 - xminusyscaled.^a).^(b - 1) ./ max(xminusy); 
Hkern(:,end) = zeros(1, length(rs)); %override last row to keep from NaNs. 
Hkern(end, end) = 1; 

Grow = (1-s) .* ny .* Hkern; 

keyboard 

nx = trapz(rs, Div') + trapz(rs, Grow') ; 


end
