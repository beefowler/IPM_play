
%This function performs IPM projection for one time step, given size
%distribution at t0, parameters, and precise light data. 

%Inputs: 
%Realsize = the center points for the size bins for cell volume 
%Nt0 = the counts of cells in each size bin, defined by Realsize, at t0
%pars = parameter values 
     % pars(1) = Alph, dispersal parameter for dirichlet. dist. in division
     % pars(2) = Phi, proportion of mother cell that goes to smaller of two
     % daughter cells, 0 < Phi <= 1/2
     % pars(3) = a,  for s probability of division, a >= 1
     % pars(4) = b,  for s probability of division, b>0
     % pars(5) = phi, for growth kernel 
     % pars(6) = psi, for growth kernel 
     % pars(7) = ro,  for growth kernel 
%E = value of incidient radiation for the time of interest. SCALAR
     
%Output: 
%Nt1 = expected value of count of cells in each size bin defined 
%by Realsize, after cells in Nt0 have divided or grown for one time step 

function [Nt1] = TimeStep(Realsize, Nt0, pars, E)

%label parameter values 
Alph = pars(1);
Phi = pars(2);
%a = pars(3);
%b = pars(4);
phi = pars(3);
psi = 3; %pars(6);
ro = 3; %pars(7); 

%% Division
%s = (Realsize.^a)./(Realsize.^a + b^a); %probability of dividing for every value of y in domain
s = (Realsize).^10 ./ ((Realsize).^10 + 3); 

%more parameters based on above 
Alph1 = Alph; 
Alph2 = Alph*(1-Phi)/Phi; 

xdomain = Realsize'* (1./Realsize) ; %get values of x/y 
xdomain(find(xdomain >= 1)) = 0;  %dirichlet funciton only applies to values between 0 and 1, ie when x is less than y 

%See Sweavetest.pdf for details on Division Kernel 
f1 = (1/beta(Alph1, Alph2)).*(xdomain.^(Alph1 - 1)).*((1-xdomain).^(Alph2 - 1)); %first f function
f2 = (1/beta(Alph1, Alph2)).*((1-xdomain).^(Alph1-1)).*((xdomain).^(Alph2 - 1)); %second f function
F = (f1 + f2) ./ Realsize; %Fecundity kernel, don't know why imaginaries are appearing 

Div = F .* s .* Nt0; %Division integrand

    
%% GROWTH

Betaval = 1; 
Beta = Betaval.*ones(length(Realsize),length(Realsize)); %2 is an arbitrary, but reasonable choice. Could change later. 

alphalpha = Betaval.*Realsize' + (Betaval*ro*E^phi)/(E^phi + psi^phi); %helpful term for gamma dist
alphalpha = alphalpha*ones(1,length(Realsize)); 
if alphalpha <= 0 %alphapha needs to be positive
        disp('err: alpha < 0') 
        alphalpha(find(alphalpha<0)) =0; 
end 

domain = Realsize' * ones(1,length(Realsize)) - (ones(length(Realsize), 1) * Realsize); 
K = gampdf(domain, Beta, alphalpha); %growth kernel is gamma distribution 
Grow = K .* (1-s) .* Nt0; %Growth integrand Kernel times probability of not splitting times number of cells

 
%% INTEGRATE 
Nt1 = trapz(Realsize, Div') + trapz(Realsize, Grow'); %final output will be n(x,t+1) over domain Realsize
 

end
