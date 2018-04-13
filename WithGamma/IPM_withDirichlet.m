%What we have: 
%day, 
%DielStart 
%DSBeads
%DScellresults
%DSMerged
%Eall
%volbins
%eukvols 
%Einterp 

%Create bins for IPM
Realsize= 2.^[-10:(2/40):10]; 
d = diff(Realsize)/2; 
edges = [Realsize(1) - d(1), Realsize(1:end-1)+d, Realsize(end) + d(end)]; %Convert enters to edges 
clear d 

%Get initial volumes 
y = sort(eukvols{1});
ny = histcounts(y, edges);

N_allday = [ny']; 
for hr = 2:length(eukvols);
    tempvols = histcounts(sort(eukvols{hr}), edges); 
    N_allday = [N_allday tempvols']; %store volumes observed for each hour 
end
clear tempvols 

%Parameters 
Alph = 4; 
Phi = 1/3; 
a = 1; 
b = 1; 
phi = 2; 
psi = 6; 
ro = 2; 



%% DIVISION 

s = (Realsize.^a)./(Realsize.^a + b^a); %probability of dividing for every value of y in domain, don't need to keep repeating

%Parameters 
%    Alph  and Phi, preset for now 

Alph1 = Alph; 
Alph2 = Alph*round((1-Phi)/Phi); %will always be an integer, just get rid of extra decimals 

xdomain = Realsize'* (1./Realsize) ; %get values of x/y 

%See Sweavetest.pdf for details on Division Kernel 
f1 = 1/beta(Alph1, Alph2).*xdomain.^(Alph1 - 1).*(1-xdomain).^(Alph2 - 1); %first f function
f2 = 1/beta(Alph1, Alph2).*(1-xdomain).^(Alph1-1).*(xdomain).^(Alph2 - 1); %second f function
F = (f1 + f2) .* (1./Realsize); %Fecundity kernel, don't know why imaginaries are appearing 
Div = F .* s .* ny; %Division integrand
Div(find(Div<0)) = 0; %remove negative values 

%% GROWTH

Betaval = 2; 
Beta = Betaval.*ones(401,401); %2 is an arbitrary, but reasonable choice. Could change later. 

E = 200; 

alphalpha = Betaval.*Realsize' + (Betaval*ro*E^phi)/(E^phi + psi^phi) + 1; %helpful term for gamma dist 
alphalpha = alphalpha*ones(1,401); 
if alphalpha <= 0 %alphapha needs to be positive
        disp('err: alpha < 0') 
        alphalpha(find(alphalpha<0)) =0; 
end 

domainlg = (log(Realsize)' * ones(1, length(Realsize))- (ones(length(Realsize),1)*log(Realsize))); %Matrix with values of x - y for collumns y and rows x
domain = ((Realsize)' * ones(1, length(Realsize))- (ones(length(Realsize),1)*(Realsize))); %Matrix with values of x - y for collumns y and rows x
%%%STRUGLING WITH WHETHER THIS DOMAIN SHOULD BE LOG OR NOT
%ALSO WHAT IS THE SCALE OF K
K = gampdf(domain, Beta, alphalpha); %growth kernel is gamma distribution 
Grow = K .* (1-s) .* ny; %Growth integrand Kernel times probability of not splitting times


%% INTEGRATE 

nx = trapz(Realsize, Div') + trapz(Realsize, Grow'); %final output will be n(x,t+1) over domain Realsize

