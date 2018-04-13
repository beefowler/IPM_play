%This is a first attempt at the Itegral Projection Model for PicoEukaryotes 
%Input all observations of class of interest (eukaryotes) at a given time
    %probably want these to already have been converted into volumes
%Sunlight data for the given time 
%parameters, later will be fed into function and optimized, 
    %hardcoded for now

%Euks_at_t1 can be a single cell X by 8, where X is the number of
%individual plankters observed. 8 collumns are our usual:  
 %    {'rec number'}    {'PE'}    {'FLS'}    {'CHL'}    {'SSC'}

 %    {'CHLpk1'}    {'CHLpk2'}    {'Class'}
% though in fact only SSC is necessary 


%mesh points
%want y and x on log 2 scale. 
%cells classified in meshphts(i) = k are of size <= 2^k 
%cells in bin are ~2 * size of those in previous bin 
edges = -10.02:(20/500):10;
meshpts = -10:(20/500):10;
meshpts(end) = []; 
realsize = 2.^meshpts; 

%Vols_at_t1 = cytosub_SSC2vol(Euks_at_t1(:,5)./ beadmatch); %create vector of volumes within the hour using cytosub function
y = log2(sort(Vols_at_t1)); %sort and log scale observed volumes 

ny = histcounts(y, edges); %counts of y size cells in each bin defined by meshbts 
ydata = [meshpts; ny]'; %n(y, t) over domain meshpts 

%preset parameters for now, later we will optimize 
pars = [1 1 2 1 -3]; %parameters a and b for division, psi phi and ro for growth 
a = pars(1);
b = pars(2);
phi = pars(3); 
psi = pars(4); 
ro = pars(5); 

s = (realsize.^a)./(realsize.^a + b^a); %probability of dividing for every value in domain, don't need to keep repeating

function [xdata] = IMPtimestep (ydata, E) 
xdata = zeros(1,2); 
    for x = meshpts
        pt1 = division(ydata, x, pars(1:2), edges) ; 
        pt2 = growth(ydata, x, E, pars, edges); 
        nx = sum(trapz(meshpts, [pt1 pt2])); %trapz integrates over each collumn separately, then we sum 
        xdata = [xdata; x nx]; %final output will be n(x,t+1) over domain meshpts  
    end %for loop x 
xdata(1,:) = []; %clear initial zeros 
end %time step function 

