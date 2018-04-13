
%calculates negative log likelihood of a set of parameters for IPM model
%given a day of observations 

%basically just does what Simulate does and then calculates negloglike 
%Inputs: all the inputs of Simulate plus observed 

%observed = counts of cells observed in each size bin for each hour of day 


function [nll] = negloglike(Realsize, pars, Einterp, Nobserved); 

%% Run Projection / Simulation  

%initialize output, start with initial population
Nt0 = Nobserved(:,1)'; 
Nday = Nt0'; 
NdayProps = Nt0' ./ sum(Nt0);

%every count is ten minutes; 
count = 1; 

while count < 144 %go until 1 day has passed
    Nt1 = TimeStep(Realsize, Nt0, pars, Einterp(count)); %project forward ten minutes
    Nt1(find(Nt1 == Inf)) = max(Nt1(find(Nt1 ~= Inf))); %Inf causes problems

    if rem(count,6) == 0 %if we are at an hour
         Nday = [Nday Nt1']; %final output will be matrix of N at each hour
         NdayProps = [NdayProps (Nt1' ./ sum(Nt1))]; %proportion of population in each bin 
    end %if statement
    
    Nt0 = Nt1; %make new size distribution Nt0 
    count = count+1;
    
end %while loop

%% Calculate NegLogLike  
logL = sum(sum(log(Nobserved .* NdayProps)));  % log(prob*observed) summed for every bin and hour

nll = -logL; %negative log likelihood to be minimized
 

if isnan(nll)
    keyboard; 
end

end 