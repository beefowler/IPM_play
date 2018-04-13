%Given initial cell sizes, parameters and Einterp for 1 day, this simulation projects
%the IPM model through 24 hours of TimeStep function. 

%Inputs: 
%Realsize = the center points for the size bins for cell volume 
%Nt0 = the counts of cells in each size bin, defined by Realsize, at t0
%pars = parameter values 
     % pars(1) = Alph, dispersal parameter for dirichlet. dist. in division
     % pars(2) = Phi, proportion of mother cell that goes to smaller of two
     % daughter cells, 0 < Phi <= 1/2
     % pars(3) = a,  shape parameter for s probability of division, a >= 1
     % pars(4) = b,  position parameter for s probability of division, b>0
     % pars(5) = phi, shape parameter for growth kernel, how average amount
     % of growth varies with light 
     % pars(6) = psi, position parameter for growth kernel, how much is an
     % average amount of growth in one time step
     % pars(7) = ro,  max level for growth kernel, when E is super
     % saturated, how much can a cell grow 
%Einterp = 1 by 145 vector of light values for 1 day, every 10 minutes
%Nobserved = optional input, used for getting sample size of simulated data 

%Output: 
%Nday = Matrix of expected value of counts of cells in each size bin (rows) 
%defined by Realsize, for each hour of day (collumn)
%NdayProps = matrix of expected values divided by expected total population
%size so each size bin has a proportion instead of a count
%mu = daily growth rate, according to simulation 
%Sample = Matrix of actual counts of cells in each size bin sampled from
%NdayProps with sample size according to Nobserved 

function [Nday, NdayProps, mu, Sample] = Simulate(Realsize, ny, pars, Einterp, Nobserved)

%initialize output, start with initial population
Nt0 = ny; 
Nday = Nt0'; 
NdayProps = Nt0' ./ trapz(Realsize, Nt0);

%every count is ten minutes; 
count = 1; 

while count < 144 %go until 1 day has passed
    Nt1 = TimeStep_g(Realsize, Nt0, pars, Einterp(count)); %project forward ten minutes
    if rem(count,6) == 0 %if we are at an hour
         Nday = [Nday Nt1']; %final output will be matrix of N at each hour
         NdayProps = [NdayProps (Nt1' ./ trapz(Realsize, Nt1))]; %proportion of population in each bin 
    end %if statement
    
    Nt0 = Nt1; %make new size distribution Nt0 
    count = count+1;
    
end %while loop

%get sample from NdayProps only if Nobserved has been given
if exist('Nobserved', 'var')
    NdayProbs = NdayProps ./ sum(NdayProps); 
    Sample = mnrnd(sum(Nobserved)', NdayProbs')'; %sample from multinomial with probabilites of NdayProps
else 
    Sample = NaN; %just in case user asks for Sample output without giving sample sizes
end 

%calculate growth rate based on starting and end size of population 
mu = log(sum(Nday(:, end))/sum(Nday(:,1))); 


end 

