
function [nll] = negloglike_simp(Realsize, pars, Nobserved); 

%% Run Projection / Simulation  

%initialize output, start with initial population
ny = Nobserved(:,1)' ; 
Nday = Nobserved(:,1); 
NdayProps = Nobserved(:,1) ./ trapz(Realsize, ny);

%every count is ten minutes; 
count = 1; 

while count < 144 %go until 1 day has passed
    Nt1 = TimeStep(ny, Realsize, pars, count/6); %project forward ten minutes
    if rem(count,6) == 0 %if we are at an hour
         Nday = [Nday Nt1']; %final output will be matrix of N at each hour
         NdayProps = [NdayProps (Nt1' ./ trapz(Realsize, Nt1))]; %proportion of population in each bin 
    end %if statement
    
    ny = Nt1; %make new size distribution Nt0 
    count = count+1;
    
end %while loop

NdayProbs = NdayProps ./ sum(NdayProps); 

% %% Calculate NegLogLike 
% NdayProbs = NdayProbs*100; 
% temp = NdayProbs .* Nobserved; 
% ind = find(temp ~=0); 
% temp = log(temp(ind)); 
% nll = -sum(temp); 
% 
% keyboard 


temp = Nobserved .* Nday; 
temp = temp(find(temp~=0)); 
temp = log(temp); 
logL = sum(temp);  % log(prob*observed) summed for every bin and hour

nll = -logL; %negative log likelihood to be minimized



if isnan(nll)
    keyboard; 
end
 


end 