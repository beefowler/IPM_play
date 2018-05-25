
function [nll] = negloglike_simp(Realsize, pars, Nobserved, Nt0); 

%% Run Projection / Simulation  

%later, initialize output, start with initial population
%ny = Nobserved(:,1)' ; 
%Nday = Nobserved(:,1); 
%NdayProps = Nobserved(:,1) ./ trapz(Realsize, ny);

%for now let's just have it preset
ny = Nt0;
Nday = Nt0'; 
NdayProps = Nday ./ trapz(Realsize, ny);


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
%Let's cut down size of Nobserved to scale probabilities larger
%factor = sum(sum(Nobserved)) / 1000; 
%Observed1000 = Nobserved / factor; 

%try percentages of hourly populations instead of portions of all
%observations. 
Observed1000 = Nobserved ./ sum(Nobserved) * 100; 

temp = Observed1000 .* log(NdayProbs); 
nll = -sum(sum(temp)); %negative log likelihood to be minimized


if isnan(nll)
    keyboard; 
end
 

end 