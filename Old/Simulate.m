
function [Nday, NdayProps, mu, Sample] = Simulate(ny, rs, pars, Nobserved)

%initialize output, start with initial population
Nday = ny'; 
NdayProps = ny' ./ trapz(rs, ny);

%every count is ten minutes; 
count = 1; 

while count < 144 %go until 1 day has passed
    Nt1 = TimeStep(ny, rs, pars, count/6); %project forward ten minutes
    if rem(count,6) == 0 %if we are at an hour
         Nday = [Nday Nt1']; %final output will be matrix of N at each hour
         NdayProps = [NdayProps (Nt1' ./ trapz(rs, Nt1))]; %proportion of population in each bin
    end %if statement
    
    ny = Nt1; %make new size distribution Nt0 
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
mu = log(sum(Nday(:, end))./sum(Nday(:,1))); 

end 
