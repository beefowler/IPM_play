function [Nday, NdayProps, mu, Sample] = Simulate2(ny, rs, pars, Nobserved)
Nt0 = ny
Nday = ny' 
NdayProps = ny' ./ trapz(rs, ny)

count = 1; 
while count < 144; 
    Nt1 = TimeStep2(Nt0, rs, pars, count/6)
    
    if rem(count,6) == 0 %if we are at an hour
         Nday = [Nday Nt1']; %final output will be matrix of N at each hour
         NdayProps = [NdayProps (Nt1' ./ trapz(rs, Nt1))]; 
    end 
    
    Nt0 = Nt1; 
    count = count +1; 
    
end 

%get sample from NdayProps only if Nobserved has been given
if exist('Nobserved', 'var')
    Sample = mnrnd(sum(Nobserved)', NdayProps')'; %sample from multinomial with probabilites of NdayProps
else 
    Sample = NaN; %just in case user asks for Sample output without giving sample sizes
end 

%calculate growth rate based on starting and end size of population 
mu = log(sum(Nday(:, end))/sum(Nday(:,1))); 

end 
