%Outputs = []; 

for i = 1:10; 

[~, ~, ~, Nobserved] = Simulate_g(Nt0, Realsize, Truepars, Nobserved); 

lb = [2 .005  .005   10]; 
ub = [20 .5   10   50]; 

x0 = zeros(1,3); 
for i = 1:4
    x0(i) = lb(i) + (ub(i) - lb(i))*rand(1); 
end 


opts=optimset('TolX', 1e-8, 'Algorithm','interior-point','MaxIter', 3000,'MaxFunEvals',10000);
[optpars,optnll]=fmincon(@(pars) negloglike_simp(Realsize, pars, Nobserved, Nt0),x0,[],[],[],[],lb,ub,[], opts)


Outputs = [Outputs; optpars optnll]; 


end
