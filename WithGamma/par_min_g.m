%Outputs = []; 

for i = 1:10; 

[~, ~, ~, Nobserved] = Simulate_g(Realsize, Nt0, [15 .2 5], Einterp, Nsums); 

lb = [2 .005  .005]; 
ub = [20 .5   10]; 

x0 = zeros(1,3); 
for i = 1:3
    x0(i) = lb(i) + (ub(i) - lb(i))*rand(1); 
end 


opts=optimset('TolX', 1e-8, 'Algorithm','interior-point','MaxIter', 3000,'MaxFunEvals',10000);
[optpars,optmse]=fmincon(@(pars) mean_square_err_g(Realsize, pars, Nobserved),x0,[],[],[],[],lb,ub,[], opts)


Outputs = [Outputs; optpars optmse]; 


end
