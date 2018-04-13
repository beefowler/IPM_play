Outputs = []; 

for i = 1:20; 

[~, ~, ~, Nobserved] = Simulate(ny, rs, [4 .4], Nsums); 

lb = [2 .005]; 
ub = [20 .5]; 

x0 = zeros(1,2); 
for i = 1:2
    x0(i) = lb(i) + (ub(i) - lb(i))*rand(1); 
end 


opts=optimset('TolX', 1e-8, 'Algorithm','interior-point','MaxIter', 3000,'MaxFunEvals',10000);
[optpars,optmse]=fmincon(@(pars) mean_square_err(Realsize, pars, Nobserved),x0,[],[],[],[],lb,ub,[], opts)


Outputs = [Outputs; optpars]; 


end
