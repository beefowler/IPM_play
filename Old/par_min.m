Outputs = []; 

for i = 1:20; 

[~, ~, ~, Sample] = Simulate(ny, rs, Truepars, Nsums); 

lb = [2 .005]; 
ub = [40 .5]; 

x0 = zeros(1,2); 
for i = 1:2
    x0(i) = lb(i) + (ub(i) - lb(i))*rand(1); 
end 


opts=optimset('TolX', 1e-8, 'Algorithm','interior-point','MaxIter', 3000,'MaxFunEvals',10000);
[optpars,optmse]=fmincon(@(pars) negloglike_simp(Realsize, pars, Sample, ny),x0,[],[],[],[],lb,ub,[], opts)


Outputs = [Outputs; optpars optmse]; 


end
