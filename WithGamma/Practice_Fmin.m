


lb=[3 .2 .8 .8 1 5 2]; %negative lower bounds
ub=[5 .5 1.1 2 2.2 7 4]; %upper bounds 

x0 = zeros(1,7); 
for i = 1:7
    x0(i) = lb(i) + (ub(i) - lb(i))*rand(1); 
end 


opts=optimset('TolX', 1e-8, 'Algorithm','interior-point','MaxIter', 3000,'MaxFunEvals',10000);
[optpars,optnll]=fmincon(@(pars) negloglike(Realsize, pars, Einterp, Nobserved),x0,[],[],[],[],lb,ub,[], opts)


%[Expectation, expectedProps, my, SimSample] = Simulate(Realsize, Nobserved(:,1), optpars, Einterp, Nobserved); 

