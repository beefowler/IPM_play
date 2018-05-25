foundpars = []; 

%for i = 1:20; 
    
[~, ~, ~, Sample] = Simulate(ny, rs, Truepars, Nsums); 

%use the same data and calculate negloglike for different parameters. 
%then create a surace. 

alph_var = 1:30; 
phi_var = .1:.05:.5 ; 
data = zeros(length(alph_var), length(phi_var)); 
for i = 1:(length(alph_var))
    for j = 1:(length(phi_var))
        [nll] = negloglike_simp(rs, [alph_var(i) phi_var(j)], Sample, ny); 
        data(i,j) = nll; 
    end 
end 


minimum = min(min(data)); 
[i j] = find(data == minimum); 
alpha = alph_var(i)
phi = phi_var(j)

foundpars = [foundpars; alpha, phi]


%end 