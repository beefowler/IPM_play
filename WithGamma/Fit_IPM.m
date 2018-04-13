%Setup MultiStart


function [Results, modelfits, Starts, Expectation] = Fit_IPM(Realsize, NObserved, Einterp) 

%Results = [pars negloglike mu exitflag length(Fits)] 

ms=MultiStart('Display','off','TolX',1e-5,'UseParallel','always','StartPointsToRun','bounds-ineqs');
opts=optimset('Display','off','TolX',1e-8,'Algorithm','interior-point','UseParallel','always','MaxIter', 3000,'MaxFunEvals',10000);

%% Set Up Paramter Bounds and Starting Values 

%Define parameter bounds 
lb=-[1e-4 1e-4 1 1e-4 1 1e-4 1e-4]; %negative lower bounds
ub=[10 .5 10 5 10 4 4]; %upper bounds 

%set parameter bounds to be interpreted by fmincon
    a1=-1*eye(7); 
    a2=eye(7);
    A=zeros(14,7);
    A(1:2:13,:)=a1;
    A(2:2:14,:)=a2;

    B=zeros(14,1);
    B(1:2:13)=lb;
    B(2:2:14)=ub;

%Create random initial parameter values within bounds 
x0 = zeros(1,7);
num_starts = 40; %choose the number of starting values we want, number of trials of fmincon
x2 = zeros(num_starts, 7); %matrix to create tpoints
for i = 1:7
    x0(i) = -lb(i) + (ub(i) + lb(i))*rand(1);  %remember lb has negative lower bounds. so this is lower bound + rand * difference 
    x2(:,i) = -lb(i) + (ub(i) + lb(i))*rand(num_starts,1);
end
clear i
tpoints = CustomStartPointSet(x2); 

%% Run Multistart 
    
problem = createOptimProblem('fmincon','x0',x0,'objective',@(pars) negloglike(Realsize, pars, Einterp, NObserved),'Aineq',A,'bineq',B,'options',opts);
[xmin,fmin,exitflag,~,soln] = run(ms,problem,tpoints);
    
    
%Investigate Solution 
    %open up the soln structure:
    temp=zeros(num_starts,10);
    start_points=zeros(num_starts,7);
    temp=zeros(num_starts,10);
    c=1;
    
    for j=1:length(soln)
        %check to see if all start points led to an individual solution or
        %not (MultiSTart will only return unique solutions)
        g=cell2mat(soln(j).X0);
        if length(g)==7 %only one start_point led to that solution
            start_points(c,:)=g;
            temp(c,1:7)=soln(j).X;
            temp(c,8)=soln(j).Fval;
            %temp(c,9)=growth_rate(Einterp,volbins,dayCOUNTS,temp(c,1:13),hr1,hr2);
            temp(c,10)=soln(j).Exitflag;
            c=c+1;
        else
            num=length(g)/7;
            start_points(c:c+num-1,:)=squeeze(reshape(g',1,7,num))';
            temp(c:c+num-1,1:7)=repmat(soln(j).X,num,1);
            temp(c:c+num-1,8)=repmat(soln(j).Fval,num,1);
            %temp(c:c+num-1,9)=repmat(growth_rate(Einterp,volbins,dayCOUNTS,temp(c,1:13),hr1,hr2),num,1);
            temp(c:c+num-1,10)=repmat(soln(j).Exitflag,num,1);
            c=c+num;
        end
    end
    %just in case have rows left as zeros
    qq=find(temp(:,1)~=0);
    temp=temp(qq,:);


    modelfits=[temp(1:10)];
    start_points=start_points(qq,:);
    Starts=start_points;

    %let's now ask, in the first batch run, did the solver "converge"?
    [sortlogL ii]=sort(modelfits(:,8));
    if abs(sortlogL(min(5,size(sortlogL,1)))-sortlogL(1)) < icsTol  %Did likelihood converge within tolerance
        flag1 = 0; %Yes it did converge
    else
        disp(num2str(sortlogL(1:min(5,size(sortlogL,1)))))
        flag1 = 1; %no it didn't 
    end;

       %did parameters converge within a tolerance 
    partol=max(modelfits(ii(1:min(5,size(sortlogL,1))),1:7))-min(modelfits(ii(1:min(5,size(sortlogL,1))),1:7));
    if sum(abs(partol) < tolvec)==7 || sum((abs(partol./modelfits(ii(1),1:7)) < 0.05))==7 %either the modelfits are within an absolute tolerance or within a relative tolerance
        flag2 = 0; %yes they did converge 
    else
        flag2 = 1; %no they didn't converge 
    end

    disp(['flag1 = ' num2str(flag1) ' flag2=' num2str(flag2)])

    k=1; %batch number
%     while ((flag1 || flag2) | size(modelfits,1) <= 20) && k <= 5 %run MultiStart 4 more times D: 
% 
%         disp(['k: ' num2str(k)])
%         k=k+1;
% 
%         %Create random initial parameter values within bounds 
%         x0 = zeros(1,7);
%         num_starts = 40; %choose the number of starting values we want, number of trials of fmincon
%         x2 = zeros(num_starts, 7); %matrix to create tpoints
%         for i = 1:7
%           x0(i) = -lb(i) + (ub(i) + lb(i))*rand(1);  %remember lb has negative lower bounds. so this is lower bound + rand * difference 
%           x2(:,i) = -lb(i) + (ub(i) + lb(i))*rand(num_starts,1);
%         end
%         clear i
%         tpoints = CustomStartPointSet(x2); 
% 
%         
%     
%         problem = createOptimProblem('fmincon','x0',x0,'objective',@(pars) negloglike(Realsize, pars, Einterp, NObserved),'Aineq',A,'bineq',B,'options',opts);
%         [xmin,fmin,exitflag,~,soln] = run(ms,problem,tpoints);
%     
%     
%         %Investigate Solution
%         %open up the soln structure:
%         temp=zeros(num_starts,10);
%         start_points=zeros(num_starts,7);
%         temp=zeros(num_starts,10);
%         c=1;
%         
%         for j=1:length(soln)
%             %check to see if all start points led to an individual solution or
%             %not (MultiSTart will only return unique solutions)
%             g=cell2mat(soln(j).X0);
%             if length(g)==7 %only one start_point led to that solution
%                 start_points(c,:)=g;
%                 temp(c,1:7)=soln(j).X;
%                 temp(c,8)=soln(j).Fval;
%                 %temp(c,9)=growth_rate(Einterp,volbins,dayCOUNTS,temp(c,1:13),hr1,hr2);
%                 temp(c,10)=soln(j).Exitflag;
%                 c=c+1;
%             else
%                 num=length(g)/7;
%                 start_points(c:c+num-1,:)=squeeze(reshape(g',1,7,num))';
%                 temp(c:c+num-1,1:7)=repmat(soln(j).X,num,1);
%                 temp(c:c+num-1,8)=repmat(soln(j).Fval,num,1);
%                 %temp(c:c+num-1,9)=repmat(growth_rate(Einterp,volbins,dayCOUNTS,temp(c,1:13),hr1,hr2),num,1);
%                 temp(c:c+num-1,10)=repmat(soln(j).Exitflag,num,1);
%                 c=c+num;
%             end
%         end
%         %just in case have rows left as zeros
%         qq=find(temp(:,1)~=0);
%         temp=temp(qq,:);
%         start_points=start_points(qq,:);
%         
%         modelfits=[modelfits; temp(1:10];
%         Starts=[Starts; start_points];
%         
%         %okay, in this batch run, did the solver "converge"?
%         [sortlogL ii]=sort(modelfits(:,8));
%        
%         if abs(sortlogL(5)-sortlogL(1)) < icsTol
%             flag1 = 0;
%         else
%             disp(num2str(sortlogL(1:5))) %should be 5, but occassionally get less than 5 solver runs returned...
%             flag1 = 1;
%         end;
% 
%         partol=max(modelfits(ii(1:5),1:7))-min(modelfits(ii(1:5),1:7));
%         if sum(abs(partol) < tolvec)==7 || sum((abs(partol./modelfits(ii(1),1:7)) < 0.05))==7 %either the modelfits are within an absolute tolerance or within a relative tolerance
%             flag2 = 0;
%         else
%             flag2 = 1;
%         end
%         disp(['flag1 = ' num2str(flag1) ' flag2=' num2str(flag2)])
% 
%     end  %while loop

    
    %Assign Outputs
    [s jj]=sort(modelfits(:,8));
    xmin=modelfits(jj(1),1:7);
    fmin=modelfits(jj(1),8);
    exitflag=modelfits(jj(1),9);
    
    %%  Evaluate 
    
    %Simulate according to results of optimization and calculate growth rate

    [Expectation, ~, mu] = Simulate(Realsize, NObserved(:,1), xmin, Einterp); 

    Results = [xmin fmin mu exitflag length(modelfits)];  %this is how the optimization results are stored in model output files 

    

end 


