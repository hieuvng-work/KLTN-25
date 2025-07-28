function [powerAllocation,waterLevel] = functionWaterfilling(totalPower,lambdaInv)
Nt = length(lambdaInv); %Extract number of singular directions
lambdaInvSorted = sort(lambdaInv,'ascend'); %Sort lambdaInv in ascending order

alpha_candidates = (totalPower+cumsum(lambdaInvSorted))./(1:Nt)'; %Compute different values on the Lagrange multiplier (i.e., waterlevel) given that 1,2,...,Nt of the singular directions get non-zero power
optimalIndex = alpha_candidates-lambdaInvSorted(1:end,1)>0 & alpha_candidates-[lambdaInvSorted(2:end,1); Inf]<0; %Find the true Lagrange multiplier alpha by checking which one of the candidates that only turns on the singular directions that are supposed to be on
waterLevel = alpha_candidates(optimalIndex); %Extract the optimal Lagrange multiplier (i.e., waterlevel in the waterfilling analogy)

powerAllocation = waterLevel-lambdaInv; %Compute power allocation
powerAllocation(powerAllocation<0) = 0; %Make sure that inactive singular directions receive zero power

end
