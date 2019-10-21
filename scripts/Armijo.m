function [xk,armEvals] = Armijo(xk,fcn,pk,c,r,f,Mx,My,lambda)
% Perform Armijo line search for the function fcn, from the point xk, using
% the update direction pk. We use gradk and c to check the Armijo condition.
% We use r to determine the backtracking rate.
% We assume the initial alpha=1 is used for the backtracking algorithm.

alpha = 1;
armEvals = 2; %we always evaluate fk and fgrad at least once
[fk, fgrad] = fcn(xk,f,Mx,My,lambda);
fkNew = fcn(xk + alpha.*pk,f,Mx,My,lambda);
while fkNew > fk + (c*alpha*fgrad'*pk) %continue changing alpha until armijo condition is met
    fprintf('i = %d', armEvals);
    alpha = r*alpha;
    fkNew = fcn(xk + alpha.*pk,f,Mx,My,lambda);
    armEvals = armEvals + 1;
end

xk = xk + alpha.*pk; %update xk

return





