function [xk,fval,fgrad,normGrad,numEvals] = SteepestDescent(xk,f,Mx,My,lambda,fcn,iterations)
% Perform a single iteration of the Steepest Descent method
% and report back the updated value and how many times the
% line search method within this routine was called.
c = 10^-4;
tol = 10^-8;
r = .5;
[fval,fgrad] = fcn(xk,f,Mx,My,lambda);
normGrad = 0;
numEvals = 2; %First iteration has 2 extra evaluations
pk = -fgrad; %Compute descent direction

for i=1:iterations
    
  if abs(norm(pk)) > tol %Do line search if pk is above tolerance
    [xk,armEvals] = Armijo(xk,fcn,pk,c,r,f,Mx,My,lambda);
    [fval,fgrad] = fcn(xk,f,Mx,My,lambda); %Get updated fval
    pk = -fgrad;
    normGrad = norm(fgrad);
    numEvals = 1 + armEvals + numEvals;
    fprintf('Number of Iterations = %d \n', i);
    fprintf('----- \n');
    fprintf('Function value = %.5f \n', fval);
    fprintf('Gradient norm = %.5f \n', normGrad);
    fprintf('Number of Evaluations = %d \n', numEvals);
    
    numEvals = 0; %Reset numEvals after each iteration
  end

end

return
    
