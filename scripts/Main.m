function [xk,fval,fgrad,normGrad,numEvals] = Main(xk,fcn,lnsrch,iterations)
    % Perform iterations of the Steepest Descent or DFP method
    % and reports back the the final point, value, gradient,
    % norm and number of evaluations. 
    % Note that xk must be a column vector if it is vector-valued.

    if strcmp(lnsrch, 'Steepest Descent') 
        [xk,fval,fgrad,normGrad,numEvals] = SteepestDescent(xk,fcn,iterations);
    end
    
    if strcmp(lnsrch, 'DFP')
       [xk,fval,fgrad,normGrad,numEvals] = DFP(xk,fcn,iterations);
    end
    
    return
