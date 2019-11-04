function [fval,fgrad] = l1optfunction(xk, f, Mx, My, lambda)
% Function for image optimization with l1 norm
% xk: Vector we want to become the original image
% f: Noisy image vector
% Mx: Cost matrix with respect to x axis
% My: Cost matrix with respect to y axis
% lambda: Cost multiplier

x = Mx*xk;
y = My*xk;
xgrad = lambda*(Mx'*Mx)*xk;
ygrad = lambda*(My'*My)*xk;

fval = (.5*norm(xk-f)^2) + (lambda*norm(x,1)) + (lambda*norm(y,1));
fgrad = xk - f + xgrad + ygrad;

endfunction
