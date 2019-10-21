function [fval,fgrad] = l1optfunction(xk, f, Mx, My, lambda)
% Function for image optimization with l1 norm
% xk: Vector we want to become the original vector
% f: Noisy image vector
% Mx: Cost matrix with respect to x axis
% My: Cost matrix with respect to y axis
% lambda: Cost multiplier

x = Mx*xk;
y = My*xk;
xgrad = zeros(length(xk),1);
ygrad = zeros(length(xk),1);
for i=1:length(xk)
  if x(i) < 0
    xgrad(i) = -1;
  elseif x(i) > 0
    xgrad(i) = 1;
  else 
    xgrad(i) = 0;
  endif

  if y(i) < 0
    ygrad(i) = -1;
  elseif y(i) > 0
    ygrad(i) = 1;
  else 
    ygrad(i) = 0;
  endif
endfor

fval = (.5*norm(xk-f)^2) + (lambda*norm(x,1)) + (lambda*norm(y,1));
fgrad = xk - f + lambda*xgrad + lambda*ygrad;

endfunction
