function Main(image, dest, noise, mu, iter, iterstep, flags)
  % Remove noise from 2d images. Calculates expected value
  % for on the image with iter points, starting at lambda and going up by iterstep
  % image: The original image to have noise added to it
  % dest: Where to store the denoised image
  % noise: Percentage to be used for gaussian noise adding
  % mu: The 'cost' for u not being piecewise const
  % iter: How many expected value points to plot
  % iterstep: How much to increase lambda by each iteration
  % flags: 0 for L2, 1 for Split Bregman

  tol = 5*10^-3;
  im=imread(image);
  im=im2double(im);
  im = im(:);
  matsize = sqrt(length(im));
  rows = length(im);
  I = speye(rows);
  xaxis = mu + (0:iter-1)*iterstep;

  % Add noise 
  imnoise = addnoise(noise, im);

  % Create x 'cost' matrix for optimization
  dia1 = ones(rows,1);
  dia1(rows+1:rows+matsize,1) = 0;
  dia2 = ones(rows+matsize,1);
  Mx = spdiags([-1*dia1 dia2],[0 matsize+1],rows,rows);

  % Create y 'cost' matrix for optimization
  dia1 = ones(rows,1);
  dia1(rows+1,1) = 0;
  dia2 = ones(rows+1,1);
  My = spdiags([-1*dia1 dia2],[0 matsize+1],rows,rows);
  for i=1:matsize-1
    x = i*matsize;
    My(x,x) = 0;
    My(x,x+1) = 0;
  endfor
  
  % Test different lambda 
  for i=1:iter
    
    % L2 Optimization
    if flags == 0
      A = I + lambda*Mx'*Mx + lambda*My'*My;
      uk = A\imnoise;
    endif
      
    % Split Bregman Optimization
    if flags == 1
      uk = imnoise;
      lambda = 2*mu;
      [ukprev, dxk, dyk, bxk, byk] = deal(zeros(length(uk),1));
      while norm(uk-ukprev,2)/norm(uk,2) > tol
        prev = uk;
        uk = gauseidel(uk, dxk, dyk, bxk, byk, imnoise, lambda, mu);
        ##printf("Gaus\n");
        ukprev = prev;
        gradx = (1/256)*Mx*uk;
        grady = (1/256)*My*uk;
        dxk = shrink(gradx + bxk, 1/lambda);
        ##printf("Shrink1\n");
        dyk = shrink(grady + byk, 1/lambda);
        ##printf("Shrink2\n");
        bxk = bxk + gradx - dxk;
        ##printf("bkx\n");
        byk = byk + grady - dyk;
        ##printf("bky\n");
      endwhile
    endif
    
##    hold on;
##    plot(uk);
##    plot(im);
##    legend('Calc Opt','Orig');

    % Change image vector to matrix
    expected = vec2mat(uk, sqrt(length(uk)))';
   % Export image
    if ~exist(dest, 'dir')
      mkdir(dest);
    endif
    fullName = fullfile(dest, [num2str(mu, '%.3f') '.tif']); 
    expected = mat2gray(expected);
    imwrite(expected, fullName);
    
    % Add graph point
    expval(i) = norm(expected(:)-im,2)/norm(im,2);
    
    % Increment lambda
    mu = mu + iterstep;
  endfor
  
  % Create and export graph
  plot(xaxis, expval);
  saveas(gcf, [dest '\plot.png']);
  
  
  function ret = shrink(x, gam)
    % Vectorization
    n = length(x);
    ret(1:n,1) = (x(1:n,1)./abs(x(1:n,1))).*max(abs(x(1:n,1))-gam,0);
    ret(isnan(ret)) = 0; 
    
##    Linear for loop
##    for i=1:length(x)
##      if (x(i) == 0)
##        ret(i) = 0;
##      else
##        ret(i) = x(i)/abs(x(i)) * max(abs(x(i))-gam,0);
##      endif
##    endfor
##    ret = ret';
  endfunction
  
  function ret = gauseidel(uk, dxk, dyk, bxk, byk, noisyim, lambda, mu)
    a = length(uk);
    sa = sqrt(a);
    uk = vec2mat(uk, sa)';
    dxk = vec2mat(dxk, sa)';
    dyk = vec2mat(dyk, sa)';
    bxk = vec2mat(bxk, sa)';
    byk = vec2mat(byk, sa)';
    noisyim = vec2mat(noisyim, sa)';
    
    temp(2:sa-1,2:sa-1) = uk(3:sa,2:sa-1)+uk(1:sa-2,2:sa-1)+uk(2:sa-1,3:sa)+uk(2:sa-1,1:sa-2)...
      +dxk(1:sa-2,2:sa-1)-dxk(2:sa-1,2:sa-1)+dyk(2:sa-1,1:sa-2)-dyk(2:sa-1,2:sa-1)...
      -bxk(1:sa-2,2:sa-1)+bxk(2:sa-1,2:sa-1)-byk(2:sa-1,1:sa-2)+byk(2:sa-1,2:sa-1);
    temp(1,2:sa-1) = uk(2,2:sa-1)+uk(1,3:sa)+uk(1,1:sa-2)...
      -dxk(1,2:sa-1)+dyk(1,1:sa-2)-dyk(1,2:sa-1)...
      +bxk(1,2:sa-1)-byk(1,1:sa-2)+byk(1,2:sa-1);
    temp(sa,2:sa-1) = uk(sa-1,2:sa-1)+uk(sa,3:sa)+uk(sa,1:sa-2)...
      +dxk(sa-1,2:sa-1)-dxk(sa,2:sa-1)+dyk(sa,1:sa-2)-dyk(sa,2:sa-1)...
      -bxk(sa-1,2:sa-1)+bxk(sa,2:sa-1)-byk(sa,1:sa-2)+byk(sa,2:sa-1);
    temp(2:sa-1,1) = uk(3:sa,1)+uk(1:sa-2,1)+uk(2:sa-1,2)...
      +dxk(1:sa-2,1)-dxk(2:sa-1,1)-dyk(2:sa-1,1)...
      -bxk(1:sa-2,1)+bxk(2:sa-1,1)+byk(2:sa-1,1);
    temp(2:sa-1,sa) = uk(3:sa,sa)+uk(1:sa-2,sa)+uk(2:sa-1,sa-1)...
      +dxk(1:sa-2,sa)-dxk(2:sa-1,sa)+dyk(2:sa-1,sa-1)-dyk(2:sa-1,sa)...
      -bxk(1:sa-2,sa)+bxk(2:sa-1,sa)-byk(2:sa-1,sa-1)+byk(2:sa-1,sa);
    temp(1,1) = uk(2,1)+uk(1,2)-dxk(1,1)-dyk(1,1)+bxk(1,1)+byk(1,1);
    temp(1,sa) = uk(2,sa)+uk(1,sa-1)-dxk(1,sa)+dyk(1,sa-1)-dyk(1,sa)+bxk(1,sa)-byk(1,sa-1)+byk(1,sa);
    temp(sa,1) = uk(sa-1,1)+uk(sa,2)+dxk(sa-1,1)-dxk(sa,1)-dyk(sa,1)-bxk(sa-1,1)+bxk(sa,1)+byk(sa,1);
    temp(sa,sa) = uk(sa-1,sa)+uk(sa,sa-1)+dxk(sa-1,sa)-dxk(sa,sa)+dyk(sa,sa-1)-dyk(sa,sa)-bxk(sa-1,sa)...
      +bxk(sa-1,sa)-byk(sa,sa-1)+byk(sa,sa);
    ret = (lambda/(mu+(4*lambda)))*temp + (mu/(mu+4*lambda))*noisyim;
    ret = ret(:);
  endfunction
  
endfunction
