function relerr = Main(image, dest, noise, mu, iter, iterstep, flags)
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
  yaxis = mu + (0:iter-1)*iterstep;

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
  
  % L2 Optimization
  if flags == 0
    % Test different mu
    for i=1:iter
      A = I + mu*Mx'*Mx + mu*My'*My;
      uk = A\imnoise;
      exportimL2(uk,dest,mu);
      mu = mu + iterstep;
      relerr(i) = norm(uk-im,2)/norm(im,2);
    endfor
    graph2d(xaxis,relerr,dest);
  endif
    
  % Split Bregman Optimization
  if flags == 1 || flags == 2
    startval = mu;
    % Test different mu
    for i=1:iter
      lambda = startval;
      % Test different lambda
      for j=1:iter
        uk = imnoise;
        loops = 0;
        [ukprev, dxk, dyk, bxk, byk] = deal(zeros(length(uk),1));
        
        while norm(uk-ukprev,2)/norm(uk,2) > tol
          loops = loops + 1;
          prev = uk;
          uk = gauseidel(uk, dxk, dyk, bxk, byk, imnoise, lambda, mu);
          ukprev = prev;
          gradx = (1/256)*Mx*uk;
          grady = (1/256)*My*uk;
          % Anisotropic
          if flags == 1
            dxk = shrink(gradx + bxk, 1/lambda);
            dyk = shrink(grady + byk, 1/lambda);
          % Isotropic
          else
            sk = sqrt(abs(gradx+bxk).^2+abs(grady+byk).^2);
            dxk = max(sk - 1/lambda,0).*((gradx+bxk)./sk);
            dyk = max(sk - 1/lambda,0).*((grady+byk)./sk);
          endif
          bxk = bxk + gradx - dxk;
          byk = byk + grady - dyk;
        endwhile
        relerr(i,j) = norm(uk-im,2)/norm(im,2);
        totaliters(i,j) = loops;
        exportimbreg(uk, dest, mu, lambda);
        lambda = lambda + iterstep;
      endfor
      mu = mu + iterstep;
    endfor
    graphheat(xaxis, yaxis, relerr, 'relative error', [dest '/relerr.png']);
    graphheat(xaxis, yaxis, totaliters, 'total iterations', [dest '/totaliters.png']);
  endif

    
    
  function exportimL2(imvec,dest,mu)  
    % Change image vector to matrix
    expected = vec2mat(imvec, sqrt(length(imvec)))';
    % Export image
    if ~exist(dest, 'dir')
      mkdir(dest);
    endif
    fullName = fullfile(dest, [num2str(mu, '%.3f') '.tif']); 
    expected = mat2gray(expected);
    imwrite(expected, fullName);
  endfunction    
  
  function graph2d(x,y,dest)
    plot(x, y);
    saveas(gcf, [dest '\plot.png']);
  endfunction
  
  function exportimbreg(imvec,dest,mu,lambda)
    % Change image vector to matrix
    expected = vec2mat(imvec, sqrt(length(imvec)))';
    % Export image
    if ~exist(dest, 'dir')
      mkdir(dest);
    endif
    fullName = fullfile(dest, ['mu' num2str(mu, '%.3f') 'lam' num2str(lambda, '%.3f') '.tif']); 
    expected = mat2gray(expected);
    imwrite(expected, fullName);
  endfunction
  
  function graphheat(x,y,vals,name,dest)
    imagesc(x,y,vals);
    colormap('hot');
    colorbar
    xlabel('lambda');
    ylabel('mu');
    title(name);
    set(gca,'YDir','normal');
    saveas(gcf, dest);
  endfunction
  
  function ret = shrink(x, gam)
    % Vectorization
    n = length(x);
    ret = (x./abs(x)).*max(abs(x)-gam,0);
    %ret(1:n,1) = (x(1:n,1)./abs(x(1:n,1))).*max(abs(x(1:n,1))-gam,0);
    ret(isnan(ret)) = 0; 
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
