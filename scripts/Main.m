function main(image, dest, noise, mu, iter, iterstep)
  % Remove noise from 2d images. Calculates expected value
  % for on the image with iter points, starting at lambda and going up by iterstep
  % image: The original image to have noise added to it
  % dest: Where to store the denoised image
  % noise: Percentage to be used for gaussian noise adding
  % mu: The 'cost' for u not being piecewise const
  % iter: How many expected value points to plot
  % iterstep: How much to increase lambda by each iteration

  tol = 5*10^-3;
  lambda = 2*mu;
  im=imread(image);
  im=im2double(im);
  im = im(:);
  rows = size(im,1);
  I = speye(rows);
  xaxis = mu + (0:iter-1)*iterstep;

  % Add noise 
  imnoise = addnoise(noise, im);

  % Create x 'cost' matrix for optimization
  Mx = sparse(rows,rows);
  for i = 1:rows-256
    Mx(i,i) = -1;
    Mx(i,i+256) = 1;
  endfor

  % Create y 'cost' matrix for optimization
  My = sparse(rows,rows);
  for i = 1:rows-1
    if (mod(i,256) == 0)
      My(i,i) = 0;
      My(i,i+1) = 0;
    else 
      My(i,i) = -1;
      My(i,i+1) = 1;
    endif
  endfor
  
  % Test different lambda 
  for i=1:iter
    
##    % L2 Optimization
##    A = I + lambda*Mx'*Mx + lambda*My'*My;
##    uk = A\imnoise;
      
    % Split Bregman Optimization
    uk = imnoise;
    [ukprev, dxk, dyk, bxk, byk] = deal(zeros(length(uk),1));
    while norm(uk-ukprev,2)/norm(uk,2) > tol
      norm(uk-ukprev,2)/norm(uk,2)
      temp = uk;
      uk = gauseidel(uk, dxk, dyk, bxk, byk, imnoise, lambda, mu);
      printf("Gaus\n");
      ukprev = temp;
      gradx = Mx*uk;
      grady = My*uk;
      dxk = shrink(gradx + bxk, 1/lambda);
      printf("Shrink1\n");
      dyk = shrink(grady + byk, 1/lambda);
      printf("Shrink2\n");
      bxk = bxk + gradx - dxk;
      printf("bkx\n");
      byk = byk + grady - dyk;
      printf("bky\n");
    endwhile

    % Change image vector to matrix
    expected = vec2mat(uk, sqrt(length(uk)));
    expected = expected';
      
   % Export image
    if ~exist(dest, 'dir')
      mkdir(dest);
    endif
    fullName = fullfile(dest, [num2str(lambda, '%.3f') '.tif']); 
    expected = mat2gray(expected);
    imwrite(expected, fullName);
    
    % Add graph point
    expval(i) = norm(expected(:)-im,2)/norm(im,2);
    
    % Increment lambda
    lambda = lambda + iterstep;
  endfor
  
  % Create and export graph
  plot(xaxis, expval);
  saveas(gcf, [pwd '\' dest '\plot.png']);
  
  
  function ret = shrink(x, gam)
    for i=1:length(x)
      if (x(i) == 0)
        ret(i) = 0;
      else
        ret(i) = x(i)/abs(x(i)) * max(abs(x(i))-gam,0);
      endif
    endfor
    ret = ret';
  endfunction
  
  function ret = gauseidel(uk, dxk, dyk, bxk, byk, noisyim, lambda, mu)
    a = length(uk);
    sa = sqrt(a);
    for i = 1:a
      if (i <= sa)
        ukup = 0;
        dxkup = 0;
        bxkup = 0;
      else
        ukup = uk(i-sa);
        dxkup = dxk(i-sa);
        bxkup = bxk(i-sa);
      endif
      
      if (i > a-sa)
        ukdown = 0;
      else 
        ukdown = uk(i+sa);
      endif
      
      if (mod(i,sa) == 1)
        ukleft = 0;
        dykleft = 0;
        bykleft = 0;
      else
        ukleft = uk(i-1);
        dykleft = dyk(i-1);
        bykleft = byk(i-1);
      endif
      
      if (mod(i,sa) == 0)
       ukright = 0;
      else
       ukright = uk(i+1);
      endif
      ret(i) = (lambda/(mu+(4*lambda)))*(ukup+ukdown+ukright+ukleft+dxkup-dxk(i)+dykleft-dyk(i)-bxkup+bxk(i)-bykleft+byk(i))+(mu/(mu+4*lambda))*noisyim(i);
    endfor
    ret = ret';
  endfunction
  
endfunction
