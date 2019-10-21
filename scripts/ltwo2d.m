function ltwo2d(image, dest, name, noise, lambda)
% Remove noise from 2d images using l2 norm.
% image: The original image to have noise added to it
% dest: Where to store the denoised image
% name: What to call the denoised image
% noise: Percentage to be used for gaussian noise adding
% lambda: The 'cost' for u not being piecewise const

im=imread(image);
im=im2double(im);
im = im(:);

% Add noise to image
v = noise*var(im);
I_noisy = imnoise(im, 'gaussian', 0, v);
I_noisy=255.*I_noisy;

rows = size(I_noisy,1);
I = speye(rows);

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
  end
endfor

% Optimize image
A = I + lambda*Mx'*Mx + lambda*My'*My;
x = A\I_noisy;

% Change image vector to matrix
u = vec2mat(x, sqrt(length(x)));
u = u';


% Export image
if ~exist(dest, 'dir')
  mkdir(dest);
end
fullName = fullfile(dest, name); 
u = mat2gray(u);
imwrite(u, fullName);

endfunction