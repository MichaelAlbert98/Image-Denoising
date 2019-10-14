function ltwo(image, dest, name, noise, lambda)
% Remove noise from 1d image using l2 norm.
% image: The original image to have noise added to it
% noise: Percentage to be used for gaussian noise adding
% lambda: The 'cost' for u not being piecewise const
im=imread(image);
im=im2double(im);

% Add noise to image
v = noise*var(im(:));
I_noisy = imnoise(im, 'gaussian', 0, v);
I_noisy=255.*I_noisy;

rows = size(I_noisy,2);
I = eye(rows);
M = zeros(rows);

% Create 'cost' matrix for optimization
for i = 1:rows-1
  M(i,i) = -1;
  M(i,i+1) = 1;
endfor
M(i+1,i+1) = 0;

% Optimize image line by line (potentially not correct?)
for j = 1:rows
  x = (I+(lambda*M'*M))^-1*I_noisy(j,:)';
  x = x';
  u(j,:) = x;
endfor

% Export image
if ~exist(dest, 'dir')
  mkdir(dest);
end
fullName = fullfile(dest, name); 
u = mat2gray(u);
imwrite(u, fullName);

endfunction