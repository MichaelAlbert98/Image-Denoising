function [orig, noisyim] = addnoise(noise, image, flags)
  % noise: Percentage to be used for gaussian noise adding
  % image: The image to add noise to
  % noise: Percentage to be used for gaussian noise adding
  
  orig=imread(image);
  orig=im2double(orig);
  
  if flags == 0
    % Add Gaussian noise to image
    orig = orig(:);
    v = noise*var(orig);
    noisyim = imnoise(orig, 'gaussian', 0, v);
  elseif flags == 1
    % Add salt and pepper noise to image
    noisyim = imnoise(orig, 'salt & pepper', noise);
    orig = orig(:);
    noisyim = noisyim(:);
  else
    % Add noise to image based on distance from center
    len = length(orig);
    cent = len/2;
    for i=1:len
      for j=1:len
        dist = sqrt((cent-i)^2+(cent-j)^2);
        noise = -.5 + rand;
        noisyim(i,j) = orig(i,j) + (dist/(len))*rand;
      endfor
    endfor
    orig = orig(:);
    noisyim = noisyim(:);
  endif
  
endfunction
