function [orig, noisyim] = AddNoise(noise, image)
  % noise: Percentage to be used for gaussian noise adding
  % image: The image to add noise to
  % noise: Percentage to be used for gaussian noise adding
  
  orig=imread(image);
  orig=im2double(orig);
  orig = orig(:);
  
  % Add noise to image
  v = noise*var(orig);
  noisyim = imnoise(orig, 'gaussian', 0, v);
  
endfunction
