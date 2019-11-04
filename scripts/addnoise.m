function noisyim = addnoise(noise, im)

  % Add noise to image
  v = noise*var(im);
  noisyim = imnoise(im, 'gaussian', 0, v);
  noisyim = 255.*noisyim;
  
endfunction
