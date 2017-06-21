function insulin_threshold = calc_insulin_threshold(img)
  img = imgaussfilt(double(img),10);
  [f,xi] = ksdensity(insulin(:),'bandwidth',1);
  spacing = xi(2)-xi(1); % Store spacing between sliding window jumps of ksdensity

  % calc second derivative of the ksdensity curve
  dx=gradient(f,spacing);
  dxx=gradient(dx,spacing);
  plot(xi,dxx)

  % The threshold is set as the x value of the maximum peak to the right of the max peak of second derivative of ksdensity of insulin pixels
  [pks, locs] = findpeaks(dxx);
  locs = ((locs-1)*spacing)+min(xi); % Set correct range of x values. The derivative necessitates the -1.
  [max_peak,max_peak_loc]=max(f);
  max_peak_loc = ((max_peak_loc-1)*spacing)+min(xi);
  pks = pks(locs>max_peak_loc);
  locs = locs(locs>max_peak_loc);
  [max_peak,max_peak_loc]=max(pks);
  thresh = locs(max_peak_loc);

end