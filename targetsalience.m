function sal = targetsalience(thisMap,v,nbins)
% Calculate relative target salience for the region delimited by the 
% polygon with the vertices V, using the saliency map THISMAP. The
% calculation uses histograms with NBINS bins. Returns the relative 
% salience of the target region, SAL.

% Normalise the input saliency map
thisMap = thisMap/max(thisMap(:));
thisMap(thisMap < 0) = 0;

% Mask out the target
mask = roipoly(thisMap,v(:,2),v(:,1));
targetPixels = thisMap(mask);
backgroundPixels = thisMap(~mask);

% Find the difference between cumulative histograms for the target and 
% background, as a measure of relative target salience
binWidth = 1/nbins;
bins = linspace(binWidth/2,1-binWidth/2,nbins);
mothHist = hist(targetPixels,bins);
mothHist = mothHist/sum(mothHist);
cumMothHist = cumsum(mothHist);
backgroundHist = hist(backgroundPixels,bins);
backgroundHist = backgroundHist/sum(backgroundHist);
cumBackgroundHist = cumsum(backgroundHist);
sal = sum(cumBackgroundHist - cumMothHist)/nbins;
%plot(bins,cumMothHist,'r-',bins,cumBackgroundHist,'g:');
