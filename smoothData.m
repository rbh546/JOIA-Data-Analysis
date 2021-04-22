function smoothedData = smoothData( data,windowSize )
%SMOOTHDATA simple smoothing function (moving avg.).
%   SMOOTHEDDATA = SMOOTHDATA( DATA,WINDOWSIZE ) smooths DATA with
%   corresponding WINDOWSIZE and outputs the smoothed data SMOOTHEDDATA.
%
%   NOTE: if DATA is a matrix, the function proceeds by column, i.e., one 
%   column at a time (considers data is ordered by column).

%   INPUTS:
%   DATA: Unsmoothed data.
%   WINDOWSIZE: Size of the averaging window.
% 
%   EXAMPLE:
%   smoothedData = smoothData( data,10 );
% 
%   Martin Richard
%   C-CORE, St. John's, NL CANADA

%   Copyright 2011-2012 C-CORE

window = ones(windowSize,1)/windowSize; 
smoothedData = nan(size(data));
for i=1:size(data,2)
    smoothedData(:,i) = convn(data(:,i),window,'same');
end
smoothedData(1:windowSize/2,:)          = nan; % unreliable smoothe portion
smoothedData(end-windowSize/2+1:end,:)  = nan; % unreliable smoothe portion