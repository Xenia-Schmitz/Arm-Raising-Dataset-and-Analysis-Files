function [pathLength dist_sq] = calPathlength(x, y, tOn, tOff, tStep)
%
%Input
%
%x: x raw data
%y: y raw data
%tOn: beginning of calculation
%tOff: end of calculation
%tStep: downsampling factor
%
%Output
%pathLength: path length (excursion) between tOn and tOff
%Length:
if length(x)~=length(y),
    error('The length of x is different from the length of y.')
end

if nargin == 2,
    tOn = 1;
    tOff = length(x);
    tStep = 1;
end

newX = x(tOn:tStep:tOff);
newY = y(tOn:tStep:tOff);

diffX = diff(newX);
diffY = diff(newY);

diff_sq_X = diffX.^2;
diff_sq_Y = diffY.^2;

xy_diff_sq = diff_sq_X + diff_sq_Y;

dist_sq = sqrt(xy_diff_sq);

pathLength = sum(dist_sq);





