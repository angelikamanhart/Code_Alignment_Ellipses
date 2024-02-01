function [x,y] = evenlySpacedPoints(N,sq,min_dist)

% Define how many points you want to place.
maxPoints = N;
% Define counter for how many actually get placed.
numPoints = 0;
% Define fail safe, how many iterations you want to keep trying before quitting.
maxIterations = 1000 * maxPoints;
loopCounter = 1;
% Define how close points can come to other points before being rejected.
minClosestDistance = min_dist;
% Declare arrays to hold the x and y coordinate values.
x = nan(1, numPoints);
y = nan(1, numPoints);
while numPoints < maxPoints && loopCounter < maxIterations
  % Get a random point.
  xPossible = sq*rand();
  yPossible = sq*rand();
  if numPoints == 0
    % First point automatically is valid.
    numPoints = numPoints + 1;
    x(numPoints) = xPossible;
    y(numPoints) = yPossible;
    continue;
  end
  % Find distances between this point and all others.
  distances = sqrt((x-xPossible) .^ 2 + (y - yPossible) .^ 2);
  if min(distances) >= minClosestDistance
    % It's far enough away from all the other points, or it's the first point.
    % Include it in our list of acceptable points.
    numPoints = numPoints + 1;
    x(numPoints) = xPossible;
    y(numPoints) = yPossible;
  end
  % Increment the loop counter.
  loopCounter = loopCounter + 1;
end
% Crop away points we were not able to get.
x = x(1:numPoints); % Crop to valid points.
y = y(1:numPoints); % Crop to valid points.

end
