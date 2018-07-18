function [r, theta, phi] = cartToSph(X, Y, Z)
% Transform cartesian coordinates to spherical.
%
%   CARTTOSPH(x, y, z)
%   
%   X, Y, Z - cartesian coordinates
%
%   Returns spherical coordinates:
%   0 < r < oo
%   0 < theta < pi
%   0 < phi < 2pi

r = sqrt( X.^2 + Y.^2 + Z.^2 );
theta = acos( Z./r );
phi = myAtan2(Y, X);
end

function myatan2 = myAtan2(Y, X)
% The atan2 function, which outcome is positive for X and Y coordinates
myatan2 = zeros(size(Y));
for i = 1 : length(X)
    x = X(i);
    y = Y(i);
    if x > 0
        result = 2 * atan( y / (sqrt(x^2 + y^2) + x));
    elseif (x <= 0) && (y ~= 0)
        result = 2 * atan( (sqrt(x^2 + y^2) - x) / y);
    elseif (x < 0) && (y == 0)
        result = pi;
    else
        result = nan;
    end
    if result < 0
        result = result + 2 * pi;
    end
    myatan2(i) = result;
end
end