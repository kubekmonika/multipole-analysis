function [x, y, z, w, r, theta, phi] = getCoordinates(n_points, R)
% Calculates coordinates of the sphere built from a given number
% of points and having a given radius.
% 
%   GETCOORDINATES(n_points, R)
%   
%   n_points - number of points on the sphere: { 6, 14, 26, 38, 50, 74, 86, 110, 
%       146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 
%       1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 },
%       if provided with other value then it is rounded to the nearest one
%   R - radius [nm]
%
%   x, y, z - spherical coordinates
%   w - weights for each point (necessary for integration)
%   r - radius [m]
%   theta, phi - spherical coordinates of the points [rad]

% number of points
n = getNPoints(n_points);
% cartesian coordinates
leb = getLebedevSphere(n);
% weights
w = leb.w;
[az, el, ~] = cart2sph(leb.x, leb.y, leb.z);
% transform to sphere with a given radius
r = R * 1e-9;
R = ones(size(az)) * R;
[x, y, z] = sph2cart(az, el, R);
[~, theta, phi] = cartToSph(x, y, z);
end

function n = getNPoints(n_points)
% Return number of points needed for calculating getLebedevSphere
points = [ 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, ...
    350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, ...
    3470, 3890, 4334, 4802, 5294, 5810 ];
if ismember(n_points, points)
    n = n_points;
else
    diff = min(abs(points - n_points));
    if ismember(n_points - diff, points)
        n = n_points - diff;
    elseif ismember(n_points + diff, points)
        n = n_points + diff;
    else
        error('Something went wrong: wrong value of n_points')
    end
end
end