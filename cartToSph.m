function [r, theta, phi] = cartToSph(X, Y, Z)
% Funkcja przeliczająca współrzędne z kartezjańskich na sferyczne.
%
%   CART2SPH(x, y, z)
%   
%   X, Y, Z - wspołrzędne sferyczne
%
%   Zwraca współrzędne:
%   0 < r < oo
%   0 < theta < pi
%   0 < phi < 2pi

r = sqrt( X.^2 + Y.^2 + Z.^2 );
theta = acos( Z./r );
phi = myatan2(Y, X);
end

function A = myatan2(Y, X)
% Funkcja atan2 dająca dodatni wynik dla {x,y} będących wartościami
% skalarnymi.
A = zeros(size(Y));
for i = 1 : length(X)
    x = X(i);
    y = Y(i);
    if x > 0
        a = 2 * atan( y / (sqrt(x^2 + y^2) + x));
    elseif (x <= 0) && (y ~= 0)
        a = 2 * atan( (sqrt(x^2 + y^2) - x) / y);
    elseif (x < 0) && (y == 0)
        a = pi;
    else
        a = nan;
    end
    if a < 0
        a = a + 2 * pi;
    end
    A(i) = a;
end
end