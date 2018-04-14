function [x, y, z, az, el, r, theta, phi] = sferawspl(R, az_n, el_n)
% Tworzy sferę o zadanym promieniu i zwraca jej wspolrzedne
% w ukladzie kartezjanskim i sferycznym.
% 
%   [x, y, z, az, el, r, theta, phi] = SFERAWSPL(R, az_n, el_n)
% 
%   R - promien sfery
%   az_n - liczba punktow na jakie dzieli sie kat azymutalny
%   el_n - liczba punktow na jakie dzieli sie kat zenitalny
% 
%   x, y, z - wspolrzedne kartezjanskie
%   az, el, r - wspolrzedne sferyczne
%   theta, phi - katy sferyczne

% wspolrzedne azymutalne
az = ones(az_n, el_n);
theta = linspace(0, 2*pi, az_n);
for i = 1 : az_n
    az(:, i) = theta(i);
end
az = az(:);
% wspolrzedne zenitalne
el = ones(az_n, el_n);
phi = linspace(-pi/2, pi/2, el_n);
for i = 1 : el_n
    el(i, :) = phi(i);
end
el = el(:);
% wspolrzedne promienia
r = ones(az_n * el_n, 1) * R;
% przechodzimy do wspolrzednych kartezjanskich
[x, y, z] = sph2cart(az, el, r);
end