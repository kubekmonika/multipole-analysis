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

% wiersze - phi - el
% kolumny - theta - az

%% wspolrzedne azymutalne
theta = linspace(1e-4, 2*pi-1e-4, az_n);
az = theta .* ones(el_n, az_n);
az = az(:);
% wspolrzedne zenitalne
phi = linspace(-pi/2, pi/2, el_n);
el = phi' .* ones(el_n, az_n);
el = el(:);
% wspolrzedne promienia
r = ones(az_n * el_n, 1) * R;
% przechodzimy do wspolrzednych kartezjanskich
[x, y, z] = sph2cart(az, el, r);
end