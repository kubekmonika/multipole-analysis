function [x, y, z, az, el, rad, theta, phi] = sferawspl(R, n_theta, n_phi)
% Tworzy sferę o zadanym promieniu i zwraca jej wspolrzedne
% w ukladzie kartezjanskim i sferycznym.
% 
%   [x, y, z, az, el, rad, theta, phi] = SFERAWSPL(R, n_theta, n_phi)
% 
%   R - promien sfery
%   n_theta - liczba punktow na jakie dzieli sie kat zenitalny
%   n_phi - liczba punktow na jakie dzieli sie kat azymutalny
% 
%   x, y, z - wspolrzedne kartezjanskie
%   az, el, rad - wspolrzedne sferyczne
%   theta, phi - katy sferyczne

% wiersze - phi - az
% kolumny - theta - el

% %% wspolrzedne zanitalne
% theta = linspace(0, 2*pi, n_theta + 1) + 0.1; 
% theta = theta(1: n_theta);
% az = theta .* ones(n_phi, n_theta);
% az = az(:);
% % wspolrzedne azymutalne
% phi = linspace(-pi/2, pi/2, n_phi + 1) + 0.001;
% phi = phi(1: n_phi);
% el = phi' .* ones(n_phi, n_theta);
% el = el(:);
% % wspolrzedne promienia
% rad = ones(n_phi * n_theta, 1) * R;
% % przechodzimy do wspolrzednych kartezjanskich
% [x, y, z] = sph2cart(az, el, rad);

%% 
% wspolrzedne zanitalne
% theta = linspace(0, pi, n_theta + 1); 
theta = linspace(-pi/2, pi/2, n_theta); 
el = theta .* ones(n_phi, n_theta);
el = el(:);
% wspolrzedne azymutalne
phi = linspace(0, 2*pi, n_phi + 1);
phi = phi(1: n_phi);
az = phi' .* ones(n_phi, n_theta);
az = az(:);
% wspolrzedne promienia
rad = ones(n_phi * n_theta, 1) * R;
% przechodzimy do wspolrzednych kartezjanskich
[x, y, z] = sph2cart(az, el, rad);
end