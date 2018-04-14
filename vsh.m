function har = vsh(kind, n, m, theta, phi, r, k)
% Wektorowe harmoniki sferyczne
% 
%   VSH(kind, n, m, theta, phi, r, k)
% 
%   kind - typ funkcji: 'M' lub 'N'
%   n - stopień (degree)
%   m - rząd (order)
%   theta - katy azymutalne (az)
%   phi - katy zenitalne (el)
%   r - promien, odleglosc od srodka ukladu wspolrzednych
%   k - wektor falowy w danym osrodku
% 
%   Zwraca tablice o wymiarach (phi x theta x 3), gdzie 
%   ostatni wymiar to skladowe (r, az, el) wektora dla 
%   danych wartosci phi i theta.

h = @sphhankel;
dh = @sphhankelderivative;
% upewniamy sie ze argumenty sa poziomymi wektorami
if m > n
    error('Blad: m > n')
end
if size(theta, 1) > 1
    theta = theta';
end
if size(phi, 1) > 1
    phi = phi';
end
if sum([isscalar(k), isscalar(r)]) < 2
    error('Zle wartosci: k,r musza byc skalarem')
end
kr = k * r;

% modyfikujemy funkcje Legendre'a, zeby uwzglednic m<0
global lgr
if m < 0
    lgr = (-1)^abs(m) * gamma(n-abs(m)+1) / gamma(n+abs(m)+1);
else
    lgr = 1;
end

% liczymy vsh
if kind == 'M'
    har = zeros(length(phi), length(theta), 3);
    alpha = ( h(n, kr) .* exp(1j*m*phi) );
    % skladowa r
    har(:,:,1) = alpha' .* zeros(1, length(theta)); 
    % skladowa azymutalna
    har(:,:,2) = 1j * alpha' .* Pi(n, m, theta); 
    % skladowa zenitalna
    har(:,:,3) = -alpha' .* Tau(n, m, theta);
elseif kind == 'N'
    P = lgr * legendre(n, cos(theta));
    har = zeros(length(phi), length(theta), 3);
    alpha = ( h(n, kr) + kr * dh(n, kr) ) / kr * exp(1j*m*phi);
    % skladowa r
    har(:,:,1) = n * (n + 1) * h(n, kr) / kr * P(abs(m)+1, :) .* exp(1j*m*phi)';
    % skladowa azymutalna
    har(:,:,2) = alpha' .* Tau(n, m, theta);
    % skladowa zenitalna
    har(:,:,3) = 1j * alpha' .* Pi(n, m, theta);
else
    error('Zly typ funkcji')
end
end

function x = Pi(n, m, theta)
global lgr
P = lgr * legendre(n, cos(theta));
x = m ./ sin(theta) .* P(abs(m)+1, :);
end

function y = Tau(n, m, theta)
global lgr
x = cos(theta);
% wielomiany
P = lgr * legendre(n, x);
P1 = lgr * legendre(n+1, x);
% czesciowe wyrazenie na pochodna
dP = (n - m + 1) * P1(abs(m)+1, :) - (n + 1) * x .* P(abs(m)+1, :);
% wynik
y = - sin(theta) .* dP ./ (x.^2 - 1);
end

% function j = sphbessel(kind, n, z)
% % sferyczne funkcje Bessel'a
% if kind == 1
%     j = sqrt(pi ./ (2*z)) .* besselj(n+1/2, z);
% elseif kind == 2
%     j = sqrt(pi ./ (2*z)) .* bessely(n+1/2, z);
% end
% end

function h = sphhankel(n, z)
% sferyczna funkcja Hankela pierwszego rodzaju
h = sqrt(pi ./ (2*z)) .* besselh(n, 1, z);
% h = sphbessel(1, n, z) + 1j * sphbessel(2, n, z);
end

function dh = sphhankelderivative(n, z)
% pochodna sferycznej funkcji Hankela pierwszego rodzaju
dh = n ./ z .* sphhankel(n, z) - sphhankel(n+1, z);
end