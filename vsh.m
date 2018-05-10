function har = vsh(kind, l, m, theta, phi, r, k)
% Wektorowe harmoniki sferyczne
% 
%   VSH(kind, l, m, theta, phi, r, k)
% 
%   kind - typ funkcji: 'M' lub 'N'
%   l - stopień (degree)
%   m - rząd (order)
%   phi - katy azymutalne (az)
%   theta - katy zenitalne (el)
%   r - promien, odleglosc od srodka ukladu wspolrzednych
%   k - wektor falowy w danym osrodku
% 
%   Zwraca tablice o wymiarach (phi x theta x 3), gdzie 
%   ostatni wymiar to skladowe (r, el, az) wektora dla 
%   danych wartosci phi i theta.

% wiersze - phi
% kolumny - theta

h = @sphhankel;

if m > l
    error('Blad: m > l')
end
% upewniamy sie ze argumenty sa poziomymi wektorami
if size(theta, 1) > 1
    theta = theta';
end
if size(phi, 1) > 1
    phi = phi';
end
% upewniamy sie, ze k i r sa skalarami
if sum([isscalar(k), isscalar(r)]) < 2
    error('Zle wartosci: k,r musza byc skalarami')
end

kr = k * r;

%%%%%% UWAGA!
%%%%%% m ma wartosci od 0 do l
%%%%%% P jest indeksowane od 1, wiec bierzemy P(m+1) dla oznaczenia funkcji
%%%%%% stowarzyszonej Legendre'a P^m_l

% liczymy vsh
if kind == 'M'
    har = zeros(length(phi), length(theta), 3);
    alpha = ( h(l, kr) .* exp(1j*m*phi) );
    % skladowa r
    har(:,:,1) = alpha' .* zeros(1, length(theta)); 
    % skladowa zenitalna
    har(:,:,2) = 1j * alpha' .* Pi(l, m, theta); 
    % skladowa azymutalna
    har(:,:,3) = -alpha' .* Tau(l, m, theta);
elseif kind == 'N'
    P = legendre(l, cos(theta));
    har = zeros(length(phi), length(theta), 3);
% pomijamy osobne liczenie dh i od razu wlaczamy je do calego wyrazenia
%     dh = n ./ kr .* sphhankel(n, kr) - sphhankel(n+1, kr);
%     h = sqrt(pi ./ (2*kr)) .* besselh(n+1/2, 1, kr);
    alpha = ( h(l, kr) + (l .* sphhankel(l, kr) - kr*sphhankel(l+1, kr))) / kr * exp(1j*m*phi);
    % skladowa r
    har(:,:,1) = l * (l + 1) * h(l, kr) / kr * P(abs(m)+1, :) .* exp(1j*m*phi)';
    % skladowa zenitalna
    har(:,:,2) = alpha' .* Tau(l, m, theta);
    % skladowa azymutalna
    har(:,:,3) = 1j * alpha' .* Pi(l, m, theta);
else
    error('Zly typ funkcji')
end
end

function x = Pi(l, m, theta)
% funkcja Pi
P = legendre(l, cos(theta));
x = m ./ sin(theta) .* P(abs(m)+1, :);
end

function y = Tau(l, m, theta)
% funkcja Tau
x = cos(theta);
% wielomiany
P = legendre(l, x);
P1 = legendre(l+1, x);
% czesciowe wyrazenie na pochodna
dP = (l - m + 1) * P1(abs(m)+1, :) - (l + 1) * x .* P(abs(m)+1, :);
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

function h = sphhankel(l, kr)
% sferyczna funkcja Hankela pierwszego rodzaju
h = sqrt(pi ./ (2*kr)) .* besselh(l+1/2, 1, kr);
% h = sphbessel(1, n, z) + 1j * sphbessel(2, n, z);
end

% function dh = sphhankelderivative(n, kr)
% % pochodna sferycznej funkcji Hankela pierwszego rodzaju
% dh = n ./ kr .* sphhankel(n, kr) - sphhankel(n+1, kr);
% end