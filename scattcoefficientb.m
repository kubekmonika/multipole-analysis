function b = scattcoefficientb(l, k, e, R, theta, phi)
%Oblicza wspolczynnik rozproszenia b(l).
%
%   SCATTCOEFFICIENTA(l, m, E, R, theta, phi)
% 
%   l - stopien, l >= |m|
%   k - wartosc wektora falowego
%   e - rozklad pola (w postaci wektora)
%   R - promien sfery
%   phi - katy zenitalne
%   theta - katy azymutalne

b = 0;
for m = -l : 1 : l
    blm = b_lm(l, m, k, e, R, theta, phi);
    b = b + k^2 * l * (l + 1) * (blm * conj(blm));
end
end

function blm = b_lm(l, m, k, e, R, theta, phi)
% Zwraca wartosc wspolczynnika rozproszenia a_lm
assert(l >= abs(m), 'Blad: |m| > l')
% VSH
% N = vsh('N', n, m, theta, phi, R, k);
% if m >= 0
%     M = vsh('M', l, m, theta, phi, R, k);
% else
%     M = (-1)^abs(m) * conj(vsh('M', l, abs(m), theta, phi, R, k));
% end
M = vshM(m, theta, phi, R, k);

% przeksztalcamy pole
E = sphreshapefield(e, length(theta), length(phi));
% zastepujemy NaN wartoscia 0
M(isnan(M)) = 0;
E(isnan(E)) = 0;
% licznik
licznik = dot(M, E, 3) .* sin(theta);
mianownik = dot(M, M, 3) .* sin(theta);
blm = sum(licznik(:)) / sum(mianownik(:));
end