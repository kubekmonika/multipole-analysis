function [C_a, C_b] = coefforeverywavelength(E, r, eneis, N, theta, phi)
%Funkcja liczaca wspolczynniki rozproszenia dla danego pola.

% funkcja liczaca wektor falowy: k^2 = 4 pi^2 n^2 / lambda^2
% wavevector = @(n, enei) sqrt(4 * pi^2 * n^2 / (enei * 0.1^9)^2);
% E - rozklad pola, wywolanie we wspolrzednych sferycznych!! - na pewno?
% c = 299792458; %[m/s] predkosc swiatla w prozni
% N - wspolczynnik zalamania, musi byc rzeczywisty bo innych przypadkow
%     teoria w tej metodzie nie przewiduje
wavevector = @(n, enei) 2 * pi * n / ((enei) * 0.1^9);     % KS: dopisałam mnożenie przez n

% obliczamy wspolczynniki rozproszenia
n_N = length(N);
C_a = zeros(1, n_N);
C_b = zeros(1, n_N);

% petla for / parfor
for i = 1: n_N
    k = wavevector(N(i), eneis(i));
    e = E(:,:,i);
    C_a(i) = scattcoefficienta(1, k, e, r, theta, phi);
    C_b(i) = scattcoefficientb(1, k, e, r, theta, phi);
end

end
