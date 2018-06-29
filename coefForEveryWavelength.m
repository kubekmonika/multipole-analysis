function [sca, scb] = coefForEveryWavelength(E, r, eneis, N, theta, phi, w)
%Funkcja liczaca wspolczynniki rozproszenia dla danego pola.
%   
%   COEFFOREVERYWAVELENGTH(E, r, eneis, N, theta, phi, w)
% 
%   E - rozkład pola we współrzędnych sferycznych
%   r - promień sfery
%   eneis - dlugości fali
%   N - współczynnik załamania w otaczającym ośrodku; macierz wartości
%       dla każdej długości fali lub jedna wartość skalarna; wszystkie
%       wartości muszą być rzeczywiste
%   theta - współrzędne zenitalne, 0 <= theta <= pi
%   phi - współrzędne azymutalne, 0 <= phi <= 2pi
%   w - wagi dla każdego punktu

assert(isreal(N), 'N musi być rzeczywiste - teoria nie przewiduje innych przypadków')
% funkcja liczaca wektor falowy
wavevector = @(n, enei) 2 * pi * n / ((enei) * 0.1^9);

if isscalar(N)
    N = ones(size(eneis)) * N;
end

n_eneis = length(eneis);
% obliczamy wspolczynniki rozproszenia
sca = zeros(n_eneis, 1);
scb = zeros(n_eneis, 1);
% wartosc dla kazdej dlugosci fali
for i = 1: n_eneis
    k = wavevector(N(i), eneis(i));
    e = E(:,:,i);
    sca(i) = scattCoefficientA(e, k, r, theta, phi, w);
    scb(i) = scattCoefficientB(e, k, r, theta, phi, w);
end
end
