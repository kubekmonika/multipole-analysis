function [C_1, C_2] = coefForEveryWavelength(E, r, eneis, N, theta, phi, type)
%Funkcja liczaca wspolczynniki rozproszenia dla danego pola.
%   
%   COEFFOREVERYWAVELENGTH(E, r, eneis, N, theta, phi)
% 
%   E - rozkład pola we współrzędnych sferycznych
%   r - promień sfery
%   eneis - dlugości fali
%   N - współczynnik załamania w otaczającym ośrodku; macierz wartości
%       dla każdej długości fali lub jedna wartość skalarna; wszystkie
%       wartości muszą być rzeczywiste
%   theta - współrzędne zenitalne, 0 <= theta <= pi
%   phi - współrzędne azymutalne, 0 <= phi <= 2pi
%   type - typ współczynnika: 'a', 'b' lub obydwa 'ab'

% N -> wspolczynnik zalamania, musi byc rzeczywisty bo innych przypadkow
% teoria w tej metodzie nie przewiduje.
assert(isreal(N), 'N musi być rzeczywiste')
% funkcja liczaca wektor falowy
wavevector = @(n, enei) 2 * pi * n / ((enei) * 0.1^9);

if isscalar(N)
    N = ones(size(eneis)) * N;
end

n_eneis = length(eneis);
% obliczamy wspolczynniki rozproszenia
if type == 'ab'
    C_1 = zeros(n_eneis, 1);
    C_2 = zeros(n_eneis, 1);
    % wartosc dla kazdej dlugosci fali
    for i = 1: n_eneis
        k = wavevector(N(i), eneis(i));
        e = E(:,:,i);
        C_1(i) = scattcoefficienta(e, k, r, theta, phi);
        C_2(i) = scattcoefficientb(e, k, r, theta, phi);
    end
elseif type == 'a'
    C_1 = zeros(n_eneis, 1);
    C_2 = nan;
    % wartosc dla kazdej dlugosci fali
    for i = 1: n_eneis
        k = wavevector(N(i), eneis(i));
        e = E(:,:,i);
        C_1(i) = scattcoefficienta(e, k, r, theta, phi);
    end
elseif type == 'b'
    C_1 = zeros(n_eneis, 1);
    C_2 = nan;
    % wartosc dla kazdej dlugosci fali
    for i = 1: n_eneis
        k = wavevector(N(i), eneis(i));
        e = E(:,:,i);
        C_1(i) = scattcoefficientb(e, k, r, theta, phi);
    end
end
end
