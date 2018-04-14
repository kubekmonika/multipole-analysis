%% Skyrpt w ktorym obliczamy wspolczynniki rozproszenia.

%% Najpierw potrzebujemy dane ze skryptu pole_sfera.m
% R - promien sfery
% enei - dlugosc fali
% e - pole
% theta, phi - katy sferyczne
% az, el - wspolrzedne sferyczne

load('pole_sfera2.mat', 'R', 'enei', 'e', 'theta', 'phi', 'az', 'el')
%% obliczamy wektor falowy
% predkosc swiatla w prozni
c = 299792458; % [m/s]
% czestosc kolowa fali
lambda = enei * 0.1^9; % [m]
% przenikalnosc elektryczna i magnetyczna prozni
epsilon0 = 8.854187817620 * 0.1^12; % [A^2 s^4 / kg m^3]
mu0 = 1 / (epsilon0 * c^2);
% wspolczynni zalamania n^2 = |epsilon_r * mu_r|
% https://refractiveindex.info/?shelf=main&book=Si&page=Green-2008
n = 3.695; % wartosc dla krzemu, dlugosc fali 776nm
% wektor falowy k^2 = 4 pi^2 n^2 / lambda^2
k = sqrt(4 * pi^2 * n^2 / lambda^2); % [1/m]

%% obliczamy pole dla calego przekroju dlugosci fali
filename = 'C:\Users\Monika\Documents\MATLAB\analiza_multipolowa\krzem_n.txt';
[E, R, eneis, n, theta, phi] = field_sphere_si(100, 100, filename);
the_end()
%% obliczamy wspolczynniki rozproszenia
C = coefforeverywavelength(E, R, eneis, n, theta, phi);

%% wykres
plot(eneis, C)