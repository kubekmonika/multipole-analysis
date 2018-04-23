% %% Skyrpt w ktorym obliczamy wspolczynniki rozproszenia.
% 
% %% Najpierw potrzebujemy dane ze skryptu pole_sfera.m
% % R - promien sfery
% % enei - dlugosc fali
% % e - pole
% % theta, phi - katy sferyczne
% % az, el - wspolrzedne sferyczne
% 
% load('pole_sfera2.mat', 'R', 'enei', 'e', 'theta', 'phi', 'az', 'el')
% %% obliczamy wektor falowy
% % predkosc swiatla w prozni
% c = 299792458; % [m/s]
% % czestosc kolowa fali
% lambda = enei * 0.1^9; % [m]
% % przenikalnosc elektryczna i magnetyczna prozni
% epsilon0 = 8.854187817620 * 0.1^12; % [A^2 s^4 / kg m^3]
% mu0 = 1 / (epsilon0 * c^2);
% % wspolczynni zalamania n^2 = |epsilon_r * mu_r|
% % https://refractiveindex.info/?shelf=main&book=Si&page=Green-2008
% n = 3.695; % wartosc dla krzemu, dlugosc fali 776nm
% % wektor falowy k^2 = 4 pi^2 n^2 / lambda^2
% k = sqrt(4 * pi^2 * n^2 / lambda^2); % [1/m]

%% obliczamy pole dla calego przekroju dlugosci fali
% filename = 'C:\Users\Monika\Documents\MATLAB\analiza_multipolowa\krzem_n.txt';
[E, R, eneis, n, theta, phi] = field_sphere_si(50, 50);
the_end() % funkcja daje sygnal, ze obliczenia sie skonczyly

%% zmiana wspolrzednych E z kartezjanskich na sferyczne
Esph = 0 * E;
for i = 1 : size(E,3)
    % r
    Esph(:,1,i) = sqrt(E(:,1,i).^2 + E(:,2,i).^2 + E(:,3,i).^2);
    % theta
    Esph(:,2,i) = atan(E(:,2,i) ./ E(:,1,i));
    % phi
    Esph(:,3,i) = acos(E(:,3,i) ./ Esph(:,1,i));
end

%% obliczamy wspolczynniki rozproszenia
C = coefforeverywavelength(Esph, R, eneis, n, theta, phi);

%% wykres
plot(eneis, C)

%%
p = trisphere( 144, 200 );
[x, y, z, ~,~,~, theta, phi] = sferawspl(R, 50, 50);
% wavevector = @(n, enei) sqrt(4 * pi^2 * n^2 / (enei * 0.1^9)^2);
wavevector = @(n, enei) 2 * pi / (enei * 0.1^9);
%%
N = vsh('N', 1, 0, theta, phi, R, wavevector(1,600));
vshN = zeros(50 * 50, 3);
vshN(:,1) = reshape(N(:,:,1), 2500, []);
vshN(:,2) = reshape(N(:,:,2), 2500, []);
vshN(:,3) = reshape(N(:,:,3), 2500, []);
%%
M = vsh('M', 1, -1, theta, phi, R, wavevector(1,600));
vshM = zeros(50 * 50, 3);
vshM(:,1) = reshape(M(:,:,1), 2500, []);
vshM(:,2) = reshape(M(:,:,2), 2500, []);
vshM(:,3) = reshape(M(:,:,3), 2500, []);

%%
figure;
hold on
% e = E(:,:,31);
e = vshN;

plot(p)
quiver3(x, y, z, e(:,1), e(:,2), e(:,3))

xlabel( 'x (nm)' );
ylabel( 'y (nm)' );
zlabel( 'z (nm)' );
title('Pole elektryczne wokol nanoczastki')

% axis equal tight
hold off