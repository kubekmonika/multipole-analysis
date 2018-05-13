%% Skyrpt w ktorym obliczamy wspolczynniki rozproszenia.

%% Najpierw potrzebujemy dane ze skryptu pole_sfera.m
% R - promien sfery
% enei - dlugosc fali
% e - pole
% theta, phi - katy sferyczne
% az, el - wspolrzedne sferyczne

%% obliczamy pole dla calego przekroju dlugosci fali

[E, R, eneis, N, theta, phi] = field_sphere_si(25, 25);

% the_end() % funkcja daje sygnal, ze obliczenia sie skonczyly

%% obliczamy wspolczynniki rozproszenia
C = coefforeverywavelength(E, R, eneis, N, theta, phi);

%% wykres
figure
plot(eneis(5:end), C(5:end))

%% SPRAWDZAMY VSH

% wavevector = @(n, enei) sqrt(4 * pi^2 * n^2 / (enei * 0.1^9)^2);

%%
l = 2;
m = 1;
dim = 625;
wavevector = @(n, enei) 2 * pi *n / (enei * 0.1^9);
p = trisphere( 144, 200 );
[x, y, z, ~,~,~, theta, phi] = sferawspl(R, sqrt(dim), sqrt(dim));

N = vsh('N', l, m, theta, phi, R, wavevector(1,600));
vshN = zeros(dim, 3);
vshN(:,1) = reshape(N(:,:,1), dim, []);
vshN(:,2) = reshape(N(:,:,2), dim, []);
vshN(:,3) = reshape(N(:,:,3), dim, []);
%
M = vsh('M', l, m, theta, phi, R, wavevector(1,600));
vshM = zeros(dim, 3);
vshM(:,1) = reshape(M(:,:,1), dim, []);
vshM(:,2) = reshape(M(:,:,2), dim, []);
vshM(:,3) = reshape(M(:,:,3), dim, []);

dot(vshM, vshN)
clearvars 'dim' 'l' 'm'
%%
e = E(:, :, 1);

%%
figure;
hold on

quiver3(x, y, z, e(:,1), e(:,2), e(:,3))

% plot(p)
quiver3(x, y, z, vshN(:,1), vshN(:,2), vshN(:,3))
quiver3(x, y, z, vshM(:,1), vshM(:,2), vshM(:,3))

xlabel( 'x (nm)' );
ylabel( 'y (nm)' );
zlabel( 'z (nm)' );
title('Pole elektryczne wokol nanoczastki')

% axis equal tight
hold off