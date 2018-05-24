%% Skyrpt w ktorym obliczamy wspolczynniki rozproszenia.

%% Najpierw potrzebujemy dane ze skryptu pole_sfera.m
% R - promien sfery
% enei - dlugosc fali
% e - pole
% theta, phi - katy sferyczne
% az, el - wspolrzedne sferyczne

%% obliczamy pole dla calego przekroju dlugosci fali

[E, R, eneis, N, theta, phi] = field_sphere_si(20, 25);

% the_end() % funkcja daje sygnal, ze obliczenia sie skonczyly

%% obliczamy wspolczynniki rozproszenia
C = coefforeverywavelength(E, R, eneis, N, theta, phi);

%% wykres
figure
plot(eneis, abs(C))

%% SPRAWDZAMY VSH

% wavevector = @(n, enei) sqrt(4 * pi^2 * n^2 / (enei * 0.1^9)^2);

%%
l = 1;
m = 0;
dim = 500;
wavevector = @(n, enei) 2 * pi *n / (enei * 0.1^9);
p = trisphere( 144, 200 );
[x, y, z, ~,~,~, ~, ~] = sferawspl(R, 20, 25);
%%
N_ = vsh('N', l, m, theta, phi, R, wavevector(1,600));
vshN = zeros(dim, 3);
vshN(:,1) = reshape(N_(:,:,1), dim, []);
vshN(:,2) = reshape(N_(:,:,2), dim, []);
vshN(:,3) = reshape(N_(:,:,3), dim, []);
%
M_ = vsh('M', l, m, theta, phi, R, wavevector(1,600));
vshM = zeros(dim, 3);
vshM(:,1) = reshape(M_(:,:,1), dim, []);
vshM(:,2) = reshape(M_(:,:,2), dim, []);
vshM(:,3) = reshape(M_(:,:,3), dim, []);

dot(vshM, vshN)
% clearvars 'dim' 'l' 'm'

%%
vshN = zeros(dim, 3, length(eneis));
vshM = zeros(dim, 3, length(eneis));

for i = 1 : length(eneis)
    N_ = vsh('N', l, m, theta, phi, R, wavevector(1, eneis(i)));
    vshN(:,1,i) = reshape(N_(:,:,1), dim, []);
    vshN(:,2,i) = reshape(N_(:,:,2), dim, []);
    vshN(:,3,i) = reshape(N_(:,:,3), dim, []);
    %
    M = vsh('M', l, m, theta, phi, R, wavevector(1,eneis(i)));
    vshM(:,1,i) = reshape(M(:,:,1), dim, []);
    vshM(:,2,i) = reshape(M(:,:,2), dim, []);
    vshM(:,3,i) = reshape(M(:,:,3), dim, []);
end
%%
e = E(:, :, 1);
e2 = E(:, :, 15);
%%
hold on
plot(eneis, acos(squeeze(dot(E(100,:,:), E(100,:,:), 2)))./squeeze(dot(E(100,:,:), E(100,:,:), 2)))
plot(eneis, acos(squeeze(dot(E(10,:,:), E(10,:,:), 2)))./squeeze(dot(E(10,:,:), E(10,:,:), 2)))
plot(eneis, acos(squeeze(dot(E(200,:,:), E(200,:,:), 2)))./squeeze(dot(E(200,:,:), E(200,:,:), 2)))
hold off
%%
hold on
plot(eneis, squeeze(dot(vshN(100,:,:), E(100,:,:), 2)))
plot(eneis, squeeze(dot(vshN(10,:,:), E(10,:,:), 2)))
plot(eneis, squeeze(dot(vshN(200,:,:), E(200,:,:), 2)))
hold off
%%
hold on
plot(eneis, squeeze(dot(vshN(10,:,:), vshN(10,:,:), 2)))
plot(eneis, squeeze(dot(vshN(100,:,:), vshN(100,:,:), 2)))
plot(eneis, squeeze(dot(vshN(200,:,:), vshN(200,:,:), 2)))
hold off
%%
hold on
skl = 3;
plot(eneis, abs(squeeze(vshN(10,skl,:))))
plot(eneis, abs(squeeze(vshN(100,skl,:))))
plot(eneis, abs(squeeze(vshN(200,skl,:))))
hold off
%%
hold on
skl = 3;
% plot(eneis, abs(squeeze(E(10,skl,:))))
plot(eneis, abs(squeeze(E(100,skl,:))))
plot(eneis, abs(squeeze(E(200,skl,:))))
hold off
%%
hold on
% plot(eneis, abs(squeeze(E(10,skl,:))))
x = abs(conj(squeeze(E(100,1,:))).*squeeze(vshN(100,1,:))+...
    conj(squeeze(E(100,2,:))).*squeeze(vshN(100,2,:))+...
    conj(squeeze(E(100,3,:))).*squeeze(vshN(100,3,:)));
plot(eneis, x)
plot(eneis, abs(dot(squeeze(E(100,:,:)), squeeze(vshN(100,:,:)))))
% plot(eneis, abs(squeeze(E(200,skl,:)).*squeeze(vshN(200,skl,:))))
hold off

%%
hold on
for m = [-1, 0, 1]
    vshN = zeros(dim, 3, length(eneis));
    vshM = zeros(dim, 3, length(eneis));

    for i = 1 : length(eneis)
        N_ = vsh('N', l, m, theta, phi, R, wavevector(1, eneis(i)));
        vshN(:,1,i) = reshape(N_(:,:,1), dim, []);
        vshN(:,2,i) = reshape(N_(:,:,2), dim, []);
        vshN(:,3,i) = reshape(N_(:,:,3), dim, []);
        %
        M = vsh('M', l, m, theta, phi, R, wavevector(1,eneis(i)));
        vshM(:,1,i) = reshape(M(:,:,1), dim, []);
        vshM(:,2,i) = reshape(M(:,:,2), dim, []);
        vshM(:,3,i) = reshape(M(:,:,3), dim, []);
    end
    % plot(eneis, abs(squeeze(E(10,skl,:))))
    x = abs(conj(squeeze(E(100,1,:))).*squeeze(vshN(100,1,:))+...
        conj(squeeze(E(100,2,:))).*squeeze(vshN(100,2,:))+...
        conj(squeeze(E(100,3,:))).*squeeze(vshN(100,3,:)));
    plot(eneis, x)
%     plot(eneis, abs(dot(squeeze(E(100,:,:)), squeeze(vshN(100,:,:)))))
end
hold off
%%
figure;
hold on

quiver3(x, y, z, e(:,1), e(:,2), e(:,3))
quiver3(x, y, z, e2(:,1), e2(:,2), e2(:,3))

% plot(p)
% quiver3(x, y, z, vshN(:,1), vshN(:,2), vshN(:,3))
% quiver3(x, y, z, vshM(:,1), vshM(:,2), vshM(:,3))

xlabel( 'x (nm)' );
ylabel( 'y (nm)' );
zlabel( 'z (nm)' );
title('Pole elektryczne wokol nanoczastki')

% axis equal tight
hold off