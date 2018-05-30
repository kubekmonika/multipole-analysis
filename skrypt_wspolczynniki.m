%% Skyrpt w ktorym obliczamy wspolczynniki rozproszenia.

%% Najpierw potrzebujemy dane ze skryptu pole_sfera.m
% R - promien sfery
% enei - dlugosc fali
% e - pole
% theta, phi - katy sferyczne
% az, el - wspolrzedne sferyczne

%% obliczamy pole dla calego przekroju dlugosci fali
[Ekart, R, eneis, N, theta, phi] = field_sphere_si(41, 21);
r = R * 1e-9;
% theta = theta - pi/2;

%% zmiana wspolrzednych E z kartezjanskich na sferyczne
E = 0 * Ekart;
for i = 1 : size(Ekart,3)
    % r
    E(:,1,i) = sqrt(Ekart(:,1,i).^2 + Ekart(:,2,i).^2 + Ekart(:,3,i).^2);
    % theta
    E(:,2,i) = acos(Ekart(:,3,i) ./ E(:,1,i));
    % phi
    E(:,3,i) = atan(Ekart(:,2,i) ./ Ekart(:,1,i));
end
% the_end() % funkcja daje sygnal, ze obliczenia sie skonczyly
%%
EE = zeros(length(phi), length(theta), 3, length(eneis));
for j = 1 : length(eneis)
    EE(:,:,:,j) = sphreshapefield(E(:,:,j), length(theta), length(phi));
end
dotE = zeros(1, length(eneis));
for i = 1 : length(eneis)
    temp = squeeze(EE(:,:,:,i));
    temp = dot(temp, temp, 3) .* sin(theta);
    dotE(i) = sum(temp(:));% / K(i)^2;
end

% dotEkart = dotE;
%% obliczamy wspolczynniki rozproszenia
[C_a, C_b] = coefforeverywavelength(E, r, eneis, N, theta, phi);

%% wykres
od = 5;

figure
hold on
subplot(2,1,1)
yyaxis left
plot(eneis(od:end), C_a(od:end), '-')
yyaxis right
plot(eneis(od:end), C_b(od:end), '--')
legend('a1', 'b1')

subplot(2,1,2)
yyaxis left
plot(eneis(od:end), C_a(od:end) + C_b(od:end), '--')
yyaxis right
plot(eneis(od:end), dotE(od:end), '-', 'LineWidth', 2)
legend('suma', 'dot(E,E)')

%%
figure
hold on
plot(eneis(od:end), dotEkart(od:end), 'r-.', 'LineWidth', 2)
% plot(eneis(od:end), dotEsph(od:end), 'b--', 'LineWidth', 2)
% plot(eneis(od:end), dotE2sph(od:end), 'g-.', 'LineWidth', 1)
plot(eneis(od:end), dotE(od:end), 'k--', 'LineWidth', 1)
% legend('kart', 'sph')
%% wykres
figure
% plot(eneis, abs(C_a))
yyaxis left
plot(eneis, C_a, 'r-')
yyaxis right
plot(eneis, C_b, 'b-')

legend('a1', 'b1')

%%
cart2sphvec = @(az, el) [-sin(az),          cos(az),         0;...
                         -sin(el)*cos(az), -sin(el)*sin(az), cos(el);...
                          cos(el)*cos(az),  cos(el)*sin(az), sin(el)];

[~, ~, ~, az, el, ~, ~, ~] = sferawspl(R, length(theta), length(phi));

E = 0 * Ekart;
for i = 1 : size(E, 3)
    for j = 1 : size(E,1)
        E(j,:,i) = cart2sphvec(az(j), el(j)) * Ekart(j,:,i)';
    end
end
%%
X = sphreshapefield(Ekart(:,:,30), length(theta), length(phi));
dotX = dot(X,X,3);
Y = sphreshapefield(E(:,:,30), length(theta), length(phi));
dotY = dot(Y,Y,3);
figure
ax1 = subplot(2,2,1);
pcolor(dotX); colorbar; colormap(ax1, 'default');
title('norma E kart')
ax2 = subplot(2,2,2);
pcolor(dotY); colorbar; colormap(ax2, 'default');
title('norma E sph')
ax3 = subplot(2,2,[3,4]);
pcolor(dotX - dotY); colorbar; colormap(ax3, 'hot');
title('roznica norm')

arrayfun(@(s) set(s,'EdgeColor','none'), findobj(gcf,'type','surface'))
%%
E = 0 * Ekart;
for i = 1 : size(Ekart,3)
    % r
    E(:,1,i) = sqrt(real(Ekart(:,1,i)).^2 + real(Ekart(:,2,i)).^2 + real(Ekart(:,3,i)).^2);
    % theta
    E(:,2,i) = acos(real(Ekart(:,3,i)) ./ real(E(:,1,i)));
    % phi
    E(:,3,i) = atan(real(Ekart(:,2,i)) ./ real(Ekart(:,1,i)));
    
        % r
    E(:,1,i) = E(:,1,i) + 1j * sqrt(imag(Ekart(:,1,i)).^2 + imag(Ekart(:,2,i)).^2 + imag(Ekart(:,3,i)).^2);
    % theta
    E(:,2,i) = E(:,2,i) + 1j * acos(imag(Ekart(:,3,i)) ./ imag(E(:,1,i)));
    % phi
    E(:,3,i) = E(:,3,i) + 1j * atan(imag(Ekart(:,2,i) ./ imag(Ekart(:,1,i))));
end
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
e = Ekart(:, :, 1);
e2 = Ekart(:, :, 15);
%%
hold on
plot(eneis, acos(squeeze(dot(Ekart(100,:,:), Ekart(100,:,:), 2)))./squeeze(dot(Ekart(100,:,:), Ekart(100,:,:), 2)))
plot(eneis, acos(squeeze(dot(Ekart(10,:,:), Ekart(10,:,:), 2)))./squeeze(dot(Ekart(10,:,:), Ekart(10,:,:), 2)))
plot(eneis, acos(squeeze(dot(Ekart(200,:,:), Ekart(200,:,:), 2)))./squeeze(dot(Ekart(200,:,:), Ekart(200,:,:), 2)))
hold off
%%
hold on
plot(eneis, squeeze(dot(vshN(100,:,:), Ekart(100,:,:), 2)))
plot(eneis, squeeze(dot(vshN(10,:,:), Ekart(10,:,:), 2)))
plot(eneis, squeeze(dot(vshN(200,:,:), Ekart(200,:,:), 2)))
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
plot(eneis, abs(squeeze(Ekart(100,skl,:))))
plot(eneis, abs(squeeze(Ekart(200,skl,:))))
hold off
%%
hold on
% plot(eneis, abs(squeeze(E(10,skl,:))))
x = abs(conj(squeeze(Ekart(100,1,:))).*squeeze(vshN(100,1,:))+...
    conj(squeeze(Ekart(100,2,:))).*squeeze(vshN(100,2,:))+...
    conj(squeeze(Ekart(100,3,:))).*squeeze(vshN(100,3,:)));
plot(eneis, x)
plot(eneis, abs(dot(squeeze(Ekart(100,:,:)), squeeze(vshN(100,:,:)))))
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
    x = abs(conj(squeeze(Ekart(100,1,:))).*squeeze(vshN(100,1,:))+...
        conj(squeeze(Ekart(100,2,:))).*squeeze(vshN(100,2,:))+...
        conj(squeeze(Ekart(100,3,:))).*squeeze(vshN(100,3,:)));
    plot(eneis, x)
%     plot(eneis, abs(dot(squeeze(E(100,:,:)), squeeze(vshN(100,:,:)))))
end
hold off

%%
e = Ekart(:,:,1);
e_resh = sphreshapefield(e, length(theta), length(phi));
figure
contourf(sqrt(dot(e_resh, e_resh, 3)),15)
%%
wavevector = @(n, enei) 2 * pi * n / (enei * 0.1^9); 
k = wavevector(1,eneis(1));
l = 1;
m = 1;

vshN_m0 = vsh('N', l, 0, theta, phi, R, k);
vshN_m1 = vsh('N', l, 1, theta, phi, R, k);
vshN_m2 = (-1)^abs(-1) * conj(vsh('N', l, abs(-1), theta, phi, R, k));

c = scattcoefficienta(1, k, e, R, theta, phi);

%%
prodNE = dot(vshN_m1, e_resh, 3);
hold on
subplot(2,1,1);
contourf(real(prodNE))
subplot(2,1,2);
contourf(real(prodNE .* sin(theta)))
hold off

%%
prodNN = dot(vshN_m2, vshN_m2, 3);
hold on
subplot(2,1,1);
contourf(prodNN)
subplot(2,1,2);
contourf(prodNN .* sin(theta))
hold off

%%
i = 11;
k = wavevector(1, eneis(i));
disp([a_lm(1, -1, k, Ekart(:,:,i), R, theta, phi), ...
    a_lm(1, 0, k, Ekart(:,:,i), R, theta, phi),...
    a_lm(1, 1, k, Ekart(:,:,i), R, theta, phi)])

%%
% ls = linspace(0, 1, 100);
ls = linspace(min(R*K), max(R*K), 1000);
f = @real;
figure
hold on
% plot(ls, f(besselh(1, ls)));
plot(ls, f(besselh(3/2, ls)));
% plot(ls, f(besselh(2, ls)));
hold off
%%
plot(theta, Pi(1, -1, theta))

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