%%sprawdzamy nowe harmoniki sferyczne
wavevector = @(n, enei) 2 * pi * n / ((enei) * 1e-9);
%%
[x, y, z, az, el, rad, theta, phi] = sferawspl(205, 51, 41);
%%
K = eneis * 0;
for i = 1 : length(eneis)
    K(i) = wavevector(1, eneis(i));
end
%%
theta_ = theta;
phi_ = phi;
%%
theta_ = linspace(0, pi, 101); 
phi_ = linspace(0, 2*pi, 91);
%%
k = wavevector(1, eneis(10));
i = 10;
X = squeeze( vsh('N', 1, 0, theta_, phi_(i), r, k) );
Y = squeeze( vshN(0, theta_, phi_(i), r, k) );

figure
hold on
subplot(3,1,1)
plot(theta_, X(:, 1), theta_, Y(:,1), '-.')
subplot(3,1,2)
plot(theta_, dot(X,X,2), theta_, dot(Y,Y,2), '-.')
subplot(3,1,3)
plot(theta_, dot(X,X,2)' .* sin(theta_), theta_, dot(Y,Y,2)' .* sin(theta_), '-.')

%%
m = -1;
licznik = zeros(1, length(eneis));
mianownik = zeros(1, length(eneis));

for j = 1 : length(eneis)
    k = wavevector(1, eneis(j));
    e = sphreshapefield(E(:,:,j), length(theta), length(phi));
    Y = squeeze( vshM(m, theta, phi, r, k) );
%     Y1 = squeeze( vshN(m+1, theta, phi, 200, k) );
    dotXE = dot(Y,e,3) .* sin(theta);
    dotXX = dot(Y,Y,3) .* sin(theta);
    licznik(j) = sum(dotXE(:));
    mianownik(j) = sum(dotXX(:));
end
wynik = abs(licznik./mianownik) .^ 2;

figure
hold on
subplot(2,1,1)
yyaxis left
plot(eneis, mianownik)
yyaxis right
plot(eneis, licznik)
legend('mianownik', 'licznik')

subplot(2,1,2)
% plot(eneis, dot(wynik, wynik, 1), 'g-')
plot(eneis, wynik, 'g-')
legend('| licz/mian |^2')

%%
Y_ = dot(Ym0, Ym1, 3);
sum(Y_(:))
Y_ = dot(Ym1, Ym01, 3);
sum(Y_(:))
Y_ = dot(Ym0, Ym01, 3);
sum(Y_(:))
%%
j = 1; % dlugosc fali
i = 40; % kat phi

e = sphreshapefield(E(:,:,j), length(theta), length(phi));
e = squeeze( e(i,:,:));

k = wavevector(1, eneis(j));

X = squeeze( vsh('N', 1, 1, theta, phi(i), r, k) );
Y = squeeze( vshN(1, theta, phi(i), r, k) );

figure
hold on

subplot(3,2,1)
plot(theta, X(:, 1), theta, Y(:,1), '-.')
title('vshN')
legend('stare N', 'nowe N')

subplot(3,2,3)
plot(theta, dot(X,X,2), theta, dot(Y,Y,2), '-.')
title('dot(N,N)')

subplot(3,2,4)
plot(theta, dot(X,X,2)' .* sin(theta), theta, dot(Y,Y,2)' .* sin(theta), '-.')
title('dot(N,N) * sin(theta)')

subplot(3,2,5)
plot(theta, dot(X,e,2), theta, dot(Y,e,2), '-.')
title('dot(N,E)')

subplot(3,2,6)
plot(theta, dot(X,e,2)' .* sin(theta), theta, dot(Y,e,2)' .* sin(theta), '-.')
title('dot(N,E) * sin(theta)')

%%
EE = zeros(length(phi), length(theta), 3, length(eneis));
for j = 1 : length(eneis)
    EE(:,:,:,j) = sphreshapefield(Ekart(:,:,j), length(theta), length(phi));
end

dotE = zeros(1, length(eneis));
for i = 1 : length(eneis)
    temp = squeeze(EE(:,:,:,i));
    temp = dot(temp, temp, 3);
    dotE(i) = sum(temp(:));
end
figure
hold on
yyaxis left
plot(eneis(15:end), dotE(15:end), '--', 'LineWidth', 2)
%%
% figure
yyaxis right
hold on
% 61 x 71
for pkt = 1 : 15 : 50
    pktE = zeros(1, length(eneis));
    for i = 1 : length(eneis)
        temp = squeeze(EE(10,pkt,:,i));
        pktE(i) = dot(temp, temp);
    end
    plot(eneis, pktE)
end

%%
h = @(l, kr) sqrt(pi ./ (2*kr)) .* besselh(l+1/2, 1, kr); 
h_appr = @(l, kr) (-1j)^l * exp(1j*kr) / (1j * kr);

H = 0 * K;
H_appr = 0 * K;
for i = 1 : length(K)
    H(i) = h(1, K(i)*r);
    H_appr(i) = h_appr(1, K(i)*r);
end

figure
plot(K, H, 'r-', K, H_appr, 'b--')

%%
enei = eneis(50);
r = R * 1e-9;
k = wavevector(1, enei);

X = squeeze( vshN(1, theta, phi, r, k) );

dotX = dot(X, X, 3);
dotXsin = dotX .* sin(theta + pi/2);

cmaxr = max(real(X(:)));
cminr = min(real(X(:)));
cmaxi = max(imag(X(:)));
cmini = min(imag(X(:)));

figure

ax1 = subplot(3,6,[1,2]);
pcolor(theta, phi, real(X(:,:,1))); colorbar;
caxis(ax1, [cminr, cmaxr]);

ax2 = subplot(3,6,[3,4]);
pcolor(theta, phi, real(X(:,:,2))); colorbar;
caxis(ax2, [cminr, cmaxr]);
title('Real X')

ax3 = subplot(3,6,[5,6]);
pcolor(theta, phi, real(X(:,:,3))); colorbar;
caxis(ax3, [cminr, cmaxr]);

ax4 = subplot(3,6,[7,8]);
pcolor(theta, phi, imag(X(:,:,1))); colorbar; colormap(ax4, 'autumn');
caxis(ax4, [cmini, cmaxi]);

ax5 = subplot(3,6,[9,10]);
pcolor(theta, phi, imag(X(:,:,2))); colorbar; colormap(ax5, 'autumn');
caxis(ax5, [cmini, cmaxi]);
title('Imag X')

ax6 = subplot(3,6,[11,12]);
pcolor(theta, phi, imag(X(:,:,3))); colorbar; colormap(ax6, 'autumn');
caxis(ax6, [cmini, cmaxi]);

ax7 = subplot(3,6,[13,14,15]);
pcolor(theta, phi, dotX); colorbar; colormap(ax7, 'summer'); 
title('dot X')

ax8 = subplot(3,6,[16,17,18]);
pcolor(theta, phi, dotXsin); colorbar; colormap(ax8, 'winter');
title('dot X * sin')

arrayfun(@(s) set(s,'EdgeColor','none'), findobj(gcf,'type','surface'))

%%
enei = eneis(50);
r = R * 1e-9;
k = wavevector(1, enei);

X = squeeze( vshN(0, theta, phi, r, k) );
% Ershp = sphreshapefield(E(:,:,50), length(theta), length(phi));
Y = squeeze( vshM(1, theta, phi, r, k) );

dotX = real(dot(X, Y, 3));
dotXsin = real(dotX .* sin(theta));

% dotX = real(dotY ./ dotX);
% dotXsin = real(dotX .* sin(theta));

cmaxr = max(real(X(:)));
cminr = min(real(X(:)));
cmaxi = max(imag(X(:)));
cmini = min(imag(X(:)));

figure

ax1 = subplot(3,6,[1,2]);
pcolor(theta, phi, real(X(:,:,1))); colorbar;
caxis(ax1, [cminr, cmaxr]);

ax2 = subplot(3,6,[3,4]);
pcolor(theta, phi, real(X(:,:,2))); colorbar;
caxis(ax2, [cminr, cmaxr]);
title('Real X')

ax3 = subplot(3,6,[5,6]);
pcolor(theta, phi, real(X(:,:,3))); colorbar;
caxis(ax3, [cminr, cmaxr]);

ax4 = subplot(3,6,[7,8]);
pcolor(theta, phi, imag(X(:,:,1))); colorbar; colormap(ax4, 'autumn');
caxis(ax4, [cmini, cmaxi]);

ax5 = subplot(3,6,[9,10]);
pcolor(theta, phi, imag(X(:,:,2))); colorbar; colormap(ax5, 'autumn');
caxis(ax5, [cmini, cmaxi]);
title('Imag X')

ax6 = subplot(3,6,[11,12]);
pcolor(theta, phi, imag(X(:,:,3))); colorbar; colormap(ax6, 'autumn');
caxis(ax6, [cmini, cmaxi]);

ax7 = subplot(3,6,[13,14,15]);
pcolor(theta, phi, dotX); colorbar; colormap(ax7, 'summer'); 
title('dot X')

ax8 = subplot(3,6,[16,17,18]);
pcolor(theta, phi, dotXsin); colorbar; colormap(ax8, 'winter');
title('dot X * sin')

arrayfun(@(s) set(s,'EdgeColor','none'), findobj(gcf,'type','surface'))

%%
for i = -1:1:1
    for j = -1:1:1
        X = squeeze( vshN(i, theta, phi, r, k) );
        Y = squeeze( vshM(j, theta, phi, r, k) );
        dotX = dot(Y, X, 3);
        dotXsin = dotX .* sin(theta);
        disp({i, j, sum(dotXsin(:))})
    end
end

%%
wavevector = @(n, enei) 2 * pi * n / ((enei) * 1e-9);
K = eneis * 0;
for i = 1 : length(eneis)
    K(i) = wavevector(1, eneis(i));
end
%%
i = 1;
k = K(i);
e = E(:,:,i);
[a, A] = scattcoefficienta(1, k, e, r, theta, phi);
[b, B] = scattcoefficientb(1, k, e, r, theta, phi);

%%
[A(3) - A(1); 1j * (A(1) + A(3)); -sqrt(2) * A(2)]

%    0.0017 + 0.3107i  -0.0016 + 0.0013i   0.4333 - 0.2123i

%%
figure;
hold on
[x, y, z, ~,~,~, ~, ~] = sferawspl(R, length(theta), length(phi));
N1 = vshN(-1, theta, phi, r, k);
N2 = vshN(0, theta, phi, r, k);
N3 = vshN(1, theta, phi, r, k);
N_ = N3 - N1;
N_x = reshape(N_(:,:,1),[],1);
N_y = reshape(N_(:,:,2),[],1);
N_z = reshape(N_(:,:,3),[],1);
% quiver3(x, y, z, e(:,1), e(:,2), e(:,3))
quiver3(x, y, z, N_x, N_y, N_z)

xlabel( 'x (nm)' );
ylabel( 'y (nm)' );
zlabel( 'z (nm)' );