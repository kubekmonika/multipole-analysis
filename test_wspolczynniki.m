%% Skrypt w ktorym obliczamy wspolczynniki rozproszenia.
%%
E(isnan(E)) = 0;
dotEkart = squeeze( sum( squeeze( dot(Ekart, Ekart, 2) ) .* theta, 1) );
dotE = squeeze( sum( squeeze( dot(E, E, 2) ) .* theta, 1 ) );

%%
figure
hold on

plot(eneis, dotEkart , 'r-.', 'LineWidth', 2)
% plot(eneis(od:end), dotEsph(od:end), 'b--', 'LineWidth', 2)
% plot(eneis(od:end), dotE2sph(od:end), 'g-.', 'LineWidth', 1)
plot(eneis, dotE, 'k--', 'LineWidth', 1)
% legend('kart', 'sph')

%%
e = Ekart(:, :, 5);
k = K(5);
e2 = Ekart(:, :, 15);

%%
figure;
hold on

plot3(x, y, z, 'Ob')
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

%%
figure
hold on

[sca, scb] = coefForEveryWavelength(E, r, eneis, N, theta, phi, 1);
yyaxis left
plot(eneis, sca , 'k-', 'LineWidth', 3)
plot(eneis, scb , 'b-', 'LineWidth', 3)
[sca, scb] = coefForEveryWavelength(E, r, eneis, N, theta, phi, 4*pi/146);
yyaxis right
plot(eneis, sca , 'r--', 'LineWidth', 2)
plot(eneis, scb , 'g--', 'LineWidth', 2)
[sca, scb] = coefForEveryWavelength(E, r, eneis, N, theta, phi, leb.w);
% yyaxis right
plot(eneis, sca , 'p-.', 'LineWidth', 1)
plot(eneis, scb , 'o-.', 'LineWidth', 1)
hold off

%% 
k = K(16);
e = E(:,:,16);
[~, A] = scattCoefficientA(e, k, r, theta, phi, leb.w);

%%
[ A(3) - A(1); 1j*(A(3)+A(1)); -sqrt(2)*A(2)]