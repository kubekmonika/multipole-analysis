%% Skrypt w ktorym obliczamy wspolczynniki rozproszenia.

%% Dane
% R - promien sfery
% enei - dlugosc fali
% e - pole
% theta, phi - katy sferyczne
% az, el - wspolrzedne sferyczne

%%
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
e = Ekart(:, :, 1);
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