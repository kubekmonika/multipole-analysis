%% Dane
% R - promien sfery
% eneis - dlugosc fali
% Ekart - pole we wspolrzednych kartezjanskich
% E - pole we wspolrzednych sferycznych
% theta, phi - katy sferyczne
% az, el - wspolrzedne sferyczne
% N - wspolczynnik zalamania w osrodku otaczajacym czastke

%% obliczamy wspolrzedne
n_points = 150;
R = 105;
[x, y, z, w, r, theta, phi] = getCoordinates(n_points, R);

%% obliczamy pole dla calego przekroju dlugosci fali
dir = [0, 0, 1];
pol = [0, 1, 0];
[Ekart, eneis, N] = field_sphere_si(dir, pol, x, y, z);

%%
% transformacja wspolrzednych
% na podstawie: https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates#Spherical_coordinate_system

cart2sphvec = @(theta, phi) [sin(theta)*cos(phi), sin(theta)*sin(phi),  cos(theta);...
                             cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta);...
                                       -sin(phi),            cos(phi),          0];

E = 0 * Ekart;
for i = 1 : size(E, 3)
    for j = 1 : size(E,1)
        E(j,:,i) = cart2sphvec(theta(j), phi(j)) * Ekart(j,:,i)';
    end
end
E(isnan(E)) = 0;

%
% dotEkart = squeeze( sum( squeeze( dot(Ekart, Ekart, 2) ) .* leb.w, 1) );
dotE = squeeze( sum( squeeze( dot(E, E, 2) ) .* leb.w, 1 ) );

% obliczamy wspolczynniki rozproszenia
% N = 1;
[C_a, C_b] = coefForEveryWavelength(E, r, eneis, N, theta, phi, leb.w);

%% wykres
od = 5;
figure
hold on
subplot(2,1,1)
% yyaxis left
plot(eneis(od:end), C_a(od:end), '-', eneis(od:end), C_b(od:end), '--')
% yyaxis right
% plot(eneis, C_b, '--')
legend('a1', 'b1')
title('dir: Z, pol: Y')

subplot(2,1,2)
yyaxis left
C_ab = C_a + C_b;
plot(eneis(od:end), C_ab(od:end), '--')
yyaxis right
plot(eneis(od:end), dotE(od:end), '-', 'LineWidth', 2)
legend('suma', 'dot(E,E)')


%%
dotE = squeeze( sum( squeeze( dot(E, E, 2) ) .* sin(theta), 1 ) ) /146 * 4 * pi;
dotEleb = squeeze( sum( squeeze( dot(E, E, 2) ).* leb.w, 1 ) );

od = 1;
figure
hold on
% yyaxis left
plot(eneis(od:end), dotE(od:end), 'k-', 'LineWidth', 2)
% yyaxis right
plot(eneis(od:end), dotEleb(od:end), 'r--', 'LineWidth', 2)

legend('E', 'Eleb')

%%
leb_xz = squeeze( sum( squeeze( dot(E, E, 2) ).* leb.w, 1 ) );

od = 1;
figure
hold on
% yyaxis left
plot(eneis(od:end), leb_xy(od:end), 'k-', 'LineWidth', 2)
% yyaxis right
plot(eneis(od:end), leb_xz(od:end), 'r--', 'LineWidth', 2)

legend('xy', 'xz')

%%
%  nanosphere
p = trisphere( 144, 10 );
%  plot particle
plot( p, 'EdgeColor', 'b' );  hold on;
%  plot centroids in black
% plot3( p.pos( :, 1 ), p.pos( :, 2 ), p.pos( :, 3 ), 'k.' );
%  plot vertices at edges
% p.verts(:,1) = 0;
plot3( p.verts( :, 1 ), p.verts( :, 2), p.verts( :, 3 ), 'ro' );
%  plot vertices for curved boundary element integration
% plot3( p.verts2( :, 1 ), p.verts2( :, 2 ), p.verts2( :, 3 ), 'g.' );
xlabel( 'x (nm)' );
ylabel( 'y (nm)' );
zlabel( 'z (nm)' );