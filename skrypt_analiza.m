%% Dane
% R - promien sfery
% eneis - dlugosc fali
% Ekart - pole we wspolrzednych kartezjanskich
% E - pole we wspolrzednych sferycznych
% theta, phi - katy sferyczne
% az, el - wspolrzedne sferyczne
% N - wspolczynnik zalamania w osrodku otaczajacym czastke

%% obliczamy pole dla calego przekroju dlugosci fali
dir = [0, 0, 1];
pol = [1, 0, 0];
[Ekart, R, eneis, N, x, y, z] = field_sphere_si(60, dir, pol);

%% obliczamy wspolrzedne
r = R * 1e-9;
[phi, theta] = cart2sph(x,y,z);
theta = pi/2 - theta;

%%
% na podstawie: https://www.mathworks.com/help/phased/ref/cart2sphvec.html
cart2sphvec = @(phi, theta) [-sin(phi),          cos(phi),         0;...
                         -sin(theta)*cos(phi), -sin(theta)*sin(phi), cos(theta);...
                          cos(theta)*cos(phi),  cos(theta)*sin(phi), sin(theta)];

E = 0 * Ekart;
for i = 1 : size(E, 3)
    for j = 1 : size(E,1)
        E(j,:,i) = cart2sphvec(phi(j), theta(j)) * Ekart(j,:,i)';
    end
end

%%
dotEkart = squeeze( sum( squeeze( dot(Ekart, Ekart, 2) ) .* theta, 1) );
dotE = squeeze( sum( squeeze( dot(E, E, 2) ) .* theta, 1 ) );

%% obliczamy wspolczynniki rozproszenia
% N = 1;
[C_a, C_b] = coefforeverywavelength(E, r, eneis, N, theta, phi, 'ab');

%% wykres
figure
hold on
subplot(2,1,1)
yyaxis left
plot(eneis, C_a, '-')
yyaxis right
plot(eneis, C_b, '--')
legend('a1', 'b1')

subplot(2,1,2)
yyaxis left
plot(eneis, C_a + C_b, '--')
yyaxis right
plot(eneis, dotE, '-', 'LineWidth', 2)
legend('suma', 'dot(E,E)')

%%
%  nanosphere
p = trisphere( 144, 10 );
%  plot particle
% plot( p, 'EdgeColor', 'b' );  hold on;
%  plot centroids in black
% plot3( p.pos( :, 1 ), p.pos( :, 2 ), p.pos( :, 3 ), 'k.' );
%  plot vertices at edges
p.verts(:,1) = 0;
plot3( p.verts( :, 1 ), p.verts( :, 2), p.verts( :, 3 ), 'rs' );
%  plot vertices for curved boundary element integration
% plot3( p.verts2( :, 1 ), p.verts2( :, 2 ), p.verts2( :, 3 ), 'g.' );
xlabel( 'x (nm)' );
ylabel( 'y (nm)' );
zlabel( 'z (nm)' );