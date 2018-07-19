%% Data
% R - radius
% eneis - wavelength
% Ecart - elecric field in cartesian coordinates
% E - electric field in spherical coordinates
% theta, phi - spherical angles
% az, el - spherical coordinates
% N - refractive index

%% Pair of spheres
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );
% dielectric constants table
epstab = { epsconst( 1 ), epstable( 'silicon2.dat'), epstable( 'silicon2.dat' ) };

% parameters
diameter = 150;
gap      =  10;

p1 = shift( trisphere( 256, diameter ), [   diameter/2.0 + gap/2.0, 0, 0] );
p2 = shift( trisphere( 256, diameter ), [ - diameter/2.0 - gap/2.0, 0, 0 ] );

p = comparticle( epstab, { p1, p2 }, [ 2, 1; 3, 1 ], 1, 2, op );

%% calculate coordinates
n_points = 350;
R = 160;
[x, y, z, w, r, theta, phi] = getCoordinates(n_points, R);

%% Medium parameters
% wavelength
eneis = linspace(450, 720, 54);  
% refractive index - vacuume
N = ones(length(eneis), 1);

%% Wave parameters
% direction
dir = [1, 0, 0];
% polarization
pol = [0, 0, 1];
% wave
exc = planewave( dir, pol, op );

%% Solver
% solver
bem = bemsolver( p, op );

%% Calculation

% results matrix
Ecart = zeros(length(x), 3, length(eneis));

multiWaitbar( 'Calculating field', 0, 'Color', 'g', 'CanCancel', 'on' );
n_eneis = length(eneis);
for i = 1 : n_eneis
    Ecart(:,:,i) = calculateField(eneis(i), x, y, z, op, p, bem, exc);
    multiWaitbar( 'Calculating field', i / n_eneis );
end
multiWaitbar('Calculating field', 'Close');

E = eCartToSph(Ecart, theta, phi);
dotE = squeeze( sum( squeeze( dot(E, E, 2) ) .* w, 1 ) );

%% Calculating scattering coefficients
[C_a, C_b] = coefForEveryWavelength(E, r, eneis, N, theta, phi, w);

%% Figure 1
figure
hold on
subplot(2,1,1)
plot(eneis, C_a, '-', eneis, C_b, '--')
legend('a1', 'b1')

subplot(2,1,2)
yyaxis left
C_ab = C_a + C_b;
plot(eneis, C_ab, '--')
yyaxis right
plot(eneis, dotE, '-', 'LineWidth', 2)
legend('a+b', 'dot(E,E)')
xlabel( 'Wavelength (nm)', 'FontSize', 11 );

%% Figure 2
figure('pos', [50 50 700 400])
hold on
plot(eneis, C_a, '-r', 'LineWidth', 2)
plot(eneis, C_b, '-b', 'LineWidth', 2)
legend('a1', 'b1')
xlabel( 'Wavelength (nm)', 'FontSize', 13 );

%% Sphere around the structure

figure
hold on
plot(p)
plot3( x, y, z, 'b.' );

xlabel( 'x (nm)' );
ylabel( 'y (nm)' );
zlabel( 'z (nm)' );