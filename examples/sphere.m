%% Dane
% R - promien sfery
% eneis - dlugosc fali
% Ecart - pole we wspolrzednych kartezjanskich
% E - pole we wspolrzednych sferycznych
% theta, phi - katy sferyczne
% az, el - wspolrzedne sferyczne
% N - wspolczynnik zalamania w osrodku otaczajacym czastke

%% Sfera w próżni
op = bemoptions( 'sim', 'ret', 'interp', 'curv', 'refine', 2 );
%  tablica funkcji dielektrycznych
epstab = { epsconst( 1 ), epstable( 'silicon2.dat' ) };

%  nanosfera
p = trisphere( 144, 200 );
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );

%% obliczamy wspolrzedne
n_points = 170;
R = 105;
[x, y, z, w, r, theta, phi] = getCoordinates(n_points, R);

%% Definiujemy parametry dla osrodka
% dlugosci fali
eneis = linspace(450, 850, 40); 
% relacja dyspersyjna
N = ones(length(eneis), 1);

%% Fala
dir = [1, 0, 0];
pol = [0, 0, 1];

%% Solwer
% solver
bem = bemsolver( p, op );
% kierunek fali
exc = planewave( dir, pol, op );

%% Liczymy pole

% macierz wynikowa
Ecart = zeros(length(x), 3, length(eneis));

multiWaitbar( 'Obliczanie pola', 0, 'Color', 'g', 'CanCancel', 'on' );
n_eneis = length(eneis);
for i = 1 : n_eneis
    Ecart(:,:,i) = calculateField(eneis(i), x, y, z, op, p, bem, exc);
    multiWaitbar( 'Obliczanie pola', i / n_eneis );
end
multiWaitbar('Obliczanie pola', 'Close');

E = eCartToSph(Ecart, theta, phi);
dotE = squeeze( sum( squeeze( dot(E, E, 2) ) .* w, 1 ) );

%% Liczymy wspolczynniki rozproszenia
[C_a, C_b] = coefForEveryWavelength(E, r, eneis, N, theta, phi, w);

%% Wykres
figure
hold on
subplot(2,1,1)
% yyaxis left
plot(eneis, C_a, '-', eneis, C_b, '--')
% yyaxis right
% plot(eneis, C_b, '--')
legend('a1', 'b1')

subplot(2,1,2)
yyaxis left
C_ab = C_a + C_b;
plot(eneis, C_ab, '--')
yyaxis right
plot(eneis, dotE, '-', 'LineWidth', 2)
legend('a+b', 'dot(E,E)')
xlabel( 'Wavelength (nm)', 'FontSize', 11 );

%%
figure('pos', [50 50 700 400])
hold on
plot(eneis, C_a, '-r', 'LineWidth', 2)
plot(eneis, C_b, '-b', 'LineWidth', 2)
legend('a1', 'b1')
xlabel( 'Wavelength (nm)', 'FontSize', 13 );

%% Sfera otaczająca badany układ

figure
hold on
plot(p)
plot3( x, y, z, 'b.' );
xlabel( 'x (nm)' );
ylabel( 'y (nm)' );
zlabel( 'z (nm)' );
