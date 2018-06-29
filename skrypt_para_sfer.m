%% Skrypt dla pary sfer
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );
%  tablica funkcji dielektrycznych
epstab = { epsconst( 1 ), epstable( 'silicon2.dat'), epstable( 'silicon2.dat' ) };

% parametry sfer
diameter = 150;
gap      =  10;

p1 = shift( trisphere( 256, diameter ), [   diameter/2.0 + gap/2.0, 0, 0] );
p2 = shift( trisphere( 256, diameter ), [ - diameter/2.0 - gap/2.0, 0, 0 ] );

p = comparticle( epstab, { p1, p2 }, [ 2, 1; 3, 1 ], 1, 2, op );

%% obliczamy wspolrzedne
n_points = 350;
R = 160;
[x, y, z, w, r, theta, phi] = getCoordinates(n_points, R);

%% Definiujemy parametry dla osrodka
% dlugosci fali
eneis = linspace(400, 720, 64);  
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
% macierz wynikowa
Ecart = zeros(length(x), 3, length(eneis));

%% Liczymy pole
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

%%
figure('pos', [50 50 700 400])
hold on
plot(eneis, C_a, '-r', 'LineWidth', 2)
plot(eneis, C_b, '-b', 'LineWidth', 2)
legend('a1', 'b1')
xlabel( 'Wavelength (nm)', 'FontSize', 13 );