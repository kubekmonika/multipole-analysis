function [E, R, eneis, N, theta, phi] = field_sphere_si(n_theta, n_phi, dir, pol)
% Oblicza pole we wspolrzednych sferycznych dla krzemowej sfery 
% o srednicy 200nm, kolejno dla roznych dlugosci fali padajacej.
% 
%   FIELD_SPHERE_SI(n_theta, n_phi)
% 
%   n_phi - liczba punktow dla wspolrzednej azymutalnej
%   n_theta - liczba punktow dla wspolrzednej zenitalnej

% Definiujemy parametry dla osrodka w ktorym znajduje sie sfera.
% dlugosci fali
eneis = linspace(450, 850, 80);  
% relacja dyspersyjna, ogólniej n może być liczbą, 
% ale dla ośrodków z absorpcją lub nieizotropowych metoda nie działa
N = ones(1, length(eneis));  % wspolczynnik zalamania

% % kierunek fali
% dir = [ 0, 0, 1 ];
% % polaryzacja fali
% pol = [ 0, 1, 0 ];

% TWORZYMY OBIEKTY
global op p bem exc x y z
op = bemoptions( 'sim', 'ret', 'interp', 'curv', 'refine', 2 );
epstab = { epsconst( 1 ), epstable( 'silicon2.dat' ) };
%  nanosfera
p = trisphere( 144, 200 );
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );
% promien otaczajcej sfery
R = 105;
% solver
bem = bemsolver( p, op );
%kierunek fali
exc = planewave( dir, pol, op );
% wspolrzedne w ktorych liczymy pole
[x, y, z, ~, ~, ~, theta, phi] = sferawspl(R, n_theta, n_phi);
% macierz wynikowa
E = zeros(n_phi*n_theta, 3, length(eneis));

% LICZYMY POLE
multiWaitbar( 'Obliczanie pola', 0, 'Color', 'g', 'CanCancel', 'on' );
n_eneis = length(eneis);
for i = 1 : n_eneis
    E(:,:,i) = calculatefield(eneis(i));
    multiWaitbar( 'Obliczanie pola', i / n_eneis );
end
multiWaitbar('Obliczanie pola', 'Close');

end

function e = calculatefield(enei)
% Liczy pole dla zadanej dlugosci fali enei
% Zwraca pole we współrzędnych kartezjanskich
global op p bem exc x y z
sig = bem \ exc( p, enei );
emesh = meshfield( p, x, y, z, op, 'mindist', 1, 'nmax', 1200 , 'waitbar', 0);
e = emesh( sig ) ; 
end