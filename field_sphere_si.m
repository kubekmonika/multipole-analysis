function [E, eneis, N] = field_sphere_si(dir, pol, x, y, z)
% Oblicza pole we wspolrzednych sferycznych dla krzemowej sfery 
% o srednicy 200nm, kolejno dla roznych dlugosci fali padajacej.
% 
%   FIELD_SPHERE_SI(n_points, dir, pol)
% 
%   dir - kierunek fali padajacej
%   pol - polaryzacja fali padajacej
%   x, y, z - wspolrzedne punktow w ktorych liczymy pole

% Definiujemy parametry dla osrodka w ktorym znajduje sie sfera.
% dlugosci fali
eneis = linspace(450, 850, 40);  
% relacja dyspersyjna, ogólniej n może być liczbą, 
% ale dla ośrodków z absorpcją lub nieizotropowych metoda nie działa
N = ones(length(eneis), 1);  % wspolczynnik zalamania

% % kierunek fali
% dir = [ 0, 0, 1 ];
% % polaryzacja fali
% pol = [ 0, 1, 0 ];

% TWORZYMY OBIEKTY
global op p bem exc
op = bemoptions( 'sim', 'ret', 'interp', 'curv', 'refine', 2 );
epstab = { epsconst( 1 ), epstable( 'silicon2.dat' ) };
%  nanosfera
p = trisphere( 144, 200 );
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );
% solver
bem = bemsolver( p, op );
%kierunek fali
exc = planewave( dir, pol, op );
% macierz wynikowa
E = zeros(length(x), 3, length(eneis));

% LICZYMY POLE
multiWaitbar( 'Obliczanie pola', 0, 'Color', 'g', 'CanCancel', 'on' );
n_eneis = length(eneis);
for i = 1 : n_eneis
    E(:,:,i) = calculatefield(eneis(i), x, y, z);
    multiWaitbar( 'Obliczanie pola', i / n_eneis );
end
multiWaitbar('Obliczanie pola', 'Close');

end

function e = calculatefield(enei, x, y, z)
% Liczy pole dla zadanej dlugosci fali enei
% Zwraca pole we współrzędnych kartezjanskich
global op p bem exc
sig = bem \ exc( p, enei );
emesh = meshfield( p, x, y, z, op, 'mindist', 2, 'nmax', 1200 , 'waitbar', 0);
e = emesh( sig ) ; 
end