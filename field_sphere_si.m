function [E, R, eneis, N, x, y, z] = field_sphere_si(n_points, dir, pol)
% Oblicza pole we wspolrzednych sferycznych dla krzemowej sfery 
% o srednicy 200nm, kolejno dla roznych dlugosci fali padajacej.
% 
%   FIELD_SPHERE_SI(n_points, dir, pol)
% 
%   n_points - liczba punktow do triangulacji; mozliwe wartosci:
%       32 60 144 169 225 256 289 324 361 400 441 484 529 576 625
%       676 729 784 841 900 961 1024 1225 1444
%   dir - kierunek fali padajacej
%   pol - polaryzacja fali padajacej

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
% promien otaczającej sfery
R = 105;
% solver
bem = bemsolver( p, op );
%kierunek fali
exc = planewave( dir, pol, op );
% wspolrzedne w ktorych liczymy pole
sphere = trisphere(n_points, 2*R);
x = sphere.verts(:,1);
y = sphere.verts(:,2);
z = sphere.verts(:,3);
% macierz wynikowa
E = zeros(n_points, 3, length(eneis));

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