function [E, R, eneis, n, theta, phi] = field_sphere_si(az_n, el_n)
% Oblicza pole dla krzemowej sfery o srednicy 200nm, 
% kolejno dla roznych dlugosci fali padajacej 
% enei=(350, 950)nm co 10 nm.
% 
%   FIELD_SPHERE_SI(az_n, el_n)
% 
%   az_n - liczba punktow dla wspolrzednej azymutalnej
%   el_n - liczba punktow dla wspolrzednej zenitalnej

% IMPORTUJEMY DLUGOSCI FALI
% filename = '/home/monika/HDD/MATLAB/analiza_multipolowa/krzem_n.txt';
% delimiter = '\t';
% formatSpec = '%f%f%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
% fclose(fileID);
% eneis = dataArray{:, 1}' * 1000; %[nm]
% n = dataArray{:, 2}';
% clearvars filename delimiter formatSpec fileID dataArray ans;

% definiujemy parametry dla osrodka w ktorym znajduje sie sfera
eneis = linspace(450, 950, 100);
n = ones(1, length(eneis));

% TWORZYMY OBIEKTY
global op p bem exc x y z
op = bemoptions( 'sim', 'ret', 'interp', 'curv', 'refine', 2 );
epstab = { epsconst( 1 ), epstable( 'silicon2.dat' ) };
%  nanosfera
p = trisphere( 144, 200 );
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );
% promien otaczajcej sfery
R = 110;
% solver
bem = bemsolver( p, op );
%kierunek fali
exc = planewave( [ 1, 0, 0 ], [ 0, 0, 1 ], op );
% wspolrzedne w ktorych liczymy pole
[x, y, z, ~, ~, ~, theta, phi] = sferawspl(R, az_n, el_n);
% macierz wynikowa
E = zeros(az_n*el_n, 3, length(eneis));

% LICZYMY POLE
multiWaitbar( 'Obliczanie pola', 0, 'Color', 'g', 'CanCancel', 'on' );
N = length(eneis);
for i = 1 : N
    E(:,:,i) = calculatefield(eneis(i));
    multiWaitbar( 'Obliczanie pola', i / N );
end
multiWaitbar('Obliczanie pola', 'Close');
end

function e = calculatefield(enei)
% Liczy pole dla zadanej dlugosci fali enei
global op p bem exc x y z
sig = bem \ exc( p, enei );
emesh = meshfield( p, x, y, z, op, 'mindist', 1, 'nmax', 5000 , 'waitbar', 0);
e = emesh( sig ) ; 
end