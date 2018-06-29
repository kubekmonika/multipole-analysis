function [x, y, z, w, r, theta, phi] = getCoordinates(n_points, R)
% Funkcja oblicza współrzędne sfery składającej się z podanej ilości 
% punktów oraz mającej zadany promień.
% 
%   GETCOORDINATES(n_points, R)
%   
%   n_points - liczba punktów na sferze: { 6, 14, 26, 38, 50, 74, 86, 110, 
%       146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 
%       1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 },
%       jeśli zostanie podana inna wartość to zostanie zaokrąglona do
%       najbliższej możliwej wartości
%   R - promień sfery [nm]
%
%   x, y, z - współrzędne kartezjańskie
%   w - wagi dla każdego z punktów, potrzebne do całkowania
%   r - promień sfery [m]
%   thtea, phi - współrzędne sferyczne punktów

% liczba punktów
n = getNPoints(n_points);
% współrzędne kartezjańskie
leb = getLebedevSphere(n);
w = leb.w;
[az, el, ~] = cart2sph(leb.x, leb.y, leb.z);
R = ones(size(az)) * R;
[x, y, z] = sph2cart(az, el, R);
% współrzędne sferyczne
r = R(1) * 1e-9;
[~, theta, phi] = cartToSph(x, y, z);
end

function n = getNPoints(n_points)
% Funkcja zwracająca liczbę punktów potrzebną do obliczenia getLebedevSphere
points = [ 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, ...
    350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, ...
    3470, 3890, 4334, 4802, 5294, 5810 ];
if ismember(n_points, points)
    n = n_points;
else
    diff = min(abs(points - n_points));
    if ismember(n_points - diff, points)
        n = n_points - diff;
    elseif ismember(n_points + diff, points)
        n = n_points + diff;
    else
        error('Coś poszło nie tak: zła wartość n_points')
    end
end
end