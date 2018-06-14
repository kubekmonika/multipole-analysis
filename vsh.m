function [V1, V2] = vsh(m, theta, phi, r, k)
% Liczy wektorowe harmoniki sferyczne (VSH) pierwszego stopnia
% na sferze o promieniu r i w punktach o wpółrzędnych określonych 
% przez kąty theta i phi.
% 
%   VSH(m, theta, phi, r, k)
%   
%   m - rząd, |m|<=1
%   theta - współrzędne zenitalne, 0 <= theta <= pi
%   phi - współrzędne azymutalne, 0 <= phi <= 2pi
%   r - promień sfery
%   k - wartość wektora falowego
% 
% Zwraca obie VSH w kolejności: [N, M].

z = k * r;
global h dh
if m >= 0
    h = Hankel(z);
    dh = Hankel(z, true);
    V1 = vshN(m, theta, phi, z);
    V2 = vshM(m, theta, phi);
else
    m = abs(m);
    h = conj(Hankel(z));
    dh = conj(Hankel(z, true));
    V1 = (-1)^m * factorial(1-m) / factorial(1+m) * ...
        vshN(m, theta, phi, z);
    V2 = (-1)^m * factorial(1-m) / factorial(1+m) * ...
        vshM(m, theta, phi);
end
end

function M = vshM(m, theta, phi)
% Liczy pierwszą harmonikę M
global h
M = zeros(length(phi), 3);
M(:,2) = 1j * exp(1j * m * phi) .* Pi(m, theta);
M(:,3) = - exp(1j * m * phi) .* Tau(m, theta);
M = M * h;
end

function N = vshN(m, theta, phi, z)
% Liczy drugą harmonikę N
global h dh
N = zeros(length(phi), 3);
alpha = exp(1j * m * phi) * (h / z + dh);
N(:,1) = 2 * exp(1j * m * phi) .* P(m, theta) * h / z;
N(:,2) = alpha .* Tau(m, theta);
N(:,3) = 1j * alpha .* Pi(m, theta);
end

function p = Pi(m, theta)
% Liczy funkcje Pi 
assert(m>=0, 'Error: m < 0 w funkcji Pi')
if m == 0
    p = zeros(length(theta), 1);
elseif m == 1
    p = -ones(length(theta), 1);
else
    assert(false, 'Error: m > 1 w funkcji Pi')
end
end

function t = Tau(m, theta)
% Liczy funkcje Tau
assert(m>=0, 'Error: m < 0 w funkcji Tau')
if m == 0
    t = -sin(theta);
elseif m == 1
    t = -cos(theta);
else
    assert(false, 'Error: m > 1 w funkcji Tau')
end
end

function h = Hankel(z, derivative)
% Liczy funkcje Hankela pierwszego rodzaju lub jej pochodna dla n=1
if nargin == 1
    h = - exp(1j * z) * (z + 1j) / z^2;
elseif (nargin == 2) && (derivative == true)
    h = exp(1j * z) * (2 * z - 1j * z^2 + 2j) / z^3;
end
end 

function p = P(m, theta)
% Liczy funkcje Legendre'a pierwszego rodzaju
if m == 0
    p = cos(theta);
elseif m == 1
    p = - abs(sin(theta));
else
    assert(false, 'Error: m > 1 w funkcji P')
end
end