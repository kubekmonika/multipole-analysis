function N = vshN(m, theta, phi, r, k)
% Liczy VSH-N dla l=1
z = k * r;
global h dh
if m >= 0
    h = Hankel(z);
    dh = Hankel(z, true);
    N = vsh(m, theta, phi, z);
else
    m = abs(m);
    h = conj(Hankel(z));
    dh = conj(Hankel(z, true));
    N = (-1)^m * factorial(1-m) / factorial(1+m) * ...
        vsh(m, theta, phi, z);
end
end

function N = vsh(m, theta, phi, z)
% Liczy VSH_N
global h dh
N = zeros(length(phi), length(theta), 3);
alpha = exp(1j * m * phi') * (h / z + dh);
N(:,:,1) = 2 * exp(1j * m * phi') .* P(m, theta) * h / z;
N(:,:,2) = alpha .* Tau(m, theta);
N(:,:,3) = 1j * alpha .* Pi(m, theta);

% if m == 0
%     alpha = ones(length(phi), 1) * exp(1j * z) / (z^3);
%     N(:,:,1) = -2 * alpha .* cos(theta) * (z + 1j);
% elseif m == 1
%     alpha = exp(1j * phi') * exp(1j * z) / (z^3);
%     N(:,:,1) = -2 * alpha .* sin(theta) * (z + 1j);
%     N(:,:,2) = alpha .* cos(theta) * (z - 1j*z^2 + 1j);
%     N(:,:,3) = 1j * alpha .* ones(1, length(theta)) * (z - 1j*z^2 + 1j);
% end
end

function p = Pi(m, theta)
% Liczy funkcje Pi 
assert(m>=0, 'Error: m < 0 w funkcji Pi')
if m == 0
    p = zeros(1, length(theta));
elseif m == 1
    p = ones(1, length(theta));
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
    p = sin(theta);
else
    assert(false, 'Error: m > 1 w funkcji P')
end
end