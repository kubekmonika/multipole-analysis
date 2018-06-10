function M = vshM(m, theta, phi, r, k)
% Liczy VSH-M dla l=1
z = k * r;
if m >= 0
    M = vsh(m, theta, phi) * Hankel(z);
else
    m = abs(m);
    M = (-1)^m * factorial(1-m) / factorial(1+m) * conj(Hankel(z)) *...
        vsh(m, theta, phi);
end
end

function M = vsh(m, theta, phi)
% Liczy składowe VSH_M
M = zeros(length(phi), 3);
M(:,2) = 1j * exp(1j * m * phi) .* Pi(m, theta);
M(:,3) = exp(1j * m * phi) .* Tau(m, theta);

% if m == 0
%     M(:,:,3) = - ones(length(phi), 1) .* sin(theta) * exp(1j * z) ...
%         * (z + 1j) / z^2;
% elseif m == 1
%     alpha = exp(1j * phi') * exp(1j * z) * (z + 1j) / z^2;
%     M(:,:,2) = -1j * alpha .* ones(1, length(theta));
%     M(:,:,3) = alpha .* cos(theta);
% end
end

function p = Pi(m, theta)
% Liczy funkcje Pi 
assert(m>=0, 'Error: m < 0 w funkcji Pi')
if m == 0
    p = zeros(length(theta), 1);
elseif m == 1
    p = ones(length(theta), 1);
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