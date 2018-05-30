function N = vshN(m, theta, phi, r, k)
% Liczy VSH-N dla l=1
z = k * r;
if m >= 0
    N = vsh(m, theta, phi, z);
else
    N = (-1)^abs(m) * conj(vsh(abs(m), theta, phi, z));
end
end

function N = vsh(m, theta, phi, z)
N = zeros(length(phi), length(theta), 3);
if m == 0
    alpha = ones(length(phi), 1) * exp(1j * z) / (z^3);
    N(:,:,1) = -2 * alpha .* cos(theta) * (z + 1j);
elseif m == 1
    alpha = exp(1j * phi') * exp(1j * z) / (z^3);
    N(:,:,1) = -2 * alpha .* sin(theta) * (z + 1j);
    N(:,:,2) = alpha .* cos(theta) * (z - 1j*z^2 + 1j);
    N(:,:,3) = 1j * alpha .* ones(1, length(theta)) * (z - 1j*z^2 + 1j);
end
end