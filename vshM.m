function M = vshM(m, theta, phi, r, k)
% Liczy VSH-M dla l=1
z = k * r;
if m >= 0
    M = vsh(m, theta, phi, z);
else
    M = (-1)^abs(m) * conj(vsh(abs(m), theta, phi, z));
end
end

function M = vsh(m, theta, phi, z)
M = zeros(length(phi), length(theta), 3);
if m == 0
    M(:,:,3) = - ones(length(phi), 1) .* sin(theta) * exp(1j * z) ...
        * (z + 1j) / z^2;
elseif m == 1
    alpha = exp(1j * phi') * exp(1j * z) * (z + 1j) / z^2;
    M(:,:,2) = -1j * alpha .* ones(1, length(theta));
    M(:,:,3) = alpha .* cos(theta);
end
end