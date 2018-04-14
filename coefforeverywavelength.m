function C = coefforeverywavelength(E, R, eneis, n, theta, phi)
%Funkcja liczaca wspolczynniki rozproszenia dla danego pola.

% funkcja liczaca wektor falowy
wavevector = @(n, enei) sqrt(4 * pi^2 * n^2 / (enei * 0.1^9)^2);

% obliczamy wspolczynniki rozproszenia
N = length(n);
C = zeros(1, N);

% pasek czasu
if isempty(gcp('nocreate'))
    parpool('local', 4);
end
dq = parallel.pool.DataQueue;
wb = waitbar(0, 'Calculating...');
wb.UserData = [0 N];
afterEach(dq, @(varargin) iIncrementWaitbar(wb));

% petla
parfor i = 1: N
    k = wavevector(n(i), eneis(i));
    e = E(:,:,i);
    C(i) = scattcoefficienta(1, k, e, R, theta, phi);
end
close(wb);
end

function iIncrementWaitbar(wb)
% Pasek czasu dla parfor
ud = wb.UserData;
ud(1) = ud(1) + 1;
waitbar(ud(1) / ud(2), wb);
wb.UserData = ud;
end