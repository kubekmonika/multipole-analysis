%% Dane
% R - promien sfery
% eneis - dlugosc fali
% Ekart - pole we wspolrzednych kartezjanskich
% E - pole we wspolrzednych sferycznych
% theta, phi - katy sferyczne
% az, el - wspolrzedne sferyczne
% N - wspolczynnik zalamania w osrodku otaczajacym czastke

%% obliczamy pole dla calego przekroju dlugosci fali
dir = [0, 0, 1];
pol = [1, 0, 0];
[Ekart, R, eneis, N, theta, phi] = field_sphere_si(51, 41, dir, pol);
r = R * 1e-9;
theta = theta + pi/2;
save('41x31_zx.mat', 'Ekart', 'R', 'r', 'eneis', 'N', 'theta', 'phi');
clear all

dir = [0, 0, 1];
pol = [0, 1, 0];
[Ekart, R, eneis, N, theta, phi] = field_sphere_si(51, 41, dir, pol);
r = R * 1e-9;
theta = theta + pi/2;
save('41x31_zy.mat', 'Ekart', 'R', 'r', 'eneis', 'N', 'theta', 'phi');

%%
cart2sphvec = @(az, el) [-sin(az),          cos(az),         0;...
                         -sin(el)*cos(az), -sin(el)*sin(az), cos(el);...
                          cos(el)*cos(az),  cos(el)*sin(az), sin(el)];

[~, ~, ~, az, el, ~, ~, ~] = sferawspl(R, length(theta), length(phi));

E = 0 * Ekart;
for i = 1 : size(E, 3)
    for j = 1 : size(E,1)
        E(j,:,i) = cart2sphvec(az(j), el(j)) * Ekart(j,:,i)';
    end
end

%%
EE = zeros(length(phi), length(theta), 3, length(eneis));
for j = 1 : length(eneis)
    EE(:,:,:,j) = sphreshapefield(E(:,:,j), length(theta), length(phi));
end
dotE = zeros(1, length(eneis));
for i = 1 : length(eneis)
    temp = squeeze(EE(:,:,:,i));
    temp = dot(temp, temp, 3) .* sin(theta);
    dotE(i) = sum(temp(:));% / K(i)^2;
end

% dotEkart = dotE;
%% obliczamy wspolczynniki rozproszenia
[C_a, C_b] = coefforeverywavelength(E, r, eneis, N, theta, phi);

%% wykres
od = 1;

figure
hold on
subplot(2,1,1)
yyaxis left
plot(eneis(od:end), C_a(od:end), '-')
yyaxis right
plot(eneis(od:end), C_b(od:end), '--')
legend('a1', 'b1')

subplot(2,1,2)
yyaxis left
plot(eneis(od:end), C_a(od:end) + C_b(od:end), '--')
yyaxis right
plot(eneis(od:end), dotE(od:end), '-', 'LineWidth', 2)
legend('suma', 'dot(E,E)')

%%
%  nanosphere
p = trisphere( 144, 10 );
%  plot particle
% plot( p, 'EdgeColor', 'b' );  hold on;
%  plot centroids in black
% plot3( p.pos( :, 1 ), p.pos( :, 2 ), p.pos( :, 3 ), 'k.' );
%  plot vertices at edges
p.verts(:,1) = 0;
plot3( p.verts( :, 1 ), p.verts( :, 2), p.verts( :, 3 ), 'rs' );
%  plot vertices for curved boundary element integration
% plot3( p.verts2( :, 1 ), p.verts2( :, 2 ), p.verts2( :, 3 ), 'g.' );
xlabel( 'x (nm)' );
ylabel( 'y (nm)' );
zlabel( 'z (nm)' );