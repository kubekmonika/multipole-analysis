function f = sphreshapefield(e, n_theta, n_phi)
% Zmienia wymiary pola wektora na takie, aby byly
% zgodne z wymiarami VSH.
% 
%   SPHRESHAPEFIELD(e, n_theta, n_phi)
% 
%   e - pole wektorowe o wymiarze (:,3), [r, theta, phi]
%   n_theta - dlugosc wektora z katami theta
%   n_phi - dlugosc wektora z katami phi
% 
%   Zwraca pole jako tablice o wymiarach (n_phi x n_theta x 3), 
%   gdzie:
%   - pierwszy wymiar odpowiada wartosciom dla kolejnych katow phi
%   - drugi wymiar odpowiada wartosciom dla kolejnych katow theta
%   - trzeci wymiar to kolejno wpolrzedne sferyczne r, el, az

f = zeros(n_phi, n_theta, 3);
f(:,:,1) = reshape(e(:,1), n_phi, n_theta);
f(:,:,2) = reshape(e(:,2), n_phi, n_theta);
f(:,:,3) = reshape(e(:,3), n_phi, n_theta);
end