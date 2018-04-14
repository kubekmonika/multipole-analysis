function f = sphreshapefield(e, theta_n, phi_n)
% Zmienia wymiary pola wektora na takie, aby byly
% zgodne z wymiarami VSH.
% 
%   SPHRESHAPEFIELD(e, theta_n, phi_n)
% 
%   e - pole wektorowe o wymiarze (:,3)
%   theta_n - dlugosc wektora z katami theta
%   phi_n - dlugosc wektora z katami phi
% 
%   Zwraca pole jako tablice o wymiarach (phi_n x theta_n x 3), 
%   gdzie:
%   - pierwszy wymiar odpowiada wartosciom dla kolejnych katow phi
%   - drugi wymiar odpowiada wartosciom dla kolejnych katow theta
%   - trzeci wymiar to kolejno wpolrzedne sferyczne r, az, el

f = zeros(phi_n, theta_n, 3);
f(:,:,1) = reshape(e(:,1), phi_n, theta_n);
f(:,:,2) = reshape(e(:,2), phi_n, theta_n);
f(:,:,3) = reshape(e(:,3), phi_n, theta_n);
end