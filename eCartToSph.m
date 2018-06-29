function E = eCartToSph(Ecart, theta, phi)
% Transformacja pola elektrycznego ze wspolrzednych kartezjanskich do
% sferycznych.
% 
% Na podstawie: 
% https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates#Spherical_coordinate_system
% 

cart2sphvec = @(theta, phi) [sin(theta)*cos(phi), sin(theta)*sin(phi),  cos(theta);...
                             cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta);...
                                       -sin(phi),            cos(phi),          0];

E = 0 * Ecart;
for i = 1 : size(E, 3)
    for j = 1 : size(E,1)
        E(j,:,i) = cart2sphvec(theta(j), phi(j)) * Ecart(j,:,i)';
    end
end
E(isnan(E)) = 0;
end