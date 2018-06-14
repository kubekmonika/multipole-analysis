%%sprawdzamy nowe harmoniki sferyczne
wavevector = @(n, enei) 2 * pi * n / ((enei) * 1e-9);

%%
K = eneis * 0;
for i = 1 : length(eneis)
    K(i) = wavevector(1, eneis(i));
end

%%
for i = -1:1:1
    for j = -1:1:1
        [X, ~] = vsh(i, theta, phi, 10, 10);
        [~, Y] = vsh(j, theta, phi, 10, 10);
        X = squeeze( X );
        Y = squeeze( Y );
        X(isnan(X)) = 0;
        Y(isnan(Y)) = 0;
        dotX = dot(X, Y, 2);
        dotXsin = dotX .* sin(theta);
        disp({i, j, sum(dotXsin(:))})
    end
end
%%
M = vshN(1, theta, phi, 10, 10);

figure;
hold on

plot3(x, y, z, 'Ob')
quiver3(x, y, z, M(:,1), M(:,2), M(:,3))
% quiver3(x, y, z, e2(:,1), e2(:,2), e2(:,3))

% plot(p)
% quiver3(x, y, z, vshN(:,1), vshN(:,2), vshN(:,3))
% quiver3(x, y, z, vshM(:,1), vshM(:,2), vshM(:,3))

xlabel( 'x (nm)' );
ylabel( 'y (nm)' );
zlabel( 'z (nm)' );
%%
i = 1;
k = K(i);
e = E(:,:,i);
[a, A] = scattcoefficienta(1, k, e, r, theta, phi);
[b, B] = scattcoefficientb(1, k, e, r, theta, phi);

%%
