
% Paramètres
a = 0; b = 1; c = 0; d = 1;  % Rectangle [0,1]x[0,1]
Lx = 10; Ly = 10;  % Nombre de subdivisions
alpha = 1; beta = 1;  % Paramètres de la fonctionnelle
epsilon = 1e-6; maxIter = 1000;

% Génération de la triangulation
[P, C, D] = maillagerectangle(a, b, c, d, Lx, Ly);

% Résolution du problème de minimisation
U = descente(P, C, D, @g, alpha, beta, epsilon, maxIter);

% Visualisation de la surface minimale
trisurf(C, P(:, 1), P(:, 2), U);
xlabel('X');
ylabel('Y');
zlabel('U');
title('Surface minimale');