
%% This code is developed by HUSSEIN ABDUELHALEIM HUSSEIN MUHAMMED (https://github.com/Jaguar101-jr)
%% B.ScH. University of Khartoum (Sudan) 2015, and M.Sc. from UPC-Qingdao (China) 2022.
%% A Finite-Difference program to solve the 2D Laplacian Equation (i.e. Heat Equ.) with plotting the
%% Analytical; Numerical solutions and Error.

clc;
clear;
close all;

%% define computational domain
a = 3*pi;
b = pi;
nx = 100;
ny = 50;
dx = a/nx;
dy = b/ny;
x = 0:dx:a;
y = 0:dy:b;
u = zeros(ny+1, nx+1);

%% Analytical solution
% define the CFD analytical solution loop
for  i = 1:ny+1
    for j = 1:nx+1
        u(i,j) = sin(x(j))./sin(a).*sinh(y(i))./sinh(b);
    end
end
maxu = max(u,[],'all');
nu = u/maxu;
% visualize the solution
figure()
contourf(nu,200,'linecolor','non');
title('Analytical Solution');
xlabel('x')
ylabel('y')
colormap(jet(256));
colorbar;
caxis([-1,1]);

%% solve the numerical solution
Co = 1/(2*(dx^2+dy^2));    %laplacian coeff.
U = zeros(ny+1, nx+1);

% Define the initial guess
U = zeros(ny+1, nx+1);

% Define the Boundary conditions
U (1,:) = 0;
U (ny+1,:) = sin(x)/sin(a);
U (:, 1) = 0;
U (:, nx+1) = sinh(y)/sinh(b);
% define the CFD numerical solution loop
for k = 1:1000  %%number of iterations
for  i = 2:ny
    for j = 2:nx
        U(i,j) = Co*(dx^2*(U(i+1,j) + U(i-1,j)) + dy^2* (U(i,j+1) +U(i,j-1)));
    end
end
end

maxu = max(U,[],'all');
nU = U/maxu;
% visualize the solution
figure()
contourf(nU,200,'linecolor','non');
title('Numerical Solution');
xlabel('x')
ylabel('y')
colormap(jet(256));
colorbar;
caxis([-1,1]);

%% ESTIMATE THE ERROR
E = nU - nu;
figure()
contourf(E,200,'linecolor','non');
title('Error Convergence');
xlabel('x')
ylabel('y')
colormap(jet(256));
colorbar;
caxis([-1,1]);


