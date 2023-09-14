%BSM
% Parameters
T = 1; %final time
r = 0.15; %risk-free interest rate
sigma = 0.2; %volatility
K = 5; %strike price

% Define grid parameters
N = 50;    % Number of price steps
M = 5000; % Number of time steps
S_max = 10; % Maximum value of S

% Calculate step sizes
dt = T / M;
dS = S_max / N;

% Initialize price and time grids
S = 0:dS:S_max;
t = 0:dt:T;

% Initialize option value matrix
V = zeros(N+1, M+1);

% Set boundary conditions
V(:, M+1) = max(S(:) - K, 0);  
V(1, :) = 0;
V(N+1, :) = S_max - K*exp(-r*(T-t(:)));

for j = M:-1:1
    for n = 2:N
            V(n, j) = V(n, j+1) - dt*(r*V(n, j+1) - r*S(n)*(V(n+1, j+1)-V(n-1, j+1))/(2*dS) - 0.5*sigma^2*S(n)^2*(V(n+1, j+1)-2*V(n, j+1)+V(n-1, j+1))/(dS^2));
    end
end
 
% Perform backward finite difference method
%for j = M:-1:1
%    for n = 2:N
%        V(n, j) = V(n-1, j+1)*0.5*dt*(-r*n+sigma^2*n*n) + V(n, j+1)*(1-dt*(sigma^2*n*n+r)) + V(n+1, j+1)*0.5*dt*(r*n+sigma^2*n*n);
%        
        % Enforce non-negativity constraint
        %V(n, j+1) = max(V(n, j+1), 0);
%    end
%nd

% Generate the 3D graph
[S_mesh, T_mesh] = meshgrid(S, t); % Corrected to use S and t
V_mesh = V';
figure;
h = surf(S_mesh, T_mesh, V_mesh);
set(h, 'LineStyle', 'none');
xlabel('S');
ylabel('t');
zlabel('V');
title('3D Graph of V(S, t)');