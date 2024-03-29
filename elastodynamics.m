% MATLAB codes for Finite Element Analysis
% Elastodynamics Analysis
% Author: Kaung Myat San
% Time: 01:48 A.M Time 21/01/2024

% read the input files
filenames = {'nodes.dat', 'elements.dat' 'materials.dat', ...
    'options.dat', 'bscforce.dat', 'bcsdisp.dat'};
for i = 1:numel(filenames); load(filenames{i}); end;

% set up global material properties and constants
E = materials(1,1);
nu = materials(2,1);
rho = materials(3,1);
dampK = materials(4,1);
dampM = materials(5,1);
dimension = options(1,1);
thickness = options(2,1);
timePeriod = options(4,1);
dt = options(5,1);
probeNode = options(6,1);

% bookkeeping
n_nodes = size(nodes, 1);
n_elements = size(elements, 1);
n_bcsforce = size(bcsforce, 1);
n_csdisp = size(bcsdisp, 1);
n_timeSteps = timePeriod/dt+1;

% set up empty matrices
U = zeros(n_nodes * dimension, n_timeSteps);
V = zeros(n_nodes * dimension, n_timeSteps);
A = zeros(n_nodes * dimension, n_timeSteps);

K = CompK(nodes, elements, materials); % compute global K matrix
M = CompM(nodes, elements, rho); % compute global M matrix
C = dampM * M + dampK * K; % compute global C matrix
F = CompF(nodes, elements, thickness, bcsforce); % compute global F vector

% time marching using Newmark Scheme
beta = 0.25
gamma = 0.5;
LHS = M*(1.0/beta*dt*dt)) + C*(gamma/(beta*dt)) + K;
penalty = max(max(abs(LHS))) * 1e+6;
A(:,1) = M\F; % initialization of acceleration
for t=2:n_timeSteps
    [U,V,A] = Newmark(t, dt, beta, gamma, dimension, ...
        K,M,C,F,bcsdisp, penalty, U, V, A);
    if rem(t*100, (n_timeSteps-1)*5)==0
        fprintf('%d %%\n', floor(t*100/(n_timeSteps-1)));
    end
end

% save the displacement results in file, plot
Uout = U(probeNode * 2,:);
save -ascii -double Uout2.dat Uout
disp('Vertical dispalcement results stored in Uout.dat');
Uout(200:205)*1000;
plot(Uout, 'LineWidth',2);
axis([0 2100 -1.2e-3 0]);
xlabel('Time step');
ylabel('Dispalcement (m)');
set(gca, 'fontsize', 16);

% Compute mass matrix
function M = CompM(nodes, elements, rho)
    n_nodes = size(nodes, 1);
    n_elements = size(elements, 1);
    n_nodes_per_element = size(elemets, 2) - 1;
    M = zeros(n_nodes * 2, n_nodes*2);
    Nv = zeros(8, 2);
    [gauss_points, gauss_weights] = GetQuadGauss(2,2);
    [N, Nx, Ny] = CompNDNatPointsQuad4(gauss_points(:,1), gauss_points(:,2));
end

% for-loop: compute M matrix: loop over all the elements
for e=1:n_elements
    me = zeros(n_nodes_per_element*2, n_nodes_per_element*2);
    [element_nodes, node_id_map] = SetElementNodes9e, nodes, elements);
    % for-loop: comoute element mass matrix me
    for g=1:size(gauss_points, 1)
        J = CompJacobian2DatPoint(element_nodes, Nx(:,g), Ny(:,g));
        detJ = det(J);
        for p=1:4
            Nv(2*p-1,1) = N(p,g);
            Nv(2*p,2)=N(p,g);
        end
        me = me+Nv*Nv'*detJ*rho*gauss_weights(g);
    end
    M = AssembleGlobalMatrix(M, me, node_id_map, 2); % assemble global M
end

% Compute global force vector
function F=CompF(nodes, elements, thickness, bcsforce)
    F = zeros(size(nodes, 1)*(size(nodes, 2)-1),1);
    n_force_nodes = size(bcsforce, 1);
    % for-loop: apply point forces
    for i=1:n_force_nodes
        row = 2*bcsforce(i, 1)-1;
        F(row, 1) = F(row, 1)+ bcsforce(i, 2)/thickness;
        F(row+1, 1) = F(row+1, 1) + bcsforce(i, 3)/thickness;
    end
end

% Time integration using Newmark scheme Trapezoid rule
function [U,V,A] = Newmark(t, dt, beta, gamma, dim, K, M, C, F, bcsdisp, penalty, U, V, A)
    n_nodes = size(F,1)/dim;
    % get u,v,a vectors
    u1 = U(:, t-1);
    vel1 = V(:, t-1);
    accel1 = A(:, t-1);
    
    
    % compute the LHS matrix and RHS vector
    LHS = M*(1.0/(beta*dt*dt)) + C*(gamma/(beta*dt)) + K;
    rhsvec = ul + vel1*dt + accel1*((.5-beta)*dt*dt);
    RHS = F + M*(rhsvec*(1.0/(beta*dt*dt))) + C*(rhsvec * ...
        (gamma/(beta*dt)) - vel1 - accel1*(dt*(1.0 - gamma)));
    
    % for-loop: apply displacement BC using the penalty method
    for j=1:size(bcsdisp, 1);
        nid = bcsdisp(j, 1);
        k = bcsdisp(j, 2);
        RHS(dim*(nid-1)+k, 1) = bcsdisp(j,3)*penalty;
        LHS(dim*(nid-1)+k, dim*(nid-1)+k) = penalty;
    end
    
    U(:, t)= LHS\RHS; % solve for the displacement
    % calculated acceleration and velocity
    A(:, t) = (U(:, t) - U(:, t-1))/(beta*dt*dt) ...
        -V(:, t-1)/(beta*dt) - A(:, t-1)*(0.5-beta)/beta;
    V(:, t) = V(:, t-1) + A(:, t-1)*(1.0-gamma)*dt + A(:, t)*dt*gamma;
end
