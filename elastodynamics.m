% MATLAB codes for Finite Element Analysis
% Elastodynamics Analysis
% Author: Kaung Myat San
% Location: England, the United Kingdom
% Time: 01:48 A.M London Time 21/01/2024
% next 3 lines: read the input files
filenames = {'nodes.dat', 'elements.dat' 'materials.dat', ...
    'options.dat', 'bscforce.dat', 'bcsdisp.dat'};
for i = 1:numel(filenames); load(filenames{i}); end;

% next 9 lines: set up global material properties and constants
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

% next 5 lines: bookkeeping
n_nodes = size(nodes, 1);
n_elements = size(elements, 1);
n_bcsforce = size(bcsforce, 1);
n_csdisp = size(bcsdisp, 1);
n_timeSteps = timePeriod/dt+1;

% next 3 ines: set up empty matrices
U = zeros(n_nodes * dimension, n_timeSteps);
V = zeros(n_nodes * dimension, n_timeSteps);
A = zeros(n_nodes * dimension, n_timeSteps);

K = CompK(nodes, elements, materials); % compute global K matrix
M = CompM(nodes, elements, rho); % compute global M matrix
C = dampM * M + dampK * K; % compute global C matrix
F = CompF(nodes, elements, thickness, bcsforce); % compute global F vector

% next 12 lines: time marching using Newmark Scheme
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

% next 9 lines: save the displacement results in file, plot
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