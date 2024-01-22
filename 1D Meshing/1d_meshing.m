% MATLAB codes for Finite Element Analysis
% 1 Dimensional Meshing Program
% Author: Kaung Myat San
% Location: England, the United Kingdom
% Time: 8:16 P.M London Time 22/01/2024

% Generate 1-D mesh
% Input: start_x, end_x: x-coordinates of the starting and ending nodes
% Input: n_elements: number of elements
% Output: saved data files: nodes.dat and elements.dat

function CreateMesh1D(start_x, end_x, n_elements)
    dl = (end_x - start_x)/n_elements;  % element length
    nodes = zeros(n_elements+1, 2); % empty nodes matrix
    nodes(:, 1) = 1:n_elements+1;   % 1st column: global node index
    nodes(:, 2) = start_x: dl:end_x;    % 2nd column: x-coordinates
    % save to nodes.dat
    dlmwrite('nodes.dat', nodes, 'delimiter', '\t', 'precision', '%.6f');
    elements = zeros(n_elements, 3);    % empty elements matrix
    elements (:,1) = 1:n_elements;      % 1st column: element index
    elements (:,2) = 1:n_elements;      % global node index of 1st node
    elements (:,3) = 2:n_elements+1;    % global node index of 2nd node
    % save to elements.dat
    dlmwrite('elements.dat', elements, 'delimiter', '\t', 'precision', '%d');
end

% Compute the value of the 1-D linear shape function and its first derivative at a set of input points.
% Input: satrt_x, end_x: starting and ending points of the element
% Input: x_vector: a list of input positions where the shape function and
% its derivatice should be evaluated.
% Output: N, Nx: shape function and itsx-derivative are stored in N and Nx,
% respectively, in the format shown as follows
% N = [N1(x1) N1(x20 N1(x3) ...
%    N2(x1) N2(x2) N2(x3) ...]
% Nx = [N'1(x1) N'1(x2) N'1(x3) ...
%       N'2(x1) N'2(x2) N'2(x3) ...]
function [N, Nx] = CompElementShapeLinear1D(start_x, end_x, x_vector)
    n = size(x_vector, 1);      % obtain the size of x_vector
    [N, Nx] = deal(zeros(2, n));    % setup empty N and Nx
    for i=1:n   % loop over each point x
        N(1, i) = (x_vector(i) - end_x)/(start_x - end_x); %N1(x)
        N(2, i) = (x_vector(i) - start_x)/(end_x - start_x); % N2(x)
        Nx(1,i) = 1/(start_x - end_x);  % N'1(x)
        Nx(2, i) = 1/(end_x - start_x); % N'2(x)
    end
end

% Assemble element stiffness matrix k and force vector f into the global
% matrix K and vector F
% Input: K, F: (empty) global stiffness matrix and force vector
% Input: eid: element index in the mesh
% Input: elements: the element matrix loaded from "elements.dat"
% Input: k,f: local stiffness matrix and force vector
% Outpyt: K, F: Assembled global stiffness matrix and force vector
function [K, F] = AssembleGlobalLinear1D(K, F, eid, elements, k, f)
    gids = [elements(eid, 2) elements(eid, 3)]';    % global indices of the 1st and 2nd nodes of element eid
    for m = 1:2                             % loop over rows
        F(gids(m)) = F(gids(m))+f(m);       % f(m) added onto global F
        for n=1:2                           % loop over columns
            K(gids(m), gids(n)) = K(gids(m), gids(n)+k(m, n);   % k(m,n) added onto K
        end
    end
end

% Function to compute the value of the 1-D linear isoparametric shape functions and their first derivatives. 
% The input vector xi_vector contains a list of points [xi] in the master element (between -1 and 1) where the shape functions should be evaluated. 
% The values of the shape functions and their x-derivative are stored in N and Nx, respectively, in the format shown as follows.
% N = [N1(xi1) N1(xi2) N1(xi3) ... 
%       N2(xi1) N2(xi2) N2(xi3) ... ]
% Nx = [N'1(xi1) N'1(xi2) N'1(xi3) ... 
%       N'2(xi1) N'2(xi2) N'2(xi3) ... ]
function [N, Nx] = CompShapeLinear1D(N,Nx,xi_vector)
    n = size(xi_vector, 1);
    for i=1:n
        N(1,i) = (1-xi_vector(i))/2.0;      % N1(xi)
        N(2,i) = (1+xi_vector(i))/2.0;      % N2(xi)
        Nx(1:2, i) = [-0.5 0.5]';           % Nx(xi)
    end
end

% Compute the Jacobian of an isoparametric element mapping
% Input: x_vector: vector containing x-coordinates of nodes
% Input: Nx, shape function derivatives, from "CompShapeLinear1D"
% Ouput: J, Jacobian vector J(xi_1, xi_2, ...)^T
function J = CompaJacobian1D(x_vector, Nx)
    J=(x_vector'*Nx)';      % compute vector of Jacobian
end