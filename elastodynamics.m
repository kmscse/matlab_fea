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
