%% User-friendly input script for setting input parameters
% Alternatively, you can call the makeSurface function directly as:
% makeSurface(salinity, slope, C_e, beta, g, k_s, saltname)
%
% where:

%   salinity  - array
%   slope     - scalar (double)
%   C_e       - scalar (double)
%   beta      - scalar (double)
%   g         - scalar (double)
%   k_s       - scalar (double)
%   saltname  - character array (e.g., 'NaCl')
%   path      - character array (e.g., 'path/to/directory')

salinity = [10,25,50,75,100,125,150,175,200,225];    % Salinity Values to Simulate Using SOFTBALL (ppt)
slope    = -0.0913;                                  % Liquidus Slope for Binary Salt Solution (k/ppt)
C_e      = 230;                                      % Eutectic Concentration (ppt)
beta     = 7.6*10^-4;                                % Haline Contraction Coefficient (1/ppt)
g        = 1.32;                                     % Gravity (m/s^2)
k_s      = 1.5*10^-9;                                % Salt Diffusivity in H20 (m^2/s)
saltname = 'NaCl';                                   % Name of Salt (for file naming)
path     = '/thayerfs/home/f0049bv';                 % Location of SOFTBALL's 'mushy-layer' directory

makeSurface(salinity,slope,C_e,beta,g,k_s,saltname,path); % Run makeSurface for these inputs

%%
% When run, this script will generate:
% - A saved figure to check the goodness of fit of the surface to the simulated data
% - A file named 'All_values_array.mat' containing the raw data in a loadable cell array where each row has simulations for each salinity value.
% - A file named 'Surf(saltname).mat' containing the input parameters and a function object called 'SOFTsurf'
%
% The 'SOFTsurf' object includes a function handle 'S_ice', which can be used to estimate interpolated ice salinity values.
% This is done using the syntax: SOFTsurf.S_ice(S_oc, dT_dz), where:
%   S_oc   = ocean salinity value
%   dT_dz  = thermal gradient at the freezing front

