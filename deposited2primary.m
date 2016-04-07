function [primary_spectrum,resnorm,residual,exitflag,output,lambda] = deposited2primary(drm,deposited_spectrum)

% The function requires Optimization and Global Optomization Toolboxes.

% Primary spectrum reconstruction from the deposited one using
% an appropriate detector response matrix (DRM) can be considered as a sparse 
% linear least-squares optimization problem with bounds. 
% Primary spectrum is nonnegative. Deposited spectrum is in primary spectrum units
% per unit area (cm^2)
%
%

if nargin ~=2,
    error('[primary_spectrum,resnorm,residual,exitflag,output,lambda] = deposited2primary(drm,deposited_spectrum)');
end

if ~ ((size(drm,1) == size(drm,2)) && (size(drm,1) == numel(deposited_spectrum)))
    error('Invalid dimensions. DRM must be NxN matrix, where N - number of elements in deposited spectrum. ');
end


% Allocate space
primary_spectrum = zeros(size(deposited_spectrum));

% Lower bounds, specified as a vector and represents
% the lower bounds elementwise
lb = zeros(size(primary_spectrum));   

% Solves the linear system drm*primary_spectrum = deposited_spectrum 
% in the least-squares sense so that the solution is always in the range
% lb < primary_spectrum
options = optimoptions('lsqlin','MaxIterations',100);
[primary_spectrum,resnorm,...
    residual,exitflag,output,lambda] = lsqlin(drm,deposited_spectrum,[],[],[],[],lb,[],[],options);


