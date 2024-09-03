function [abl] = precalculate_abl(nX, nY, nT, boundary_factor)

if nargin < 4
    boundary_factor = 1/5; % Default value if not provided
end

% Ensure boundary factor is within a valid range
if boundary_factor <= 0 || boundary_factor >= min([nX, nY, nT])
    error('boundary_factor must be greater than 0 and less than the minimum of nX, nY, and nT');
end

disp(['Precalculating absorbing boundary layer...'])
abl_tmp = ablvec(nX, max(round(nX * boundary_factor), 1)) * ablvec(nY, max(round(nY * boundary_factor), 1))';
abl_vec = ablvec(nT, max(round(nT * boundary_factor), 1));
abl = zeros(nX, nY, nT); % Pre-allocate abl array for efficiency
for k = 1:nT
    abl(:, :, k) = abl_tmp * abl_vec(k);
end

abl = single(abl); % Convert abl to single precision
disp(['done.'])
end



