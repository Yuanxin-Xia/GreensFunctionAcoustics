function P = Greensfunction(nxx,nyy,nzz,f_range,x,y,z,x0,y0,z0,lx,ly,lz,...
                            beta_x_0,beta_x_L,beta_y_0,beta_y_L,beta_z_0,beta_z_L,c)
tic

    % Determine number of frequency points
    f_len = length(f_range);

    % Ensure the boundary condition parameters are vectors of length f_len
    if isscalar(beta_x_0)
        beta_x_0 = repmat(beta_x_0, f_len, 1);
    else
        beta_x_0 = beta_x_0(:);
    end
    if isscalar(beta_x_L)
        beta_x_L = repmat(beta_x_L, f_len, 1);
    else
        beta_x_L = beta_x_L(:);
    end
    if isscalar(beta_y_0)
        beta_y_0 = repmat(beta_y_0, f_len, 1);
    else
        beta_y_0 = beta_y_0(:);
    end
    if isscalar(beta_y_L)
        beta_y_L = repmat(beta_y_L, f_len, 1);
    else
        beta_y_L = beta_y_L(:);
    end
    if isscalar(beta_z_0)
        beta_z_0 = repmat(beta_z_0, f_len, 1);
    else
        beta_z_0 = beta_z_0(:);
    end
    if isscalar(beta_z_L)
        beta_z_L = repmat(beta_z_L, f_len, 1);
    else
        beta_z_L = beta_z_L(:);
    end

    % Precompute Mode-Related Quantities
    V  = lx * ly * lz;  % Room volume

    % Generate a 3D grid of mode indices (nx, ny, nz)
    [nx, ny, nz] = ndgrid(0:nxx, 0:nyy, 0:nzz);
    
    % Compute Neumann factors: 1 if index is zero, otherwise 2
    neumann_x = ones(size(nx)); neumann_x(nx > 0) = 2;
    neumann_y = ones(size(ny)); neumann_y(ny > 0) = 2;
    neumann_z = ones(size(nz)); neumann_z(nz > 0) = 2;
    
    % Precompute scaling factors for the boundary loss terms
    tau_x = neumann_x / lx;
    tau_y = neumann_y / ly;
    tau_z = neumann_z / lz;
    
    % Compute the Lambda factor (reciprocal of the product of Neumann factors)
    Lambda = 1 ./ (neumann_x .* neumann_y .* neumann_z);
    
    % Compute the squared eigenvalues (Kn^2) for the modes
    Kn2 = (pi * nx / lx).^2 + (pi * ny / ly).^2 + (pi * nz / lz).^2;
    
    % Compute the mode shapes at the source and receiver locations (frequency-independent)
    source_mode = cos(pi * nx * x0 / lx) .* cos(pi * ny * y0 / ly) .* cos(pi * nz * z0 / lz);
    receiver_mode = cos(pi * nx * x / lx) .* cos(pi * ny * y / ly) .* cos(pi * nz * z / lz);
    
    % Flatten the 3D arrays into 1D vectors for use in the parallel loop
    Lambda_flat        = Lambda(:);
    Kn2_flat           = Kn2(:);
    source_mode_flat   = source_mode(:);
    receiver_mode_flat = receiver_mode(:);
    
    % Preallocate output result vector
    P = zeros(1, f_len);
    
    % Open parallel pool if not already running
    if isempty(gcp('nocreate'))
        parpool;
    end
    
    % Parallel Computation over frequency
    parfor i = 1:f_len
        f = f_range(i);
        k = 2 * pi * f / c;  % Compute the wavenumber
        
        % Compute the tau factor for the current frequency using frequency-dependent beta values
        tau_factor_i = tau_x * (beta_x_0(i) + beta_x_L(i)) + ...
                       tau_y * (beta_y_0(i) + beta_y_L(i)) + ...
                       tau_z * (beta_z_0(i) + beta_z_L(i));
        tau_factor_i_flat = tau_factor_i(:);
        
        % Compute the denominator for all modes (vectorized)
        denom = Kn2_flat - k^2 - 1i * k * tau_factor_i_flat;
        
        % Sum the contributions of all modes to get the Green's function at frequency f
        G = sum((source_mode_flat .* receiver_mode_flat) ./ (V .* Lambda_flat .* denom));
        
        P(i) = G;
    end
toc
end
