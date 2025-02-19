function [ac, Z_in] = PorousImpedance(ff, sigma, h, config, rho0, c0, d_air, method)
% PorousAbsorption computes the absorption coefficient and specific acoustic
% impedance of a porous layer using a transfer matrix method for either a 
% double-layer (porous layer backed by a rigid wall) or triple-layer (porous layer 
% with an air gap in front of a rigid wall) configuration.
% Reference: https://doi.org/10.1016/j.apacoust.2013.06.004
% Usage:
%   [ac, Z_in] = PorousAbsorption(ff, sigma, h, config, rho0, c0, d_air, method)
%
% Inputs:
%   ff     - Frequency vector [Hz]
%   sigma  - Flow resistivity [Pa·s/m^2]
%   h      - Thickness of the porous material [m]
%   config - Layer configuration: 2 for double-layer (default) or 
%            3 for triple-layer (porous + air gap + rigid wall)
%   rho0   - Air density [kg/m^3] (default: 1.204)
%   c0     - Speed of sound [m/s] (default: 343)
%   d_air  - Air gap thickness [m] (only required for config==3, default: 0.05)
%   method - Method for porous layer impedance calculation (default: 8)
%            1: Delany & Bazley (valid for 1e-2 <= X <= 1, where X=(rho0*f)/sigma)
%            2: Qunli (based on measurements of porous foams, typically 200-2000 Hz)
%            3: Miki (regression based on Delany & Bazley data; always positive absorption)
%            4: Mechel (regression method; note potential discontinuities at transitions)
%            5: Mechel (piecewise, using different coefficients for X<=0.025 and X>0.025,
%                       with coefficients chosen for glass fiber or rock fiber)
%            6: Komatsu (based on measurements from glass and rock wools; typically valid for
%                       sigma in the range 6000 to 72900 Pa·s/m^2)
%            7: Allard & Champoux (theoretical method requiring dynamic density and bulk modulus)
%            8: Modified Allard & Champoux (default; calibrated via regression)
%
% Outputs:
%   ac   - Absorption coefficient at each frequency
%   Z_in - Specific acoustic impedance at each frequency [Pa·s/m]
%
% Notes on method limitations:
%   - Method 1 (Delany & Bazley): Best for dimensionless frequency X=(rho0*f)/sigma in [1e-2, 1].
%   - Method 2 (Qunli): Developed from measurements on porous foams; typically valid from 200 to 2000 Hz.
%   - Method 3 (Miki): Provides always positive absorption coefficients; valid for a slightly
%                      larger frequency range than Method 1.
%   - Method 4 (Mechel): Uses regression with potential discontinuities at the low-high frequency transition.
%   - Method 5 (Mechel, piecewise): Uses separate regression coefficients for X<=0.025 and X>0.025.
%   - Method 6 (Komatsu): Based on specific airflow resistivity values from 6000 to 72900 Pa·s/m^2.
%   - Method 7 (Allard & Champoux): Theoretical model based on dynamic density and bulk modulus.
%   - Method 8 (Modified Allard & Champoux): Calibrated version of Method 7 using regression; default choice.
%
% Example:
%   f = linspace(100, 4000, 100)';
%   [ac, Z_in] = PorousAbsorption(f, 15000, 0.05, 2);
%

% Set default input values if not provided
if nargin < 4, config = 2; end
if nargin < 5, rho0 = 1.204; end
if nargin < 6, c0 = 343; end
if nargin < 7
    if config == 3
        d_air = 0.05;
    else
        d_air = [];
    end
end
if nargin < 8, method = 8; end

% Number of frequency points
N = length(ff);
ac   = zeros(N, 1);
Z_in = zeros(N, 1);

% Loop over frequencies
for i = 1:N
    f = ff(i);
    omega = 2 * pi * f;
    
    % Rigid backwall impedance using a simple rigid-wall absorption model
    alpha = RigidwallImpedance(f);
    R = sqrt(1 - alpha);
    Z1 = rho0 * c0 * (1 + R) / (1 - R);
    
    % Calculate porous layer impedance and propagation constant using the selected method
    switch method
        case {1,2,3,4,8}  % Methods with similar formulation
            X = (rho0 * f) / sigma;
            if method == 1
                % Delany & Bazley (Method 1)
                % Valid for: 1e-2 <= X <= 1
                c1 = 0.05710; c2 = -0.75400; c3 = 0.08700; c4 = -0.73200;
                c5 = 0.18900; c6 = -0.59500; c7 = 0.09780; c8 = -0.70000;
            elseif method == 2
                % Qunli (Method 2)
                % Typically valid for frequencies between 200 and 2000 Hz.
                c1 = 0.20900; c2 = -0.54800; c3 = 0.10500; c4 = -0.60700;
                c5 = 0.16300; c6 = -0.59200; c7 = 0.18800; c8 = -0.55400;
            elseif method == 3
                % Miki (Method 3)
                % Always positive absorption coefficients; slightly extended frequency range.
                c1 = 0.07000; c2 = -0.63200; c3 = 0.10700; c4 = -0.63200;
                c5 = 0.16000; c6 = -0.61800; c7 = 0.10900; c8 = -0.61800;
            elseif method == 4
                % Mechel (Method 4)
                % Regression method with possible discontinuities at transitions.
                c1 = 0.06080; c2 = -0.71730; c3 = 0.13230; c4 = -0.66010;
                c5 = 0.20820; c6 = -0.61930; c7 = 0.10870; c8 = -0.67310;
            elseif method == 8
                % Modified Allard & Champoux (Method 8)
                % Calibrated via regression; default method.
                c1 = 0.07290; c2 = -0.66228; c3 = 0.18700; c4 = -0.53790;
                c5 = 0.28800; c6 = -0.52600; c7 = 0.09820; c8 = -0.68500;
            end
            Z_absorber = rho0 * c0 * (1 + c1 * X^(c2) - 1i * c3 * X^(c4));
            gamma = (omega / c0) * ( c5 * X^(c6) + 1i * (1 + c7 * X^(c8)) );
            
        case 5  % Method 5: Mechel piecewise (using Glass fiber coefficients)
            X = (rho0 * f) / sigma;
            if X <= 0.025
                % Low-frequency coefficients for Glass fiber
                c1 = 0.06880; c2 = -0.70700; c3 = 0.19600; c4 = -0.54900;
                c5 = 0.39600; c6 = -0.45800; c7 = 0.13500; c8 = -0.64600;
            else
                % High-frequency coefficients for Glass fiber
                c1 = 0.02350; c2 = -0.88700; c3 = 0.08750; c4 = -0.77000;
                c5 = 0.17900; c6 = -0.67400; c7 = 0.10200; c8 = -0.70500;
            end
            Z_absorber = rho0 * c0 * (1 + c1 * X^(c2) - 1i * c3 * X^(c4));
            gamma = (omega / c0) * ( c5 * X^(c6) + 1i * (1 + c7 * X^(c8)) );
            
        case 6  % Komatsu's method
            % Based on measurements for glass and rock wools (sigma ~6000 to 72900 Pa·s/m^2)
            X = 2 - log(f/sigma);
            c1 = 0.00027; c2 = 6.20000; c3 = 0.00470; c4 = 4.10000;
            c5 = 0.00690; c6 = 4.10000; c7 = 0.00040; c8 = 6.20000;
            Z_absorber = rho0 * c0 * (1 + c1 * X^(c2) - 1i * c3 * X^(c4));
            gamma = (omega / c0) * ( c5 * X^(c6) + 1i * (1 + c7 * X^(c8)) );
            
        case 7  % Theoretical Allard & Champoux method
            % Requires calculation of dynamic density and bulk modulus.
            X = (rho0 * f) / sigma;
            P0 = 101320;  % Atmospheric pressure [Pa]
            rho_dyn = 1.2 + sqrt(-0.0364 * X^(-2) - 1i * 0.1144 * X^(-1));
            K_dyn = P0 * (1i * 29.64 + sqrt(2.82 * X^(-2) + 1i * 24.9 * X^(-1))) / ...
                    (1i * 21.17 + sqrt(2.82 * X^(-2) + 1i * 24.9 * X^(-1)));
            Z_absorber = sqrt(rho_dyn * K_dyn);
            gamma = 1i * 2 * pi * f * sqrt(rho_dyn / K_dyn);
            
        otherwise
            error('Unknown method specified.');
    end
    
    % Use transfer matrix method to compute overall impedance
    if config == 2
        % Double-layer configuration: porous layer backed by a rigid wall
        Z_in(i) = Z_absorber * (Z1 * cosh(gamma * h) + Z_absorber * sinh(gamma * h)) / ...
                  (Z1 * sinh(gamma * h) + Z_absorber * cosh(gamma * h));
    elseif config == 3
        % Triple-layer configuration: porous layer, air gap, then rigid wall.
        Z_air = rho0 * c0;
        gamma_air = 1i * 2 * pi * f / c0;
        Z_s2 = Z_air * (Z1 * cosh(gamma_air * d_air) + Z_air * sinh(gamma_air * d_air)) / ...
               (Z1 * sinh(gamma_air * d_air) + Z_air * cosh(gamma_air * d_air));
        Z_in(i) = Z_absorber * (Z_s2 * cosh(gamma * h) + Z_absorber * sinh(gamma * h)) / ...
                  (Z_s2 * sinh(gamma * h) + Z_absorber * cosh(gamma * h));
    else
        error('Unsupported configuration. Use 2 for double-layer or 3 for triple-layer.');
    end
    
    % Convert to specific surface impedance
    Z_in(i) = Z_in(i)./(rho0 * c0);
    % Compute absorption coefficient
    ac(i) = 1 - abs((Z_in(i) - 1) / (Z_in(i) + 1))^2;
end
function alpha = RigidwallImpedance(f)
    % Experimental validation and uncertainty quantification in wave-based computational room acoustics
    if f < 125*2^(-1/2) % below 125 octave band
        alpha = 0.01;
    elseif f >= 125*2^(-1/2) && f <= 125*2^(1/2)
       alpha =   0.01;
    elseif f >= 250*2^(-1/2) && f <= 250*2^(1/2)
       alpha =   0.01;
    elseif f >= 500*2^(-1/2) && f <= 500*2^(1/2)
       alpha =  0.01;
    else 
        alpha = 0.01;
    end
end


end
