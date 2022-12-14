clearvars ; close all ; clc ;
% constants ===============================================================
k_B = 1.3806503e-23 ; % Boltzman constant [ m ^ 2 kg / s ^ 2 / K ]
T = 300 ; % Temperature [ K ]

% inputs ==================================================================
a = 1.4e-6 ; % radius [ m ]
viscosity = 0.89e-3 ;   % [ Pa s = kg / m / s ]
mobility  = 1 / ( 6 * pi * viscosity * a ) ; % Mobility [ s / kg ]
rho = 1045 ; % Density of BangLabs silica spheres [ kg / m ^ 3 ]
mass = 4 * pi * a ^ 3 * rho / 3 ; % mass of sphere [ kg ]

% settings ================================================================
frames = 1 + 1e6 ; % Number of frames, plus one s.t. finally array is even  
dt = 1e-4; % time step, >> than the momentum relaxation time, m * D / k_B / T ~ 1e -7 [ s ]
repeats = 10 ; % Number of repated calculations * 'frames' for each seperation

% inits ===================================================================
positions_array = zeros ( frames , repeats ) ; % positions in 1 D
dwell_times_array = { 1 , repeats } ;
mod = zeros ( frames , 1 ) ; % positions in 1 D
x_ps = zeros( frames - 1 , 1 ) ;
vel = 0 ;

% params ==================================================================
k = 8e-7 ;
U = 6 * k_B * T ;
d = .5e-6 ;
trap_centre = d ; % Roots of quartic = trap centres [ m ]
schmitt_trigger = trap_centre * 0.1 ; % Schmitt threshold for filtering

% Start simulation ========================================================

sigma = sqrt( 2 * k_B * T * dt / mobility ) ; % [ m ^ 2 kg^ 2 / s ^ 2 ] this is equivalent to <n1 n2 > see correlatedNoise.m

mod_freq = .4 * pi ; 
mod_amp = .4*1e3 * sqrt( 2 * k_B * T * dt / mobility ) ;
tic
parfor q = 1 : repeats % loop simulation for 'repeats'

  particle_position = zeros ( frames , 1 ) ; % positions in 1 D
  noise  = sigma * randn( frames - 1 , 1 ) ; % correlated excitation
  
  clc ; disp(['Iteration ' num2str(q)]) % Display current iteration
  
  for n  = 1 : frames -1 % Calculate Langevin
    % modulation = mod_amp * sin( mod_freq * dt * n ) ; %sin( pi * particle(n) / 2 / d + mod_freq * dt * n  ) 
    % V = modulation + 
    
    % Arbitrary quartic
    % V = a_q * particle( n ) - b_q * particle( n ) ^ 3 ; % Quartic trap force
    
    % Physical quartic
    V =   - 4 * U * particle_position( n ) .* ( ( particle_position( n ) / d ) .^ 2 - 1 ) / d^2 ; %quartic
    
    % Double Gaussian 
    % V = - exp( - k * ( - d + particle_position( n ) )^2 / 2 / U ) * k * ( - d + particle_position( n ) ) - exp( - k * ( d + particle_position( n ) )^2 / 2 / U ) * k * ( d + particle_position( n ) ) ;
    
    dparticle = mobility * V * dt + mobility * noise( n , 1 ) ;
    
    % dvel = ( mobility * noise( n , 1 ) - dt * V - dt * vel / mobility ) / mass ;
		% vel = vel + dvel ;
		% dparticle = dt * vel ;
    particle_position( n + 1 ) = particle_position( n ) + dparticle ; % kick the particle position by dp
    % mod( n ) = modulation ;
    % fax( n ) = faxen  ;
  end
  
  positions_array( : , q ) = particle_position ;
  [ filtered_signal , dwell_times ] = schmitt_filter( particle_position , frames , dt , schmitt_trigger )

  dwell_times_array{ : , q } = dwell_times
  
end

toc
positions_array( 1 , : ) = [ ] ; % Trim starting zeros

% Calculate power spectrum=======================================================================================
% Calculate power spectrum for entire data set and smooth using blocking
x_fft = fft( positions_array ) ; % Fourier transform data
x_0 = fftshift( x_fft ) ; % Shift values to centre on zero freq
x_ps = x_0 .* conj( x_0 ) ;	% Power spectra

t = dt * ones( frames - 1 , 1) ; t( 1 ) = 0 ; t = cumsum( t ) ;
PS = mean( dt )^2 * mean( x_ps , 2 ) / max( t ) ;
fNyq    =   length( PS ) / max( t ) / 2 ;                         % Nyquist frequency

freq = ( [ 1 : length( PS ) ]  / max( t ) ) ;                      % 0-centered frequency range

ind  = find( freq <= fNyq ) ;                      % only to the Nyquist f
freq = rot90( freq( ind ) ) ;                        % rotate to keap everthing n x 1
PS = PS( ind , : ) ;

loglog( freq , PS , 'o' )
%
% A_guess = 100 ;                    % Initial guess for viscosity in mPa
% tau_guess = 0.5 ;                      % Initial guess for corner frequency
% par_guess = [ A_guess , tau_guess ] ;
% fit_opt = statset('MaxIter',1e4);
% dwell_prob_par = nlinfit( xout2 , hist2 , @exponential , par_guess , fit_opt ) ;
%
%
% dwell_prob_fit = exponential( dwell_prob_par , xout2 ) ;
%  loglog( xout2, exponential( [ 100 , 0.5 ] , xout2 ) , '-r', xout2,hist2,'.k')


% Plotting ==========================================================================
range = 1 : 1e5 ;
ptsize = 1 ;
linewdth = 1 ;

%Particle position-------------------------------------------------------------------------------------------------------------------------
figure

subplot ( 2 , 2 , 1 ) ;
plot1 = plot(t( range ) , positions_array( range , 1 ) , 'k' );%, t( range ) , Sch( range ) * trap_centre, '--b' , t( range ) , mod( range ) * trap_centre / max( mod ) , '-r' ) ;

subplot ( 2 , 2 , 2 ) ;
[ position_hist , position_hist_x ] = hist( positions_array , sqrt( frames - 1 ) ) ;
position_hist = mean( position_hist , 2 ) ;
potential = -log( position_hist / sum( position_hist ) ) ;

plot( position_hist_x , potential ,...
  '-' ,...
  'LineWidth', linewdth ,...
  'MarkerEdgeColor','k',...
  'MarkerFaceColor','g',...
  'MarkerSize', ptsize ) ;



