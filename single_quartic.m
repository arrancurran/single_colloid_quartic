clearvars ; close all ; clc ;
% constants ===============================================================
k_B = 1.3806503e-23 ; % Boltzman constant [ m ^ 2 kg / s ^ 2 / K ]
T = 300 ; % Temperature [ K ]

% inputs ==================================================================
a = 1.2e-6 ; % radius [ m ]
viscosity = 0.89e-3 ;   % [ Pa s = kg / m / s ]
mobility  = 1 / ( 6 * pi * viscosity * a ) ; % Mobility [ s / kg ]
rho = 1045 ; % Density of BangLabs silica spheres [ kg / m ^ 3 ]
mass = 4 * pi * a ^ 3 * rho / 3 ; % mass of sphere [ kg ]
% settings ================================================================
frames = 1 + 1e6 ; % Number of frames, plus one s.t. finally array is even  
dt = 1e-3; % time step, >> than the momentum relaxation time, m * D / k_B / T ~ 1e -7 [ s ]
repeats = 10 ; % Number of repated calculations * 'frames' for each seperation
% inits ===================================================================
PARTICLE = zeros ( frames , repeats ) ; % positions in 1 D
mod = zeros ( frames , 1 ) ; % positions in 1 D
x_ps = zeros( frames - 1 , 1 ) ;
vel = 0 ;
hopcount_U = 0 ; % Hop count, ( Particle 1 Lower to Upper trap ) [ P1L_U P1U_L P2L_U P2U_L ]
hopcount_L = 0 ;

% params ==================================================================
k = 5e-7 ;
U = 8 * k_B * T ;
d = .4e-6 ;
trap_centre = d ; % Roots of quartic = trap centres [ m ]
schmitt = trap_centre * 0.1 ; % Schmitt threshold for filtering

% Start simulation ========================================================

sigma = sqrt( 2 * k_B * T * dt / mobility ) ; % [ m ^ 2 kg^ 2 / s ^ 2 ] this is equivalent to <n1 n2 > see correlatedNoise.m

mod_freq = .4 * pi ; 
mod_amp = .4*1e3 * sqrt( 2 * k_B * T * dt / mobility ) ;
tic
parfor q = 1 : repeats % loop simulation for 'repeats'

  particle = zeros ( frames , 1 ) ; % positions in 1 D
  noise  = sigma * randn( frames - 1 , 1 ) ; % correlated excitation
  
  clc ; disp(['Iteration ' num2str(q)]) % Display current iteration
  
  for n  = 1 : frames -1 % Calculate Langevin
    % modulation = mod_amp * sin( mod_freq * dt * n ) ; %sin( pi * particle(n) / 2 / d + mod_freq * dt * n  ) 
    % V = modulation + 
    
    % Arbitrary quartic
    % V = a_q * particle( n ) - b_q * particle( n ) ^ 3 ; % Quartic trap force
    
    % Physical quartic
    % V =   - 4 * U * particle( n ) .* ( ( particle( n ) / d ) .^ 2 - 1 ) / d^2 ; %quartic
    
    % Double Gaussian 
    V = - exp( - k * ( - d + particle( n ) )^2 / 2 / U ) * k * ( - d + particle( n ) ) - exp( - k * ( d + particle( n ) )^2 / 2 / U ) * k * ( d + particle( n ) ) ;
    
    dparticle = mobility * V * dt + mobility * noise( n , 1 ) ;
    
    % dvel = ( mobility * noise( n , 1 ) - dt * V - dt * vel / mobility ) / mass ;
		% vel = vel + dvel ;
		% dparticle = dt * vel ;
    particle( n + 1 ) = particle( n ) + dparticle ; % kick the particle position by dp
    % mod( n ) = modulation ;
    % fax( n ) = faxen  ;
  end
  
  PARTICLE( : , q ) = particle ;

  %%Schmitt filter ========================================================

  hoptime = 1 ; % dwell time for current event [ t1U t1L t2U t2L ]
  Sch    = zeros( frames - 1 , 2 ) ; % Schmitt trace for x and y
  dwell_t = zeros( frames - 1 , 1 ) ; % Dwell times
  
  % Check where the bead is at t = 0 and t = t
  if particle ( 1 ) > 0 % upper trap or lower trap at the start of the data
      Sch( 1 ) = 1 ;
  else
      Sch( 1 ) = -1 ;
  end

  if particle ( frames - 1 ) > 0 % upper trap or lower trap at the end of the data
      Sch( frames - 1 ) = 1 ;
  else
      Sch( frames - 1 ) = -1 ;
  end

  % Scan the rest of the data
  for n = 2 : frames - 2
    % In the upper trap, so look for the next hop to the lower trap
    if Sch( n - 1 ) > 0
      
      if particle( n ) > - schmitt && particle( n + 1 ) < - schmitt	% The bead has hoped from high to low
              
        Sch( n ) =  - 1 ; % Sch will transit one frame early...???
        hopcount_U = hopcount_U + 1 ; % Count hop from up to low
        dwell_time_U( hopcount_U ) =  hoptime * dt ; % Calculate dwell time for upper
        hoptime = 1 ; % Reset dwell time counter
        
        % Scan for phase
        
        if mod( n ) > 0 % phase is < pi

          if mod( n ) > mod( n + 1 ) % this means phase is < pi and > pi/2
            phase_U( hopcount_U ) = pi - asin( mod( n ) / max( mod ) ) ;
          else % phase is < pi/2 > 0
            phase_U( hopcount_U ) = asin( mod( n ) / max( mod ) ) ;
          end
            
        else % phase is > pi
          
          if mod( n ) > mod( n + 1 ) ; % this means phase is  < 3pi/2 and > pi
            phase_U( hopcount_U ) = pi + asin( - mod( n ) / max( mod ) ) ;
          else % phase is < 2pi > 3pi/2
            phase_U( hopcount_U ) = 2 * pi  - asin( - mod( n ) / max( mod ) ) ;
          end
        
        end
      
      else
        Sch( n ) = 1 ;  % Still in the upper trap
        hoptime = hoptime + 1	; % Frame counter in upper trap
      end	% Stop and move to next frame
          
    % In the lower trap, so look for the next hop to the upper trap
    else
      
      if particle( n ) < schmitt && particle( n + 1 ) > schmitt          % The bead has hoped from low to upper
        
        Sch( n ) =  1 ;
        hopcount_L = hopcount_L + 1 ;											% Count hop from low to upper
        dwell_time_L( hopcount_L ) =  hoptime * dt ;					% Calculate dwell time for lower
        hoptime = 1 ;																			% Reset dwell time counter
             
        if mod( n ) > 0 % this means phase is < pi
                 
          if mod( n ) > mod( n + 1 ) % this means phase is < pi and > pi/2
            phase_L( hopcount_L ) = pi - asin( mod( n ) / max( mod ) ) ;
           
          else % phase is < pi/2 > 0
            phase_L( hopcount_L ) = asin( mod( n ) / max( mod ) ) ;
          
          end
                  
        
        else % phase is > pi
                  
          if mod( n ) > mod( n + 1 ) % this means phase is  < 3pi/2 and > pi
            
            phase_L( hopcount_L ) = pi + asin( - mod( n ) / max( mod ) ) ;
            
          else % phase is < 2pi > 3pi/2
          
            phase_L( hopcount_L ) = 2 * pi  - asin( - mod( n ) / max( mod ) ) ;
            
          end
          
        end
        
      else
      
        Sch( n ) = - 1 ;                                                                          % Still in the lower trap
        
        hoptime =  hoptime + 1 ;
        
      end
      
    end
  
  end
  
end

PARTICLE( 1 , : ) = [ ] ; % Trim starting zeros

% Calculate power spectrum=======================================================================================
% Calculate power spectrum for entire data set and smooth using blocking
x_fft = fft( PARTICLE ) ; % Fourier transform data
x_0 = fftshift( x_fft ) ; % Shift values to centre on zero freq
x_ps = x_0 .* conj( x_0 ) ;	% Power spectra

t = dt * ones( frames - 1 , 1) ; t( 1 ) = 0 ; t = cumsum( t ) ;
PS = mean( dt )^2 * mean( x_ps , 2 ) / max( t ) ;
fNyq    =   length( PS ) / max( t ) / 2 ;                         % Nyquist frequency

freq = ( [ 1 : length( PS ) ]  / max( t ) ) ;                      % 0-centered frequency range

ind  = find( freq <= fNyq ) ;                      % only to the Nyquist f
freq = rot90( freq( ind ) ) ;                        % rotate to keap everthing n x 1
PS = PS( ind , : ) ;
toc

% Blocking data to smooth out noise
%{
block=100;                   % Blocking length
count=0;
x_ps_smooth = zeros(floor(length(PS)/block),1);
freq_smooth = zeros(floor(length(PS)/block),1);
for b=1:block:length(PS)-block;
  count=count+1;
  x_block=sum(PS(b:b+block))/block;
  freq_block=sum(freq(b:b+block))/block;
  x_ps_smooth(count)=x_block;
  freq_smooth(count)=freq_block;
end
%}
subplot(1,2,1)
loglog( freq , smooth( PS , 400 ,'loess' ) , 'o' )
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

subplot(1,2,2)
[hist1 xout1 ] = hist( PARTICLE , sqrt( frames - 1 ) ) ;
potential = -log( hist1 / sum( hist1 ) ) ;
plot( xout1, potential, '--rs','LineWidth', linewdth ,...
  'MarkerEdgeColor','k',...
  'MarkerFaceColor','g',...
  'MarkerSize', ptsize ) ;

%Particle position-------------------------------------------------------------------------------------------------------------------------
figure

[hist3 xout3 ] = hist( phase_U , 50 ) ;
[hist4 xout4 ] = hist( phase_L , 50 ) ;

subplot ( 2 , 2 , 1 ) ;
plot1 = plot(t( range ) , particle( range ) , 'k' );%, t( range ) , Sch( range ) * trap_centre, '--b' , t( range ) , mod( range ) * trap_centre / max( mod ) , '-r' ) ;

subplot ( 2 , 2 , 2 ) ;
[hist1 xout1 ] = hist( particle , sqrt( frames - 1 ) ) ;
potential = -log( hist1 / sum( hist1 ) ) ;
plot( xout1, potential, '--rs','LineWidth', linewdth ,...
  'MarkerEdgeColor','k',...
  'MarkerFaceColor','g',...
  'MarkerSize', ptsize ) ;

subplot ( 2 , 2 , 3 ) ;
[hist2 , xout2 ] = hist( dwell_time_L  , sqrt( length( dwell_time_L ) ) ) ;

plot( xout2 , hist2 , 's','LineWidth', linewdth ,...
  'MarkerEdgeColor','k',...
  'MarkerFaceColor','k',...
  'MarkerSize', ptsize ) ;

subplot( 2 , 2 , 4 )
plot(xout3 / pi, hist3 / sum( hist3 ) , 'o' , xout4 / pi, hist4 / sum( hist4 ) , 's' )
loglog( freq , PS )

