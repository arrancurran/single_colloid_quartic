function [ filtered_signal , dwell_times ] = schmitt_filter( signal , frames , dt , schmitt_trigger )
% Low pass filter with filtered_signalmitt trigers
% designed to remove Brownian noise from a colloids position
%

  % inits
  dwell_count               = 1 ;
  tran_count_h2l            = 0 ;                                            
  tran_count_l2h            = 0 ;
  
  filtered_signal           = zeros( frames - 1 , 1 ) ;                     
  dwell_times               = zeros( frames - 1 , 2 ) ;                     % Measured dwell times sperated into dwell times before a l2h and h2l
  
  % Check at t = 0 and t = t
  
  if signal ( 1 ) > 0                                                        % upper state at start
    filtered_signal( 1 ) = 1 ;          
  else                                                                       % lower state at start
    filtered_signal( 1 ) = -1 ;
  end
  
  if signal ( frames - 1 ) > 0                                               % upper state at end
    filtered_signal( frames - 1 ) = 1 ;
  else                                                                       % lower state at end
    filtered_signal( frames - 1 ) = -1 ;
  end
  
  % Scan the rest of the signal
  for n = 2 : frames - 2
    
    if filtered_signal( n - 1 ) > 0
  
      % In the high state, so look for the next transition to the low state
      
      if signal( n ) > - schmitt_trigger && signal( n + 1 ) < - schmitt_trigger
    	  
        % The signal has transitioned below the negative trigger
        
        tran_count_h2l = tran_count_h2l + 1 ;                     % Count transition from high to low
  
        filtered_signal( n ) =  - 1 ;                             % NOTE: filtered_signal will transit one frame early...???
        
        dwell_times( tran_count_h2l , 2 ) =  dwell_count * dt ;   % Calculate dwell time for high state
        
        dwell_count = 1 ;                                         % Reset dwell time counter
        
      else
  
        % Still in the high state
  
        filtered_signal( n ) = 1 ;  
  
        dwell_count = dwell_count + 1	;
  
      end
    
    else
  
      % In the low state, so look for the next transition to the high state
      
      if signal( n ) < schmitt_trigger && signal( n + 1 ) > schmitt_trigger     
        
        % The signal has transitioned above the positive trigger
  
        tran_count_l2h = tran_count_l2h + 1 ;											% Count transition from low to high
        
        filtered_signal( n ) =  1 ;
        
        dwell_times( tran_count_l2h , 1 ) =  dwell_count * dt ;	  % Calculate dwell time for low state
  
        dwell_count = 1 ;																			    % Reset dwell time counter
    
      else
  
        % Still in the low state
      
        filtered_signal( n ) = - 1 ;
        
        dwell_count =  dwell_count + 1 ;
        
      end
      
    end
  
  end

end
