% flockenspiel.m is a Matlab script written for the literature thesis 
% MODELING ACTIVE MATTER INTO MICROFLUIDIC LOGIC by M. Faase. 
% The code is adapted from the paper NOVEL TYPE OF PHASE TRANSITION IN A 
% SYSTEM OF SELF-DRIVEN PARTICLES by T. Vicsek et al.(1995).

% clear previous simulation settings
close all
clear all

pause(5)

%% Simulation settings
style    =   0;       % choose if particles are displayed as
                      % 1: arrows (quiver)
                      % 0: circles (plot)
path     =   0;       % show trajectory, leaves previous points on plot
                      % 1: yes (doesn't work well with arrows)
                      % 0: no
spd      =   0.05;    % changes simulation playback speed

ord      =   0;       % show order parameter in subplot
                      % 1: yes
                      % 0: no

 % density and noise parameters / scaling options
np       =   1.5;     % noise parameter (noise strength)
n        =   300;     % particle number (density)
r        =   1;       % radius of interaction range (alignment strength)
L        =   7;       % domain size
v        =   0.05;    % particle velocity

 % time steps
t_max    =   500;     % maximum timestep
step     =   1;       % step size of next iteration

 % periodicity
period   =   1;       % set boundary contitions
                      % 1: periodic boundary condition 
                      % 0: unlimited
                      
%% initialize model parameters
                      
 % generate initial state (x,y,th) of n particles
x        =   L * rand(1,n);          % n random numbers between 0 and L
                                     % particle x position
y        =   L * rand(1,n);          % n random numbers between 0 and L
                                     % particle y position
th       =   2*pi *(rand(1,n) - .5); % n random numbers between -pi and pi
                                     % particle angle theta

%% Simulation time

for time =   1:t_max  % iterate through time steps
    
    dist    =   pdist([x' y']);  % calculate pairwise euclidean distances 
                                 % between all particles; number of pairs  
                                 % is the number of combinations, C(n,2)
    
    % periodic boundary properties
    if period == 1
        % for x, identify position of neighbors beyond the boundary 
        % within range r
        per_x(x < r)               =   L + x(x < r);
        per_x(x > L-r)             =   x(x > L - r) - L;
        per_x(r <= x & x <= L - r) =   x(r <= x & x <= L - r);
        
        % for y, identify position of neighbors beyond the boundary 
        % within range r
        per_y(y < r)               =   L + y(y < r);
        per_y(y > L - r)           =   y(y > L - r) - L;
        per_y(r <= y & y <= L - r) =   y(r <= y & y <= L - r);
        
        period_D    =   pdist([per_x' per_y']); % calculate pairwise 
                                                % euclidean distances
                                                % between particles
                                                % over boundaries
        
        dist        =   min([dist; period_D]);  % update euclidean distance
                                                % double with found nearer 
                                                % distances 
    end
    
    sq               =   squareform(dist);      % Matrix representation 
                                                % for the distance double
                                                
    [list_1, list_2] =   find(0 < sq & sq < r); % Identify which pairs of 
                                                % particles are within 
                                                % range r of eachother
    
    for i = 1:n                                 % iterate through the 
                                                % number of particles
                                                
        list = list_1 (list_2 == i);            % call pair values that 
                                                % correspond to current 
                                                % iteration
                                                
        if ~isempty(list)
            % the angles of the found nearest neighbors are decomposed into 
            % unit vectors, which are then averaged and the average angle 
            % is recalculated
            th_a(i) = atan2(mean(sin(th(list))),mean(cos(th(list))));
                                                
        else
            % if there are no found found pairs within range r, the average 
            % angle is set to the current angle of particle i
            th_a(i) = th(i);               
        end
    end
    
    %% Position Update
    x   =    x + v*cos(th)*step;     % x and y positions are updated, 
    y   =    y + v*sin(th)*step;     % component-wise according to the 
                                     % previous angle, per time step 'step'
                                     % with velocity v
    if period == 1
        
        x(x < 0) = L + x(x < 0);     % find new neighbors beyond the 
        x(L < x) = x(L < x) - L;     % boundaries
        
        y(y < 0) = L + y(y < 0);
        y(L < y) = y(L < y) - L;
    end
    

    
    %% Quiverplot / plot
    % prepare scale-dependent unit vectors from theta
    v_x     =    cos(th);
    v_y     =    sin(th);
    % find average vector
    a_x     =    mean(v_x); 
    a_y     =    mean(v_y);
    
    % compute dot product of all vectors wrt the average
    % the dot product represents the differences in vector direction
    for i = 1:n
        diff(i) = dot([v_x(i) v_y(i)],[a_x a_y]);
    end
    
    % by averaging the absolute values of the dot product we have a value 
    % the degree of alignment of all values 
    chi = mean(abs(diff));
    
    if ord == 1 
        subplot(1,2,2);
        plot(time,chi,'.','MarkerSize',5,'Color','k')
        xlim([0 t_max]);
        ylim([0 1]);
        hold on
        subplot(1,2,1);
    end
        
    if style == 1
        quiver(x, y, v_x, v_y,'AutoScale','off')
        xlim([0 L]);
        ylim([0 L]);
        axis square
        pause(spd)
        if path == 1
            hold on
        end
    else
        plot(x , y,'.','MarkerSize',5)
        xlim([0 L]);
        ylim([0 L]);
        axis square
        pause(spd)
        if path == 1
            hold on
        end
    end
    
   

      
    %% Angle update    
    th  =    th_a + np*(rand(1,n) - 0.5); % the angles th are updated
                                          % according to the found averages
                                          % of the nearest neighbors, th_a.
                                          % There is an addition of noise 
                                          % that is limited by the scalar 
                                          % noise parameter np
end