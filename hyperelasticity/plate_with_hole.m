%% PLATE_WITH_HOLE main script to solve task E of CA2.
% The naming convention is adopted from CALFEM manual version 3.4.
% /Rostyslav Skrypnyk

close all
clear variables
format compact
clc

addpath(genpath('~/Documents/MATLAB/calfem-3/')) % Add Calfem routines.
addpath( genpath('../CA2/') ) % Add path to Magnus's routines.

%% Pre-processing
save_to_file = false; % Save plot data.
% Load mesh data: Coord, Dof, dof_fixed, dof_free, dof_prescr, Edof, Ex, Ey
load('cass2_mesh_data.mat')
N_dof = max(max(Edof(:,2:end)));
N_el = length( Edof(:,1) ); % Number of elements.
N_el_nodes = 3; % Number of nodes in each element.
N_el_dof = 2; % Number of DOFs per node.
[el_x, el_y] = coordxtr(Edof, Coord, Dof, N_el_nodes);

% Load and time stepping:
thickness = 1; % [mm], thickness of the plate.
u_max = 30; % [mm], prescribed horizontal displacement.
dt = 0.01; % [s]. Many steps are needed when truss areas are different.
total_t = 1; % [s].
N_steps = round(total_t / dt);

% Velocity boundary conditions: velocity is used to conveniently set
% prescribed value by multiplying with time step.
bc = [ [dof_fixed; zeros(size(dof_fixed))]'; % DOF, prescribed value.
       [dof_prescr; ones(size(dof_prescr)) * u_max / total_t]' ];

% Material parameters:
mu = 6.9; % [MPa].
lambda = 62.1; % [MPa].
params = [mu, lambda, -0.1*mu, mu/30];

% Tolerances:
error_tol = 1.e-3; % Error tolerance for Newton's iteration.
max_iter_Newton = 20; % Maximum number of Newton iterations.
max_iter = 100;

%% Processing
% Initialise history and output variables:
u = zeros( numel(Dof),1 ); % Global displacement vector.
du = u; % Increment of displacement.
% History containers:
u_hist = zeros( 1, N_steps ); % All top nodes move together.
horizont_force_hist = zeros( 1, N_steps );

for step=1:N_steps % Time stepping.
    du = zeros( numel(Dof),1 ); % Do not use results from previous time step.
    du( bc(:,1) ) = bc(:,2)*dt; % Reason why BC are defined as velocities.
    u_el = extract(Edof,u); % Element displacements.
    
    %% Equilibrium iteration
    for iter=1:max_iter
        %% Newton's method
        if iter <= max_iter_Newton
            du_el = extract(Edof,du); % Incremental element displacements. 

            % Initialize sparse matrix triplets I,J,V 
            % for the global stiffness matrix:
            ind = 0; 
            I = zeros(N_el*(N_el_nodes * N_el_dof)^2,1); 
            J = zeros(N_el*(N_el_nodes * N_el_dof)^2,1); 
            V = zeros(N_el*(N_el_nodes * N_el_dof)^2,1);
         
            residual_vector = zeros(N_dof,1);
            
            %% Assemble matrices
            for i=1:N_el % Loop over elements:
                [force,K] = element_routine(el_x(i,:), el_y(i,:), ...
                                            u_el(i,:), du_el(i,:), ...
                                            params, thickness);
                % Assemble global stiffness matrix and RHS vector:
                for ii = 1:N_el_nodes * N_el_dof 
                    for jj = 1:N_el_nodes * N_el_dof
                        ind = ind+1; 
                        I(ind) = Edof(i,1+ii); % 1st column in Edof is element number.
                        J(ind) = Edof(i,1+jj); 
                        V(ind) = K(ii,jj);
                    end
                end
                residual_vector(Edof(i,2:end)) = ...
                    residual_vector(Edof(i,2:end)) + force;
            end
            K_system = sparse(I,J,V);
            
            delta_du = - K_system(dof_free,dof_free) \ ...
                        residual_vector(dof_free);
            % Update displacement increment:
            du(dof_free) = du(dof_free) + delta_du;
        else
            %% False position method
            if iter == max_iter_Newton+1
                fprintf(['/// Convergence was not reached within %d Newton iterations.',...
                         ' Switched to false position method.\n'],max_iter_Newton)
                du_a = du; du_b = du;
                du_a(dof_free) = -sign(delta_du) * min(abs(delta_du),50);
                du_b(dof_free) = sign(delta_du) * min(abs(delta_du),50);
            end
            du_el_a = extract(Edof,du_a);
            du_el_b = extract(Edof,du_b);
            residual_vector_a = zeros(N_dof,1);
            residual_vector_b = zeros(N_dof,1);
            for i=1:N_el % Loop over elements:
                [force_a,~] = element_routine(el_x(i,:), el_y(i,:), ...
                                               u_el(i,:), du_el_a(i,:), ...
                                               params, thickness);
                residual_vector_a(Edof(i,2:end)) = ...
                    residual_vector_a(Edof(i,2:end)) + force_a;
                
                [force_b,~] = element_routine(el_x(i,:), el_y(i,:), ...
                                              u_el(i,:), du_el_b(i,:), ...
                                              params, thickness);
                residual_vector_b(Edof(i,2:end)) = ...
                    residual_vector_b(Edof(i,2:end)) + force_b;
            end
            % If residuals have different sign:
            if sign(residual_vector_a(dof_free)) * sign(residual_vector_b(dof_free)) < 0
                % Find the root of a line equation:
                du(dof_free) = du_b(dof_free) - ...
                    residual_vector_b(dof_free) * ( du_b(dof_free)-du_a(dof_free) ) / ...
                    ( residual_vector_b(dof_free) - residual_vector_a(dof_free) );
                % Find residual for the new displacement increment:
                du_el = extract(Edof,du);
                residual_vector_c = zeros(N_dof,1);
                for i=1:N_el % Loop over elements:
                    [force_c,~] = element_routine(el_x(i,:), el_y(i,:), ...
                                                  u_el(i,:), du_el(i,:), ...
                                                  params, thickness);
                    residual_vector_c(Edof(i,2:end)) = ...
                        residual_vector_c(Edof(i,2:end)) + force_c;
                end
                residual_vector = residual_vector_c;
                % Update bounds to hunt down the equilibrium position:
                if sign(residual_vector_b(dof_free)) * sign(residual_vector_c(dof_free)) > 0
                    du_b(dof_free) = du(dof_free);
                else
                    du_a(dof_free) = du(dof_free);
                end
            else % Search wider:
                du_a(dof_free) = 1.25 * du_a(dof_free);
                du_b(dof_free) = 1.25 * du_b(dof_free);
            end
        end
        % Check if (internal - external) forces are zero:
        if norm(residual_vector(dof_free)) <= error_tol
            fprintf('/// Step %d converged in %d iterations.\n',...
                    step, iter)
            break
        end
    end % Equilibrium iteration.
    
    %% Update variables:
    u = u + du;
    % Save history:
    u_hist(step) = u(dof_prescr(1));
    horizont_force_hist(step) = sum(residual_vector(dof_prescr));

    if iter == max_iter
        error('/// Step %d did NOT converge in %d iterations.',...
              step, max_iter)
    end
end

%% Post-processing
figure(1) % Fig.3.8(b): vertical force vs deflection
plot([0, u_hist], [0, horizont_force_hist],'o-')  
xlabel('Horizontal displacement, [mm]')
ylabel('Reaction force, [N]')

grid on
if save_to_file % Create file for LaTeX.
    f_name = ['../doc/data/force_displacement.dat'];
    f_id = fopen(f_name,'w');
    header = '# Horiz displ, [mm]   Reaction force, [N]';
    
    fprintf(f_id, '%s\n', header);
    fprintf(f_id, '%.4f           %.4f\n', [[0, u_hist]; [0, horizont_force_hist]]);
    fclose(f_id);
    
    figure(2)
    % Plot undeformed shape:
    eldraw2(Ex, Ey,[2 1 0]) % Black.

    % Plot deformed shape:
    u_el_array = extract(Edof, u);
    eldisp2(Ex, Ey, u_el_array, [1,4,0], 1);        
            
    axis off
    set(gca,'color','none') % Make background transparent.            
    print('../doc/img/mesh', '-dpng', '-r600') % Colour PNG.
end
        