function [forces] = get_intern_eq_forces(xy0, xy, params, t)
%% GET_INTERN_EQ_FORCES returns internal equivalent nodal forces for CST element according to task C of assignment 2.
%
% Input:
% xy0 -- 2d matrix of nodal intial coordinates (each row is X-Y pair).
% xy -- 2d matrix of nodal current coordinates (each row is x-y pair).
% params -- vector of material parameters: mu, lambda.
% t -- element thickness.
%
% Output:
% forces -- internal equivalent nodal forces (1x9 vector).

w = 0.5; % Weight in Gauss quadrature with 1 integration point;

[~,~,dN_dx,F]=shape_gradients(xy0,xy);
shape_grad_local = [-1 1 0;
                    -1 0 1];
dx_dksi = xy'*shape_grad_local';

% For the plane strain condition, the deformation gradient becomes:
F_temp = [F, [0; 0]];
F = [F_temp; [0 0 1]];
[S2,~]=yeoh(m_2_v9(F'*F),params);
%compute Cauchy stress (as a 3x3 matrix) as a push forward of S2
cauchy_stress=1/det(F)*F*v9_2_m(S2)*F';

[n_nodes, n_dim] = size(xy);

forces = zeros(n_nodes*n_dim,1);
for i=1:n_nodes
forces(1+(i-1)*n_dim:i*n_dim) = cauchy_stress(1:2,1)*dN_dx(1,i) + ...
                              cauchy_stress(1:2,2)*dN_dx(2,i);
end
forces = forces * w * t * det(dx_dksi);
      
end % function GET_INTERN_EQ_FORCES.