function [k_const, k_initial, k] = get_tangent_matrices(xy0, xy, params, t)
%% GET_TANGENT_MATRICES returns constitutive and initial stress matrices for CST element according to task D of assignment 2.
%
% Input:
% xy0 -- 2d matrix of nodal intial coordinates (each row is X-Y pair).
% xy -- 2d matrix of nodal current coordinates (each row is x-y pair).
% params -- vector of material parameters: mu, lambda.
% t -- element thickness.
%
% Output:
% k_const -- constitutive tangent matrix [6 x 6].
% k_initial -- initial stress tangent matrix [6 x 6].
% k -- total tangent matrix [6 x 6].

w = 0.5; % Weight in Gauss quadrature with 1 integration point;

[~,~,dN_dx,F]=shape_gradients(xy0,xy);

% Number of dimensions and number of nodes:
[n_dim, n_nodes] = size(dN_dx);
n_dof = n_dim * n_nodes; % Number of degrees of freedom.

B = zeros(3,n_dim*n_nodes);
for node=1:n_nodes
    B(:,node*n_dim-1:node*n_dim) = [dN_dx(1,node)          0;
                                    0          dN_dx(2,node);
                                    dN_dx(2,node) dN_dx(1,node)];
end
   
shape_grad_local = [-1 1 0;
                    -1 0 1];
dx_dksi = xy'*shape_grad_local';

% For the plane strain condition, the deformation gradient becomes:
F_temp = [F, [0; 0]];
F = [F_temp; [0 0 1]];
[S2,dS2_dE]=yeoh(m_2_v9(F'*F),params);
% ///Debugging: compute numerical material stiffness for comparison with
% analytical dS2_dE:
% for jj=1:9
%       Cdiff=m_2_v9(F'*F); num_pert=1.e-7; Cdiff(jj)=Cdiff(jj)+num_pert;
%       [S2diff,~]=yeoh(Cdiff,params);
%       dS2_dEdiff(:,jj)=2*(S2diff-S2)/num_pert;
% end

%compute Cauchy stress (as a 3x3 matrix) as a push forward of S2
cauchy_stress=1/det(F)*F*v9_2_m(S2)*F';

% Compute spatial 4th order elasticity tensor as Piola push forward of 
% material elasticity tensor (see eq-n 6.14 in Bonet & Wood, 2nd edition):
Fv = m_2_v9(F);
Ft=m_2_v9(F');
D = 1/det(F)*f9_open_u_9(Fv,Fv)*dS2_dE*f9_open_u_9(Ft,Ft);
D = D([1,2,4],[1,2,4]); % Extract 2D part.

k_const = zeros(n_dof,n_dof);
k_initial = zeros(n_dof,n_dof);
coef = w * t * det(dx_dksi);
for node_a = 1:n_nodes
    for node_b = 1:n_nodes
        k_const(node_a*n_dim-1 : node_a*n_dim, node_b*n_dim-1:node_b*n_dim) = ...
            coef * B(:,node_a*n_dim-1:node_a*n_dim)' * D * B(:,node_b*n_dim-1:node_b*n_dim);
        k_initial(node_a*n_dim-1 : node_a*n_dim, node_b*n_dim-1:node_b*n_dim) = ...
            coef * [dN_dx(1,node_a) dN_dx(2,node_a)] * cauchy_stress(1:n_dim,1:n_dim) * ...
            [dN_dx(1,node_b); dN_dx(2,node_b)] * eye(2);
    end
end

k = k_const + k_initial;

end % function GET_TANGENT_MATRICES.