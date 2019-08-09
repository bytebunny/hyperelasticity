function [PioKir2, dPioKir2_dE] = yeoh(C,params)
%% YEOH implements the Yeoh hyperelastic material model according to task A of assignment 2.
% Interface of the function is designed to comply with the provided test 
% program for material subroutine test_const.m.
%
% Input:
% C -- right Cauchy-Green deformation tensor (in 9-vector representation).
% params -- vector of material parameters: mu, lambda, c2, c3.
%
% Output:
% PioKir2 -- 2nd Piola-Kirchhoff stress tensor (1x9 vector).
% dPioKir2_dE -- material (Lagrangian) elasticity tensor (1x9 matrix).

mu = params(1);
lambda = params(2);
c2 = params(3);
c3 = params(4);

Cm = v9_2_m(C);
C_inv = m_2_v9( inv(Cm) );
ident=[1 1 1 0 0 0 0 0 0]';
trace_C = sum(C.*ident);
Jacobi = sqrt(det(Cm));

PioKir2 = (mu + 4*c2*(trace_C-3) + 6*c3*(trace_C-3)^2) * ident ...
          +(lambda*log(Jacobi) - mu) * C_inv;

dPioKir2_dE = 2*( (4*c2 + 12*c3*(trace_C-3)) * (ident*ident') ...
                 + (0.5*lambda - lambda*log(Jacobi) + mu) * (C_inv*C_inv') );
      
end % function YEOH.