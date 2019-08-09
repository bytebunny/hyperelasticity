clear variables
close all
format short
clc

addpath( genpath('../CA2/') ) % Add path to Magnus's routines.

% Initial and current positions from example 9.1:
xy0=[0,0; 4,0; 0,3];
xy=[2,3; 10,3; 10,9];

[N,dN_dX,dN_dx,F]=shape_gradients(xy0,xy);

dN_dX

dN_dx

% For the plane strain condition, the deformation gradient becomes:
F_temp = [F, [0; 0]];
F = [F_temp; [0 0 1]];

% Use Magnus's Neo-Hookean model to compute the 2nd Piola-Kirchhoff stress:
C = m_2_v9(F'*F);
[PK2,~] = neo_hooke(C,[3,2])