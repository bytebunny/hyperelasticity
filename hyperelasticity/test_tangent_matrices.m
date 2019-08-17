clear variables, close all, format short, clc
% Initial and current positions from example 9.1:
xy0=[0,0; 4,0; 0,3];
xy=[2,3; 10,3; 10,9];

[k_const, k_init, k_tot] = get_tangent_matrices(xy0, xy, [3,2], 1)