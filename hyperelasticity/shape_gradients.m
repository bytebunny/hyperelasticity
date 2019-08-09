function [shape_func, shape_grad0, shape_grad, ...
          deform_grad] = shape_gradients(xy0, xy)
%% SHAPE_GRADIENTS returns shape functions and different gradients according to task B of assignment 2.
% The shape functions are defined for Constant Strain Triangle (plane strain)
% element.
%
% Input:
% xy0 -- 2d matrix of nodal intial coordinates (each row is X-Y pair).
% xy -- 2d matrix of nodal current coordinates (each row is x-y pair).
%
% Output:
% shape_func -- 1d array of nodal shape functions for triangle with the
%               following node numbering:
%               3
%               1 2
% shape_grad0 -- 2d array of material derivatives of shape functions.
% shape_grad -- 2d array of spatial derivatives of shape functions.
% deform_grad -- deformation gradient.

% Since there is 1 Gauss point in CST, can plug in local coordinates
% directly:
shape_func = [1-1/3-1/3; 1/3; 1/3];

% Shape functions gradient w.r.t. local coordinates:
shape_grad_local = [-1 1 0;
                    -1 0 1];

shape_grad0 = inv(xy0'*shape_grad_local' )' * shape_grad_local;
% This many transpositions are kept in order to comply with the definitions
% of equations in the textbook.
shape_grad = inv(xy'*shape_grad_local')' * shape_grad_local;

deform_grad = xy' * shape_grad0';
      
end % function SHAPE_GRADIENTS.