clear all
close all
format short e
%%%

addpath( genpath('../hyperelasticity/') ) % Add path to my routines. /Rostyslav

%definition of time history
tmax=100; no_of_timesteps=1000;
t=linspace(0,tmax,no_of_timesteps);

%max value of engineering strain
eps_max=1.5;


%definition of material parameters
Emod=20e0; v=0.45;
mu=Emod/(2*(1+v)); lambda=Emod*v/((1-2*v)*(1+v));
para=[mu; lambda; -mu/10; mu/30 ];

%small strain stiffness tensor (linear isotropic Hooke elasticity)
ident=[1 1 1 0 0 0 0 0 0]';
C_hooke=2*mu* f9_open_u_9(ident,ident)+lambda*ident*ident';

for i=1:no_of_timesteps
   %deformation gradient  as a 3x3 matrix
   F11=1+eps_max*(i-1)/no_of_timesteps;
   Fm=... %[F11 0 0; 0 F11^(-v) 0; 0 0 F11^(-v)];
   [1+eps_max*(i-1)/no_of_timesteps 0 0; 0 1 0; 0 0 1];
   %translation to 9 vector representation
   F=m_2_v9(Fm);
   %right Cauchy Green deformation tensor as a 3x matrix
   Cm=Fm'*Fm;
   %translation to 9 vector representation
   C=m_2_v9(Cm);
   
   %call constitutive model
   %[S2,dS2_dE]=neo_hooke(C,para(1:2));
   [S2,dS2_dE]=yeoh(C,para);
   
   % compare dS2_dE and C_hooke 
   %dS2_dE 
   %C_hooke
   %pause
   
   %numerical computation of stifness tensor to check that the
   %implementation
   %for jj=1:9
   %    Cdiff=C; num_pert=1.e-7; Cdiff(jj)=Cdiff(jj)+num_pert;
   %    [S2diff]=neo_hooke(Cdiff,para);
   %    dS2_dEdiff(:,jj)=2*(S2diff-S2)/num_pert;
   %end
   %dS2_dEdiff
   %dS2_dE
   %pause
   
   %compute Cauchy stress (as a 3x3 matrix) as a push forward of S2
   sigma=1/det(Fm)*Fm*v9_2_m(S2)*Fm';
   
   %compute the sitffness c by push forward operation
   Ft=m_2_v9(Fm');
   c=1/det(Fm)*f9_open_u_9(F,F)*dS2_dE*f9_open_u_9(Ft,Ft);
   
   
   %save plot data
   sigma_out(i)=sigma(1,1); strain_out(i)=Fm(1,1)-1; 
   sigma_out2(i)=sigma(2,2);
end

figure(1)
hold on
plot(strain_out,sigma_out,'b-')
figure(2)
hold on
plot(strain_out,sigma_out2,'b-')






