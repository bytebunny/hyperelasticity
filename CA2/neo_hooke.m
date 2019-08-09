function [S2,dS2_dE]=neo_hooke(C,para)

%input: 
%C=right Cauchy Green deformation tensor 9x1
%para = mtrl parameters vector [mu ; lambda] 2x1

%output:
%S2=2nd Piola-Kirchhoff stress tensor 9x1
%dS2_dE=stiffness 9x9

mu=para(1); lambda=para(2);

invC=m_2_v9( inv(v9_2_m(C) ) );
detC=det(v9_2_m(C)); J=sqrt(detC);

%2nd order identity tensor
ident=[1 1 1 0 0 0 0 0 0]';

%eqn (6.28)
S2=mu*(ident-invC)+lambda*log(J)*invC;

%eqn (6.30)
dS2_dE=lambda*invC*invC'+2*(mu-lambda*log(J))*f9_open_u_9(invC,invC);
%0.5*( f9_open_u_9(invC,invC)+f9_open_l_9(invC,invC) );


%sym dS2_dE, otherwise otherwise FE problem will not converge ... (sym
%assumption is assumed in the derivations of the element stiffness)
%
%rows
dS2_dE(4,:)=0.5*( dS2_dE(4,:)+dS2_dE(8,:) );
dS2_dE(8,:)=dS2_dE(4,:);
dS2_dE(5,:)=0.5*( dS2_dE(5,:)+dS2_dE(9,:) );
dS2_dE(9,:)=dS2_dE(5,:);
dS2_dE(6,:)=0.5*( dS2_dE(6,:)+dS2_dE(7,:) );
dS2_dE(7,:)=dS2_dE(6,:);
%columns
dS2_dE(:,4)=0.5*( dS2_dE(:,4)+dS2_dE(:,8) );
dS2_dE(:,8)=dS2_dE(:,4);
dS2_dE(:,5)=0.5*( dS2_dE(:,5)+dS2_dE(:,9) );
dS2_dE(:,9)=dS2_dE(:,5);
dS2_dE(:,6)=0.5*( dS2_dE(:,6)+dS2_dE(:,7) );
dS2_dE(:,7)=dS2_dE(:,6);
