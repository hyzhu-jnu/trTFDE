%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   This is the matlab code for fractinal equation            %%
%%    d^alp/dt u(x,t)-d^2/dx^2 u(x,t)=f(x,t), a1<=x<=a2, t>0   %%
%%                             u(x,0)=g(x),   a1<=x<=a2        %%
%%         d/dx u(x,t)=d^(alp/2)/dt u(x,t),   x=a1             %%
%%         d/dx u(x,t)=-d^(alp/2)/dt u(x,t),   x=a2            %%
%% by the method in the paper SIAM2016 of Lv&Xu(3-alp) with    %%
%% grided mesh in time discretization.                         %% 
%% In the first time step, we just use the L1 scheme.          %%
%% the title NS_FSDE_NZB_3_1D means: NS:nonsmooth solution,    %%
%% FSDE:fractional sub-diffusion equation,                     %%
%% NZB:non zero bound                                          %%
%% 3:3-alp scheme in time, 1D: 1 dimension in space.           %%
%% Author:Hongyi Zhu                                           %%
%% The remark begin with the * means for different problem,    %%
%% you need to change the parameter.                           %%
%% this code is also applicable for equal mesh(smooth solution)%%
%% only by taking gam=1.                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NS_FSDE_NZB_3_1D(Nt,bs,exact_u)%[L2e,semiH1e,H1e,Linfty][u]
  fclose('all');
  %clear;
  %clc;
  format long 
  global tau a0 a0h aw0 aw0h wx xx Nx
  %define some parameters which need to be changed for different problem
  a1=-1.5; a2=1; %* domain of space
  T=1; %* domain of time
  alp=0.8; %* order of fractional derivative
  Nx=24; %* degree of polynomial use for space

  %Nt=10000; %* length of time step
  gam=(3-alp)/alp;%grided mesh
  tau(1)=T*(1/Nt)^gam;%step size tau(1)
  j=1:Nt; tt=T*(j/Nt).^gam; tt=[0 tt];
  j=j.^gam;
  tau(2:Nt)=tau(1).*diff(j);%step size tau(2:Nt)
  %two useful constants
  a0=gamma(3-alp)*tau(1)^(alp);
  a0h=gamma(3-alp/2)*tau(1)^(alp/2);
  aw0=a0/(2-alp);
  aw0h=a0h/(2-alp/2);
 
  %*the exact solution and the right hand side term f
  %exact_u=@(x,t)((1+t^alp+t^(1+alp))*sin(x));
  %d_exact_u=@(x,t)((1+t^alp+t^(1+alp))*cos(x));
  f=@(x,t)(0*x);
  
  %compute the GLL nodes,weights and the first-order differentiation matrix in [a1,a2]
  [xx,wx]= legslb(Nx+1);%Legendre-Gauss-Lobatto points and weights in [-1,1]
  xx=(a2-a1)/2*xx+(a2+a1)/2; %Legendre-Gauss-Lobatto for [a1,a2]
  wx=(a2-a1)/2*wx;%weights in [a1,a2]
  D=FADM(1,xx,xx);%the first-order differentiation matrix 
  
  u=zeros(Nx+1,Nt+1);%the numerical solution
  %u(:,1)=exact_u(xx,0);%the initial value
  for i=1:length(xx)
      if xx(i)<-1
          u(i,1)=0;
      elseif xx(i)<0
          u(i,1)=2*(xx(i)+1);
      elseif xx(i)<0.5
          u(i,1)=-4*xx(i)+2;
      else
          u(i,1)=0;
      end
  end

  %the first time step
  %compute coefficient matrix
  t1=zeros(Nx+1,1);
  t1(1)=aw0/aw0h;t1(Nx+1)=t1(1);
  A=diag(wx)+aw0*D'*diag(wx)*D+diag(t1);%coefficient matrix for the first step  
  F=rightside1(u,f);
  u(:,2)=pcg(A,F,1e-11,50);%*get u^{:,1}
  
  %exact_u=dlmread('uexact.txt');%size of exact_u:(Nx+1,10001)
  %bs=(size(exact_u,2)-1)/Nt;
  tmp1=u(:,2)-exact_u(:,bs+1);
  L2e=(abs(tmp1))'*diag(wx)*(abs(tmp1));
  e=L2e;%maximum L2 error
    
  %start the main scheme,from the second step on
  time=tau(1);
  t2=zeros(Nx+1,1);
  t2(1)=1;t2(Nx+1)=1;
  for k=2:Nt  
      k
      time=time+tau(k);
      %[betak,xi]=coefficient_d(alp,gam,k);%direct calculation
      [betak,xi]=coeffi_gk(alp,k,tt,tau,gam);%adaptive Gauss-Kronrod quadrature method
      [betakh,xih]=coeffi_gk(alp/2,k,tt,tau,gam);
      A=diag(wx)+a0/betak*D'*diag(wx)*D+a0/a0h*betakh/betak*diag(t2);%coefficient matrix
      F=rightside(betak,betakh,xi,xih,k,u,f,time);
      u(:,k+1)=pcg(A,F,1e-13,50);%     
      
      tmp1=u(:,k+1)-exact_u(:,k*bs+1);
      L2e=(abs(tmp1))'*diag(wx)*(abs(tmp1));
      e=max([L2e,e]);%maximum L2 error
  end  
  dlmwrite('maxL2_error',[Nt,e],'-append','precision','%7d %19.15f','newline','pc');
  return
  %dlmwrite('uexact.txt',u,'precision','%.15f');
  %return
  
  %compute error in L^2, L^{inf}, semiH^{1}, H^{1} norm
%   exact_u=dlmread('uexact.txt');
%   tmp1=u(:,end)-exact_u;
%   L2e=(abs(tmp1))'*diag(wx)*(abs(tmp1)); %L^2 error 
%   tmp=D*tmp1;
%   semiH1e=4/(a2-a1)^2*tmp'*diag(wx)*tmp;
%   H1e=L2e+semiH1e;
%   L2e=sqrt(L2e);
%   semiH1e=sqrt(semiH1e);
%   H1e=sqrt(H1e);

  %plot the solution
  %fplot(exact_u,[a1,a2])
  %hold on
  %fplot(numeri_u,[a1,a2])
  %legend('exact solution','numerical solution');
  
  %*save the time errors
  str=['timeerror-' num2str(alp)  '-' num2str(Nx) '-' num2str(T) '.txt'];
  fid = fopen(str,'a');
  t3=[Nt L2e semiH1e H1e Linfty];
  fprintf(fid,'%8d  %19.17e  %19.17e  %19.17e  %19.17e\n',t3);
  fclose(fid);

  %*save the space errors
  %str=['spaceerror-' num2str(alp)  '-' num2str(dt) '-' num2str(T) '.txt'];
  %fid = fopen(str,'a');
  %t3=[Nx L2e semiH1e H1e Linfty];
  %fprintf(fid,'%4d  %20.17e  %20.17e  %20.17e  %20.17e\n',t3);
  %fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F=rightside1(u,f) % compute right hand side for the fist iteration
  global tau aw0 wx xx
  f1=f(xx,tau(1));
  F=diag(wx)*(aw0*f1+u(:,1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F=rightside(betak,betakh,xi,xih,k,u,f,time)
  global wx xx Nx a0 a0h 
  f1=0;
  for i=1:k
      f1=f1+xi(i)*u(:,i);
  end
  F=diag(wx)*(f1+a0/betak*f(xx,time));
  f2=zeros(Nx+1,1);
  for i=1:k
      f2(1)=f2(1)+xih(i)*u(1,i);
      f2(Nx+1)=f2(Nx+1)+xih(i)*u(Nx+1,i);
  end
  F=F+a0/a0h/betak*betakh*f2;
end

%###################################################################
function [L2e,semiH1e,H1e,Linfty]=err(exact_u,d_exact_u,u,ax,bx,T,xx)
  [e,we]= legslb(127);
  ex=(bx-ax)/2*e+(bx+ax)/2;  wex=(bx-ax)/2*we;
  %*value of exact solution at ex
  EU=exact_u(ex,T);
  %value of numerical solution at ex
  PU=FADM(0,xx,ex)*u;
  %L^2-error
  L2e=wex'*(EU-PU).^2;
  %H^1-error
  EUX=d_exact_u(ex,T);
  PUX=FADM(1,xx,ex)*u;
  semiH1e=wex'*(EUX-PUX).^2;
  H1e=L2e+semiH1e;
  L2e=sqrt(L2e);semiH1e=sqrt(semiH1e);H1e=sqrt(H1e);
  %L^inf-error
  x=ax:0.01:bx;
  Linfty=max(abs(exact_u(x,T)'-FADM(0,xx,x)*u));
end