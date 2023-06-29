%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      This is the matlab code for fractinal equation          %%
%%    d^alp/dt u(x,t)-d^2/dx^2 u(x,t)=f(x,t), a1<=x<=a2, t>0    %%
%%                             u(x,0)=g(x),   a1<=x<=a2         %%
%%         d/dx u(x,t)=d^(alp/2)/dt u(x,t),   x=a1              %%
%%         d/dx u(x,t)=-d^(alp/2)/dt u(x,t),   x=a2             %%
%% by the method in the paper SIAM2016 of Lv&Xu(3-alp).         %%
%% In the first time step, sub-stepping scheme is used in order %%
%% to get the global 3-alp accuracy.                            %%
%% Author:Hongyi Zhu                                            %%
%% The remark begin with the * means for different problem,     %%
%% you need to change the parameter.                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [L2,H1,ratio_L2,ratio_H1]=FDSTFDE_NZB_3()
  %clear;
  %clc;
  format long 
  %define some parameters which need to be changed for different problem
  a1=-2; a2=2; %* domain of space
  T=1; %* domain of time
  alp=0.9; %* order of fractional derivative

  exact_u=@(x,t)(exp(-5*x.^2).*t.^(3+alp));
  %exact_u=@(x,t)((x.^4).*((pi-x).^4).*(exp(-x)*t^(3+alp)+1));%* exact solution
  %f=@(x,t)(gamma(4+alp)*(x.^4).*((pi-x).^4).*exp(-x)*t^3/6-...
   %  (x.^2).*((pi-x).^2).*(t^(3+alp)*exp(-x).*...
   % ((x.^2).*(56-16*x+x.^2)-2*pi*x.*(28-12*x+x.^2)+pi^2*(12-8*x+x.^2))+...
   % 4*(3*pi^2-14*pi*x+14*x.^2)));%* right hand term
  f=@(x,t)(exp(-5*x.^2).*t.^3.*(gamma(4+alp)/6-10*(10*x.^2-1).*t.^alp));
  
  Nt=1000; %number of step
  NxN=10:2:36;%* degree of polynomial use for space
  for nj=1:length(NxN)
      Nx=NxN(nj);
      %compute the GLL nodes,weights and the first-order differentiation matrix in [a1,a2]
      [xx,wx]= legs(Nx+1);%Legendre-Gauss points and weights in [-1,1]
      D=legsdiff(Nx+1,xx);%the first-order differentiation matrix in [-1,1]
      xx=(a2-a1)/2*xx+(a2+a1)/2; %Legendre-Gauss points for [a1,a2]
      wx=(a2-a1)/2*wx;%weights in [a1,a2]
  for ni=1:1%*
      dt=T/Nt; %* length of time step
      Nt1=ceil(dt^(alp-1));%number of step for the sub-stepping scheme
      dt1=dt/Nt1;%length of time step for the sub-stepping scheme
      u=zeros(Nx+1,Nt+1); %save numerical solution
      u1=zeros(Nx+1,Nt1+1);%save numerical solution of the sub-stepping scheme
      %u(:,1)=xx.^4.*(pi-xx).^4;%* initial valu
      u(:,1)=0;
      u1(:,1)=u(:,1);
      %four useful constant
      a0=gamma(3-alp)*dt^(alp);
      a0h=gamma(3-alp/2)*dt^(alp/2);
      aw0=gamma(2-alp)*dt1^(alp);
      aw0h=gamma(2-alp/2)*dt1^(alp/2);

      %start the sub-steping scheme to get u^{1},i.e. u^{1,Nt1}
      %compute coefficient matrix
      t1=zeros(Nx+1,1);
      t1(1)=aw0/aw0h;t1(Nx+1)=t1(1);
      A=diag(wx)+4/(a2-a1)^2*aw0*D'*diag(wx)*D+diag(t1);%coefficient matrix for the sub-stepping scheme  
      F=rightside11(u1,f,dt1,xx,wx,aw0,aw0h,Nx);%*compute right hand side for the first iteration of the sub-stepping scheme
      u1(:,2)=pcg(A,F,1e-12,50);%*get u^{1,1}
  
      b0=zeros(Nt1,1); 
      [b0]=coeff_b0(alp,Nt1);
      b02h=zeros(Nt1,1); 
      [b02h]=coeff_b0(alp/2,Nt1);
  
      for k=2:Nt1
          F=rightside12(b0,b02h,u1,k,f,dt1,xx,wx,aw0,aw0h,Nx); %*compute right hand side of the sub-stepping scheme after the first time step
          u1(:,k+1)=pcg(A,F,1e-12,50);%*get u^{1,k-1}
      end
      u(:,2)=u1(:,Nt1+1);%obtain u^{1,Nt1}

      %start the main scheme
      [a,b,c]=coeff_abc(alp,Nt);
      [ah,bh,ch]=coeff_abc(alp/2,Nt);
      beta0=c(1)+2-alp/2;
      beta0h=ch(1)+2-alp/4;

      t2=zeros(Nx+1,1);
      t2(1)=1;t2(Nx+1)=1;
      A=diag(wx)+4/(a2-a1)^2*a0/beta0*D'*diag(wx)*D+a0/a0h*beta0h/beta0*diag(t2);%coefficient matrix
      for k=2:Nt
          F=rightside(a,b,c,ah,bh,ch,k,u,f,dt,alp,xx,wx,a0,a0h,beta0,Nx);%*
          u(:,k+1)=pcg(A,F,1e-12,50);%*
      end
      %compute error in L^2, H^{1} norm 
      numeri_u=@(x)(lagrange(xx,u(:,Nt+1),x));% numerical solution
      [L2e,semiH1e,H1e,Linfty]=err(numeri_u,exact_u,a1,a2,T);
      %Nt=2*Nt;
  end
  L2(nj)=L2e;H1(nj)=H1e;
  end
  ratio_L2=diff(-log2(L2));
  ratio_H1=diff(-log2(H1));
%   figure(1)%time
%   kk=1./2.^(3:7);%*
%   loglog(kk,L2,'-*','LineWidth',2)
%   hold on
%   loglog(kk,H1,'-*','LineWidth',2)
%   hold on
%   loglog(kk,kk.^(3-alp)/exp(1.2),'--','LineWidth',1.5)
%   xlabel('Time step','fontsize',13,'FontWeight','bold');
%   ylabel('Error','fontsize',13,'FontWeight','bold');
%   legend('L^2 error','H^1 error','Slope:2.7','fontsize',13,'FontWeight','bold')%*
%   title('\alpha=0.3','fontsize',14.3,'FontWeight','bold')%*
%   return
  figure(2)%space
  semilogy(NxN,L2,'-*','LineWidth',2)
  hold on
  semilogy(NxN,H1,'-+','LineWidth',2)
  xlabel('Polynomial degree N','fontsize',13,'FontWeight','bold');
  ylabel('Error','fontsize',13,'FontWeight','bold');
  legend('L^2 error','H^1 error','fontsize',13,'FontWeight','bold')
  title('\alpha=0.9','fontsize',14.3,'FontWeight','bold')%*
  return
  
  %*save the time errors
  str=['timeerror-' num2str(alp)  '-' num2str(Nx) '-' num2str(T) '.txt'];
  fid=fopen(str,'a');
  t3=[dt L2e semiH1e H1e dt^(3-alp)];
  fprintf(fid,'%7.5f  %20.17e  %20.17e  %20.17e  %20.17e\n',t3);
  fclose(fid);
  
  %*save the space errors
  %str=['spaceerror-' num2str(alp)  '-' num2str(dt) '-' num2str(T) '.txt'];
  %fid = fopen(str,'a');
  %t3=[Nx L2e semiH1e H1e Linfty];
  %fprintf(fid,'%4d  %20.17e  %20.17e  %20.17e  %20.17e  ',t3);
  %fclose(fid);

  %plot the solution
  %fplot(exact_u,[a1,a2])
  %hold on
  %fplot(numeri_u,[a1,a2])
  %legend('exact solution','numerical solution');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F=rightside11(u1,f,dt1,xx,wx,aw0,aw0h,Nx)% compute right hand side for the fist iteration of the sub-stepping scheme
  f1=f(xx,dt1);
  f2=zeros(Nx+1,1);
  f2(1)=aw0/aw0h*u1(1,1);
  f2(Nx+1)=aw0/aw0h*u1(Nx+1,1);
  F=diag(wx)*(aw0*f1+u1(:,1))+f2;
end
%##############################################################

function F=rightside12(b0,b02h,u1,k,f,dt1,xx,wx,aw0,aw0h,Nx) %compute right hand side for other step
  f1=f(xx,k*dt1);
  if k==2
     F=diag(wx)*((1-b0(1))*u1(:,2)+b0(1)*u1(:,1)+aw0*f1);
     f11=zeros(Nx+1,1);
     f11(1)=(1-b02h(1))*u1(1,2)+b02h(1)*u1(1,1);
     f11(Nx+1)=(1-b02h(1))*u1(Nx+1,2)+b02h(1)*u1(Nx+1,1);
     F=F+aw0/aw0h*f11;
  else
     f2=0;
     for i=1:k-2
         f2=f2+(b0(i)-b0(i+1))*u1(:,k-i);
     end
     F=diag(wx)*((1-b0(1))*u1(:,k)+f2+b0(k-1)*u1(:,1)+aw0*f1);
     f22=zeros(Nx+1,1);
     for i=1:k-2
         f22(1)=f22(1)+(b02h(i)-b02h(i+1))*u1(1,k-i);
         f22(Nx+1)=f22(Nx+1)+(b02h(i)-b02h(i+1))*u1(Nx+1,k-i);
     end
     f22(1)=f22(1)+(1-b02h(1))*u1(1,k)+b02h(k-1)*u1(1,1);
     f22(Nx+1)=f22(Nx+1)+(1-b02h(1))*u1(Nx+1,k)+b02h(k-1)*u1(Nx+1,1);
     F=F+aw0/aw0h*f22;
  end
end
  
%####################################################################
function F=rightside(a,b,c,ah,bh,ch,k,u,f,dt,alp,xx,wx,a0,a0h,beta0,Nx)
  f1=0;
  d=coeff_d(a,b,c,k,alp);%d is different from theoretical! here d dosen't divide beta0
  for i=1:k
      f1=f1+d(i)*u(:,i);
  end
  F=diag(wx)*(f1+a0*f(xx,k*dt))/beta0; 
  dh=coeff_d(ah,bh,ch,k,alp/2);
  f2=zeros(Nx+1,1);
  for i=1:k
      f2(1)=f2(1)+dh(i)*u(1,i);
      f2(Nx+1)=f2(Nx+1)+dh(i)*u(Nx+1,i);
  end
  F=F+a0/a0h/beta0*f2;
end

%#################################################################
function [b0]=coeff_b0(alp,Nt1)
  for i=1:Nt1
    b0(i)=(i+1)^(1-alp)-i^(1-alp);%compute the coefficients b_j
  end
end

%################################################################
function [a,b,c]=coeff_abc(alp,Nt)
  for i=1:Nt-1
      a(i)=-3/2*(2-alp)*(i+1)^(1-alp)+(2-alp)*i^(1-alp)/2+(i+1)^(2-alp)-i^(2-alp);
      b(i)=2*(2-alp)*(i+1)^(1-alp)-2*(i+1)^(2-alp)+2*i^(2-alp);
      c(i)=-a(i)-b(i);
  end
end
    
%################################################################
function d=coeff_d(a,b,c,k,alp)
  if k==2
     d(1)=-a(1)-alp/2;
     d(2)=-b(1)+2;
  elseif k==3
     d(1)=-a(2);
     d(2)=-a(1)-b(2)-alp/2;
     d(3)=-b(1)-c(2)+2;
  else
     d(1)=-a(k-1);
     d(2)=-a(k-2)-b(k-1);
     for j=3:k-2
         d(j)=-a(k-j)-b(k-j+1)-c(k-j+2);
     end
     d(k-1)=-a(1)-b(2)-c(3)-alp/2;
     d(k)=-b(1)-c(2)+2;
  end
end
   %###################################################################
function [L2e,semiH1e,H1e,Linfty]=err(f1,f2,a1,a2,T)%f1 is the numerical solution, f2 is the exact solution
  [xx,wx]= legslb(127);
  D=legslbdiff(length(xx),xx);%the first-order differentiation matrix in [-1,1]
  xx=(a2-a1)/2*xx+(a2+a1)/2; %Legendre-Gauss-Lobatto for [a1,a2]
  wx=(a2-a1)/2*wx;%weights in [a1,a2]
  
  tmp1=f1(xx)-f2(xx,T);
  L2e=(abs(tmp1))'*diag(wx)*(abs(tmp1)); %L^2 error 
  
  tmp=D*tmp1;
  semiH1e=4/(a2-a1)^2*tmp'*diag(wx)*tmp;
  
  H1e=L2e+semiH1e;
  L2e=sqrt(L2e);
  semiH1e=sqrt(semiH1e);
  H1e=sqrt(H1e);
  
  x1=a1:0.001:a2;
  Linfty=max(abs(f1(x1)-f2(x1,T)));  
end