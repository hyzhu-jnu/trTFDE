%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is the fast algorithm matlab code for fractinal equation %%
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


function [L2,H1,ratio_L2,ratio_H1]=FDSTFDE_NZB_3_fast()
  fclose('all');
  %clear;
  clc;
  format long 
  %define some parameters which need to be changed for different problem
  a1=-2; a2=2; %* domain of space
  T=1; %* domain of time
  alp=0.1; %* order of fractional derivative
  
  bet=alp+1;beth=alp/2+1;reps=10^(-15);%*desired relative error in fast algorithm
  [nzt,nwt,Nact]= SOEappr(bet,reps,0.001,T);
  [nzth,nwth,Nacth]= SOEappr(beth,reps,0.001,T);
  exact_u=@(x,t)(exp(-5*x.^2).*t.^(3+alp));
  %exact_u=@(x,t)((x.^4).*((pi-x).^4).*(exp(-x)*t^(3+alp)+1));%* exact solution
  %f=@(x,t)(gamma(4+alp)*(x.^4).*((pi-x).^4).*exp(-x)*t^3/6-...
   %  (x.^2).*((pi-x).^2).*(t^(3+alp)*exp(-x).*...
   % ((x.^2).*(56-16*x+x.^2)-2*pi*x.*(28-12*x+x.^2)+pi^2*(12-8*x+x.^2))+...
   % 4*(3*pi^2-14*pi*x+14*x.^2)));%* right hand term
  f=@(x,t)(exp(-5*x.^2).*t.^3.*(gamma(4+alp)/6-10*(10*x.^2-1).*t.^alp));
  
  Nt=8; %number of step
  NxN=48;%* degree of polynomial use for space
  for nj=1:length(NxN)
      Nx=NxN(nj);
      %compute the GLL nodes,weights and the first-order differentiation matrix in [a1,a2]
      [xx,wx]=legs(Nx+1);%Legendre-Gauss points and weights in [-1,1]
      D=legsdiff(Nx+1,xx);%the first-order differentiation matrix in [-1,1]
      xx=(a2-a1)/2*xx+(a2+a1)/2; %Legendre-Gauss points for [a1,a2]
      wx=(a2-a1)/2*wx;%weights in [a1,a2]
  for ni=1:5
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
      F=rightside11(u1,f,Nx,dt1,aw0,aw0h,xx,wx);%*compute right hand side for the first iteration of the sub-stepping scheme
      u1(:,2)=pcg(A,F,1e-13,50);%*get u^{1,1}
  
      b0=zeros(Nt1,1); 
      [b0]=coeff_b0(alp,Nt1);
      b02h=zeros(Nt1,1); 
      [b02h]=coeff_b0(alp/2,Nt1);
  
      for k=2:Nt1
          F=rightside12(b0,b02h,u1,k,f,Nx,dt1,aw0,aw0h,xx,wx); %*compute right hand side of the sub-stepping scheme after the first time step
          u1(:,k+1)=pcg(A,F,1e-13,50);%*get u^{1,k-1}
      end
      u(:,2)=u1(:,Nt1+1);%obtain u^{1,Nt1}

      %start the main scheme
      [a,b,c,ta,tb,tc]=coeff_abc(alp,nzt,nwt,dt);
      [ah,bh,ch,tah,tbh,tch]=coeff_abc(alp/2,nzth,nwth,dt);

      t2=zeros(Nx+1,1);
      t2(1)=1;t2(Nx+1)=1;
      A=((4-alp)/2-tc)*diag(wx)+4/(a2-a1)^2*a0*D'*diag(wx)*D+...
        a0/a0h*((4-alp/2)/2-tch)*diag(t2);%coefficient matrix
      Uhist=zeros(Nx+1,Nact);
      Uhisth=zeros(Nx+1,Nacth);
      F=rightside21(ta,tb,tah,tbh,u,f,dt,Nx,alp,a0,a0h,,xx,wx);
      u(:,3)=pcg(A,F,1e-13,50);%*
      for i=1:Nact
          Uhist(:,i)=a(i)*u(:,1)-b(i)*u(:,2)+c(i)*u(:,3);
      end
      Uhisth(1,:)=ah(:)*u(1,1)-bh(:)*u(1,2)+ch(:)*u(1,3);
      Uhisth(Nx+1,:)=ah(:)*u(Nx+1,1)-bh(:)*u(Nx+1,2)+ch(:)*u(Nx+1,3);
      for k=3:Nt
          F=rightside22(ta,tb,tah,tbh,k,u,Uhist,Uhisth,Nact,Nacth,f,dt,Nx,alp,nzt,nwt,nzth,nwth,a0,a0h,,xx,wx);%*
          u(:,k+1)=pcg(A,F,1e-13,50);%*
      for i=1:Nact
          Uhist(:,i)=exp(-dt*nzt(i))*Uhist(:,i)+(a(i)*u(:,k-1)-b(i)*u(:,k)+c(i)*u(:,k+1));
      end
      Uhisth(1,:)=exp(-dt*nzth(:))'.*Uhisth(1,:)+(ah(:)*u(1,k-1)-...
                  bh(:)*u(1,k)+ch(:)*u(1,k+1))';
      Uhisth(Nx+1,:)=exp(-dt*nzth(:))'.*Uhisth(Nx+1,:)+(ah(:)*u(Nx+1,k-1)-...
                     bh(:)*u(Nx+1,k)+ch(:)*u(Nx+1,k+1))';
      end
      
      %compute error in L^2, L^{inf}, H^{1} norm
      numeri_u=@(x)(lagrange(xx,u(:,Nt+1),x));% numerical solution
      [L2e,semiH1e,H1e,Linfty]=err(numeri_u,exact_u,a1,a2,T);
      L2(ni)=L2e;H1(ni)=H1e;
      Nt=2*Nt;
  end
  end
  ratio_L2=diff(-log2(L2));
  ratio_H1=diff(-log2(H1));
  figure(1)%time
  kk=1./2.^(3:7);%*
  loglog(kk,L2,'-*','LineWidth',2)
  hold on
  loglog(kk,H1,'-*','LineWidth',2)
  hold on
  loglog(kk,kk.^(3-alp)/exp(1),'--','LineWidth',1.5)
  xlabel('Time step','fontsize',13,'FontWeight','bold');
  ylabel('Error','fontsize',13,'FontWeight','bold');
  legend('L^2 error','H^1 error','Slope:2.9','fontsize',13,'FontWeight','bold','Location','southeast')%*
  title('\alpha=0.1','fontsize',14.3,'FontWeight','bold')%*
  return
  figure(2)%space
  semilogy(NxN,L2,'-*','LineWidth',2)
  hold on
  semilogy(NxN,H1,'-+','LineWidth',2)
  legend('L^2 error','H^1 error','fontsize',13,'FontWeight','bold')
  return
  save('FDSTFDE_NZB_3_fast.mat');
  %plot the solution
  %fplot(exact_u,[a1,a2])
  %hold on
  %fplot(numeri_u,[a1,a2])
  %legend('exact solution','numerical solution');

  %*save the time errors
  str=['timeerror-' num2str(alp)  '-' num2str(Nx) '-' num2str(T) '.txt'];
  fid = fopen(str,'a');
  t3=[dt L2e semiH1e H1e Linfty dt^(3-alp)];
  fprintf(fid,'%7.5f  %19.17e  %19.17e  %19.17e  %19.17e  %19.17e\n',t3);
  fclose(fid);

  %*save the space errors
  %str=['spaceerror-' num2str(alp)  '-' num2str(dt) '-' num2str(T) '.txt'];
  %fid = fopen(str,'a');
  %t3=[Nx L2e semiH1e H1e Linfty];
  %fprintf(fid,'%4d  %20.17e  %20.17e  %20.17e  %20.17e',t3);
  %fclose(fid);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F=rightside11(u1,f,Nx,dt1,aw0,aw0h,xx,wx) % compute right hand side for the fist iteration of the sub-stepping scheme
  f1=f(xx,dt1);
  f2=zeros(Nx+1,1);
  f2(1)=aw0/aw0h*u1(1,1);
  f2(Nx+1)=aw0/aw0h*u1(Nx+1,1);
  F=diag(wx)*(aw0*f1+u1(:,1))+f2;
  
%##############################################################
function F=rightside12(b0,b02h,u1,k,f,Nx,dt1,aw0,aw0h,xx,wx) %compute right hand side for other step
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
  
%####################################################################
function  F=rightside21(ta,tb,tah,tbh,u,f,dt,Nx,alp,a0,a0h,xx,wx)
  tmp1=(1-alp)*(2-alp);
  tmph1=(1-alp/2)*(2-alp/2);
  f1=a0*f(xx,2*dt)+(2-tb-tmp1)*u(:,2)+(ta-alp/2)*u(:,1)+tmp1/2^(alp)*u(:,1);
  f2=zeros(Nx+1,1);
  f2(1)=(2-tbh-tmph1)*u(1,2)+(tah-alp/4)*u(1,1)+tmph1/2^(alp/2)*u(1,1);
  f2(Nx+1)=(2-tbh-tmph1)*u(Nx+1,2)+(tah-alp/4)*u(Nx+1,1)+tmph1/2^(alp/2)*u(Nx+1,1);
  F=diag(wx)*f1+a0/a0h*f2;

%####################################################################
function F=rightside22(ta,tb,tah,tbh,k,u,Uhist,Uhisth,Nact,Nacth,f,dt,Nx,alp,nzt,nwt,nzth,nwth,a0,a0h,xx,wx)
  tmp1=(1-alp)*(2-alp);
  f1=a0*f(xx,k*dt)+(2-tb-tmp1)*u(:,k)+(ta-alp/2)*u(:,k-1)+tmp1/k^(alp)*u(:,1);
  f2=0;
  for i=1:Nact
      f2=f2+nwt(i)*exp(-nzt(i)*dt)*Uhist(:,i);
  end
  f1=diag(wx)*(f1+alp*tmp1*dt^(alp)*f2);
  tmph1=(1-alp/2)*(2-alp/2);
  fh=zeros(Nx+1,1);
  fh(1)=(2-tbh-tmph1)*u(1,k)+(tah-alp/4)*u(1,k-1)+tmph1/k^(alp/2)*u(1,1);
  fh(Nx+1)=(2-tbh-tmph1)*u(Nx+1,k)+(tah-alp/4)*u(Nx+1,k-1)+tmph1/k^(alp/2)*u(Nx+1,1);
  th1=0;
  th2=0;
  for i=1:Nacth
      th1=th1+nwth(i)*exp(-nzth(i)*dt)*Uhisth(1,i);
      th2=th2+nwth(i)*exp(-nzth(i)*dt)*Uhisth(Nx+1,i);
  end
  fh(1)=fh(1)+alp/2*tmph1*dt^(alp/2)*th1;
  fh(Nx+1)=fh(Nx+1)+alp/2*tmph1*dt^(alp/2)*th2;
  F=f1+a0/a0h*fh;
  
%#################################################################
function [b0]=coeff_b0(alp,Nt1)
  for i=1:Nt1
    b0(i)=(i+1)^(1-alp)-i^(1-alp);%compute the coefficients b_j
  end
    
%################################################################
function [a,b,c,ta,tb,tc]=coeff_abc(alp,nzt,nwt,dt)
  tmp1=exp(-nzt*dt);
  tmp2=(dt*nzt).^2.*nzt;
  a=tmp1./tmp2-tmp1.^2./nzt-(3*tmp1.^2-tmp1)./(2*dt*nzt.^2)-tmp1.^2./tmp2;
  b=2*tmp1./tmp2-2*tmp1.^2./(dt*nzt.^2)-tmp1./nzt-2*tmp1.^2./tmp2;
  c=tmp1./tmp2-(tmp1.^2+tmp1)./(2*dt*nzt.^2)-tmp1.^2./tmp2;
  ta=alp*(1-alp)*(2-alp)*dot(nwt,a)*dt^(alp);
  tb=alp*(1-alp)*(2-alp)*dot(nwt,b)*dt^(alp);
  tc=alp*(1-alp)*(2-alp)*dot(nwt,c)*dt^(alp);
  
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