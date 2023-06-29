%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is the fast algorithm matlab code for fractinal equation        %%
%%    d^alp/dt u(x,t)-d^2/dx^2 u(x,t)=f(x,t), a1<=x<=a2, t>0            %%
%%                             u(x,0)=g(x),   a1<=x<=a2                 %%
%%         d/dx u(x,t)=d^(alp/2)/dt u(x,t),   x=a1                      %%
%%         d/dx u(x,t)=-d^(alp/2)/dt u(x,t),   x=a2                     %%
%% by the method in the paper SIAM2016 of Lv&Xu(3-alp) with             %%
%% grided mesh in time discretization.                                  %% 
%% In the first time step, we just use the L1 scheme.                   %%
%% the title NS_FSDE_NZB_3fast_1D means: NS:nonsmooth solution,         %%
%% FSDE:fractional sub-diffusion equation,                              %%
%% NZB:non zero bound                                                   %%
%% 3fast:3-alp fast scheme in time, 1D: 1 dimension in space.           %%
%% Author:Hongyi Zhu                                                    %%
%% The remark begin with the %* means for different problem,            %%
%% you need to change the parameter.                                    %%
%% The remark %@: as we don't need to store all the previous solution,  %%
%% so we can define u(Nx+1,4), where u(:,1) store the initial value all %%
%% the time, when calculating the kth step, u(:,2) u(:,3) store u(:,k-1)%%
%% u(:,k), and u(:,4) store u(:,k+1).                                   %%       
%% this code is also applicable for equal mesh(smooth solution)         %%
%% only by taking gam=1.                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [num2,varargout]=NS_FSDE_NZB_3fast_1D(Nt,alp,varargin)%[fastu][L2e,Nd]= 
  fclose('all');
  %clear;
  %clc;
  format long
  global h a0 a0h aw0 aw0h wx xx Nx
  %define some parameters which need to be changed for different problem
  a1=-1; a2=0.5; %* domain of space
  T=1; %* domain of time
  %alp=0.8; %* order of fractional derivative
  Nx=24; %* degree of polynomial use for space

  %Nt=1000;%*number of step
  gam=1;%grided mesh
  h(1)=T*(1/Nt)^gam;%step size h(1)
  j=1:Nt; tt=T*(j/Nt).^gam; tt=[0 tt];
  j=j.^gam;
  h(2:Nt)=h(1).*diff(j);%step size h(2:Nt)
  
  %two useful constants
  a0=gamma(3-alp)*h(1)^(alp);
  a0h=gamma(3-alp/2)*h(1)^(alp/2);
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

  u=zeros(Nx+1,4); %@save numerical solution from Nd to Nt
  %u(:,1)=exact_u(xx,0); % initial value
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
  %u(:,1)=exp(-5*xx.^2);
  %dlmwrite('fa_exact_u.txt',u(:,1)','precision','%.15f','newline','pc');
 
  %the first time step
  t1=zeros(Nx+1,1);
  t1(1)=aw0/aw0h;t1(Nx+1)=t1(1);
  A=diag(wx)+aw0*D'*diag(wx)*D+diag(t1);%coefficient matrix for the first step 
  F=rightside1(u,f);%compute right hand side for the first iteration
  u(:,2)=pcg(A,F,1e-12,50);%*get the numerical solution at the first time step
  %dlmwrite('fa_exact_u.txt',u(:,2)','-append','precision','%.15f','newline','pc');
  
  %[L2e,semiH1e,H1e,Linfty]=err(exact_u,d_exact_u,u(:,2),a1,a2,h(1),xx);%compute error, when there exits exact solution
  %when there is no exact solution
  %exact_u=dlmread('uexact.txt');exact_u=exact_u';%size of exact_u:(Nx+1,Nt+1)
  %bs=(size(exact_u,2)-1)/Nt;
  %tmp1=u(:,2)-exact_u(:,bs+1);
  %L2e=(abs(tmp1))'*diag(wx)*(abs(tmp1));
  %e=L2e;%maximum L2 error

  %start the main scheme
  divid=10^(-4);%we always assume divid>h(1)
  Nd=min(find((h(2:end)-divid)>=0));%h(2),...,h(Nd-1)<divid; h(Nd),...,h(Nt)>=divid
   if isempty(Nd)
      disp('too small time steps, all the time steps are computed by direct method');
      return
  else
      Nd=Nd+1;
   end
  %for the previous time steps, from 1 to Nd-1, we use the direct
  %method;from Nd to Nt, we use the fast method.
  bet=alp+1; beth=alp/2+1;
  reps=10^(-12);%*desired relative error in fast algorithm
  hh=max([h(2),divid]); [nzt,nwt,Nact]= SOEappr(bet,reps,hh,T);
  [nzth,nwth,Nacth]= SOEappr(beth,reps,hh,T);
  uu=zeros(Nx+1,Nd); %@save numerical solution from 0 to Nd-1
  uu(:,1)=u(:,1); uu(:,2)=u(:,2);
  Uhist=zeros(Nx+1,Nact);Uhisth=zeros(Nx+1,Nacth);
  time=h(1);
  t2=zeros(Nx+1,1);
  t2(1)=1;t2(Nx+1)=1;
  for k=2:Nd-1
      k
      time=time+h(k);
      %[betak,xi]=coefficient_d(alp,gam,k);%direct calculation
      [betak,xi]=coeffi_gk(alp,k,tt,h,gam);%adaptive Gauss-Kronrod quadrature method
      [betakh,xih]=coeffi_gk(alp/2,k,tt,h,gam);
      A=diag(wx)+a0/betak*D'*diag(wx)*D+a0/a0h*betakh/betak*diag(t2);%coefficient matrix
      F=rightside_d(betak,betakh,xi,xih,k,uu,f,time);
      uu(:,k+1)=pcg(A,F,1e-12,50);%   
      %dlmwrite('fa_exact_u.txt',uu(:,k+1)','-append','precision','%.15f','newline','pc');
      
      [nnzt,a,b,c]=coeffk(alp,k,gam,h,nzt,nwt,tt);
      [nnzth,ah,bh,ch]=coeffk(alp/2,k,gam,h,nzth,nwth,tt);
      for i=1:Nact
          Uhist(:,i)=nnzt(i)*Uhist(:,i)+(a(i)*uu(:,k-1)-b(i)*uu(:,k)+c(i)*uu(:,k+1));
      end
      Uhisth(1,:)=nnzth'.*Uhisth(1,:)+(ah(:)*uu(1,k-1)-...
                  bh(:)*uu(1,k)+ch(:)*uu(1,k+1))';
      Uhisth(Nx+1,:)=nnzth'.*Uhisth(Nx+1,:)+(ah(:)*uu(Nx+1,k-1)-...
                     bh(:)*uu(Nx+1,k)+ch(:)*uu(Nx+1,k+1))';
                     
      %[L2e,semiH1e,H1e,Linfty]=err(exact_u,d_exact_u,uu(:,k+1),a1,a2,h(1),xx);%compute         error, when there exits exact solution
      %when there is no exact solution
      %tmp1=uu(:,k+1)-exact_u(:,k*bs+1);
      %L2e=(abs(tmp1))'*diag(wx)*(abs(tmp1));
      %e=max([L2e,e]);%maximum L2 error
  end
  u(:,2)=uu(:,Nd-1); u(:,3)=uu(:,Nd);
  num2(:,1:Nd)=uu;
  clear uu betak xi betakh xih
  
  for k=Nd:Nt
      k
      time=time+h(k);
      [nnzt,a,b,c,yi,de,be,ta,tb,tc,cc]=coeffk(alp,k,gam,h,nzt,nwt,tt);
      [nnzth,ah,bh,ch,yih,deh,beh,tah,tbh,tch,cch]=coeffk(alp/2,k,gam,h,nzth,nwth,tt);
      A=(yi-tc)*diag(wx)+a0*D'*diag(wx)*D+a0/a0h*(yih-tch)*diag(t2);%coefficient matrix
      F=rightside_f(time,alp,k,gam,Uhist,cc,de,tb,ta,be,Uhisth,cch,deh,tbh,tah,beh,u,f);
      u(:,4)=pcg(A,F,1e-12,50);%@
      num2(:,k+1)=u(:,4);
      %dlmwrite('fa_exact_u.txt',u(:,4)','-append','precision','%.15f','newline','pc');
      for i=1:Nact
          Uhist(:,i)=nnzt(i)*Uhist(:,i)+(a(i)*u(:,2)-b(i)*u(:,3)+c(i)*u(:,4));
      end
      Uhisth(1,:)=nnzth'.*Uhisth(1,:)+(ah(:)*u(1,2)-...
                  bh(:)*u(1,3)+ch(:)*u(1,4))';
      Uhisth(Nx+1,:)=nnzth'.*Uhisth(Nx+1,:)+(ah(:)*u(Nx+1,2)-...
                     bh(:)*u(Nx+1,3)+ch(:)*u(Nx+1,4))';
                            
     %e=max([abs([exact_u(xx,time)-u(:,4)]'),e]);%e is the maximum nodal error 
     %[L2e,semiH1e,H1e,Linfty]=err(exact_u,d_exact_u,uu(:,k+1),a1,a2,h(1),xx);%compute         error, when there exits exact solution
      %when there is no exact solution
      %tmp1=u(:,4)-exact_u(:,k*bs+1);
      %L2e=(abs(tmp1))'*diag(wx)*(abs(tmp1));
      %e=max([L2e,e]);%maximum L2 error
    
      u(:,2)=u(:,3);%@
      u(:,3)=u(:,4);%@
  end
  %dlmwrite('fa_maxL2_error.txt', [Nt, e], '-append', 'precision', '%7d %19.15f', 'newline', 'pc');
  if nargin==2
      return
  else
      num=varargin{1};
      dif=num-num2(:,1:2:end);
      err=diag(dif'*diag(wx)*dif);
      err=max(sqrt(err));
      varargout{1}=err;
  end
  return
  %*plot the solution
  %fplot(exact_u,[a1,a2])
  %hold on
  %fplot(numeri_u,[a1,a2])
  %legend('exact solution','numerical solution');
  
  %*save the time errors
  str=['timeerror-' num2str(alp)  '-' num2str(Nx) '-' num2str(T) '.txt'];
  %str=['ttt.txt'];
  fid=fopen(str,'a');
  t3=[Nt e gam];
  fprintf(fid,'%7d  %19.17e %5.1f\n',t3);
  fclose(fid);

  %*save the space errors
  %str=['spaceerror-' num2str(alp)  '-' num2str(dt) '-' num2str(T) '.txt'];
  %fid = fopen(str,'a');
  %t3=[Nx L2e semiH1e H1e Linfty];
  %fprintf(fid,'%4d  %20.17e  %20.17e  %20.17e  %20.17e\n',t3);
  %fclose(fid);
end
 
%##################################################################
function F=rightside1(u,f) % compute right hand side for the fist iteration
  global h aw0 wx xx
  f1=f(xx,h(1));
  F=diag(wx)*(aw0*f1+u(:,1));
end

%##################################################################
function F=rightside_d(betak,betakh,xi,xih,k,uu,f,time)
  global wx xx Nx a0 a0h 
  f1=0;
  for i=1:k
      f1=f1+xi(i)*uu(:,i);
  end
  F=diag(wx)*(f1+a0/betak*f(xx,time));
  f2=zeros(Nx+1,1);
  for i=1:k
      f2(1)=f2(1)+xih(i)*uu(1,i);
      f2(Nx+1)=f2(Nx+1)+xih(i)*uu(Nx+1,i);
  end
  F=F+a0/a0h/betak*betakh*f2;
end
 
%##################################################################
function F=rightside_f(time,alp,k,gam,Uhist,cc,de,tb,ta,be,Uhisth,cch,deh,tbh,tah,beh,u,f)
global wx  xx Nx a0 a0h
  f1=0;
  f1=(1-alp)*(2-alp)/k^(alp*gam)*u(:,1)+Uhist*cc;
  f1=f1+(de-tb)*u(:,3)+(ta-be)*u(:,2)+a0*f(xx,time);
  F=diag(wx)*f1;
  fh=zeros(Nx+1,1);
  fh(1)=(1-alp/2)*(2-alp/2)/k^(alp/2*gam)*u(1,1)+(deh-tbh)*u(1,3)+(tah-beh)*u(1,2)+dot(Uhisth(1,:),cch');
  fh(Nx+1)=(1-alp/2)*(2-alp/2)/k^(alp/2*gam)*u(Nx+1,1)+(deh-tbh)*u(Nx+1,3)+(tah-beh)*u(Nx+1,2)+dot(Uhisth(Nx+1,:),cch');
  F=F+a0/a0h*fh;
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

%###################################################################
function [nnzt,a,b,c,varargout]=coeffk(alp,k,gam,h,nzt,nwt,tt)%[nnzt,a,b,c,yi,de,be,ta,tb,tc,cc]
  tm=(1-alp)*(2-alp)*h(1)^alp;
  %two methods for calculating yi,de,be,a,b,c
  %method 1: direct calculation
  t=[k-2 k-1 k];t=t.^gam;t=diff(t);%t(1)=(k-1)^gam-(k-2)^gam, t(2)=k^gam-(k-1)^gam
  yi=(2*t(2)+(2-alp)*t(1))/(t(2)^alp*sum(t));
  de=alp*((t(2)+(2-alp)*t(1))/(t(2)^alp*t(1)));
  be=alp*t(2)^(2-alp)/(t(1)*sum(t));
  %a=exp(-nzt*h(k)-log(h(k-1)*(h(k-1)+h(k))*nzt.^2)).*(h(k)+2./nzt-exp(-nzt*h(k-1)+log(nzt*h(k-1)*(h(k-1)+h(k))+2*h(k-1)+h(k)+2./nzt)));
  %c=exp(-nzt*h(k)-log(h(k)*(h(k-1)+h(k))*nzt.^2)).*(-h(k-1)+2./nzt-exp(-nzt*h(k-1)+log(h(k-1)+2./nzt)));
  %b=a+c-exp(-nzt*h(k)-log(nzt)).*(1-exp(-nzt*h(k-1)));
 
  %method 2: adaptive Gauss-Kronrod quadrature
  %eps=min([10^(-7),(tt(k+1)-tt(k))/10^4]);
  %yi=quadgk(@(s)tm*(tt(k+1)-s).^(-alp).*((2*s-tt(k-1)-tt(k))./(h(k)*(h(k-1)+h(k)))),tt(k),tt(k+1)-eps);
  %de=quadgk(@(s)tm*(tt(k+1)-s).^(-alp).*((2*s-tt(k-1)-tt(k+1))./(h(k-1)*h(k))),tt(k),tt(k+1)-eps);
  %de=de-(1-alp)*(2-alp)/(k^gam-(k-1)^gam)^alp;
  %be=quadgk(@(s)tm*(tt(k+1)-s).^(-alp).*((2*s-tt(k)-tt(k+1))./(h(k-1)*(h(k-1)+h(k)))),tt(k),tt(k+1)-eps);
  for i=1:length(nzt)
      a(i)=quadgk(@(s)exp(-nzt(i)*(tt(k+1)-s)).*((s-tt(k))/h(k-1).*((s-tt(k+1))/(h(k-1)+h(k)))),tt(k-1),tt(k));
      b(i)=quadgk(@(s)exp(-nzt(i)*(tt(k+1)-s)).*((s-tt(k-1))/h(k-1).*((s-tt(k+1))/h(k))),tt(k-1),tt(k));
      c(i)=quadgk(@(s)exp(-nzt(i)*(tt(k+1)-s)).*((s-tt(k-1))/h(k).*((s-tt(k))/(h(k-1)+h(k)))),tt(k-1),tt(k));
  end
  
  ta=alp*tm*dot(nwt,a);
  tb=alp*tm*dot(nwt,b);
  tc=alp*tm*dot(nwt,c);
  nnzt=exp(-nzt*h(k));
  cc=alp*tm*nwt.*nnzt;
  varargout{1}=yi; varargout{2}=de; varargout{3}=be; varargout{4}=ta; varargout{5}=tb; varargout{6}=tc; varargout{7}=cc;
end