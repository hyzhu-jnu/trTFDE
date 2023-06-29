function paint_diff_u_fastu()
  alp=0.5;
  T=10; Nt=1000; dt=T/Nt;
  Nx=24;
  [xx,wx]= legslb(Nx+1);
  
  reps=1*10^(-10);%*desired relative error in fast algorithm
  [nzt,nwt,Nact]= SOEappr(alp,reps,dt,T);
  [fastu]=NIBP_FDSTFDE_3_fast_nosub(alp,T,Nt,Nx,xx,wx,nzt,nwt,Nact)

  [u]=FDSTFDE_3_nosub(alp,T,Nt,Nx,xx,wx);

  t=linspace(0,T,Nt+1);
  xx=pi/2*xx+pi/2;%*transform xx from [-1,1] to the compute domain 

  %compute fastu-u directly, where fastu,u are two-dimentional vectors, (Nx+1)*(Nt+1)
  diff=fastu-u;
  norm(diff,Inf)

  %write fasu,u into lagrange functions 
  %sx=linspace(0,pi,3000);
  %for j=1:Nt+1
  %    diff(:,j) = Lagrange_1(xx,fastu(:,j)-u(:,j),sx);
  %end
  %max(max(abs(diff)))%find the maximun value in abs(diff)
  %min(min(abs(diff)))%find the minimun value in abs(diff)
  
  %surf(t,sx,diff);%*
  %shading interp;
  %xlabel('t');
  %ylabel('x');
  %zlabel('Error');
  %str2=['6fastu-u-' num2str(alp*100)  '-' num2str(Nx) '-' num2str(T) '-' num2str(Nt) '.eps'];%*
  %saveas(gcf,str2,'psc2');%save the picture in eps, 'psc2' makes the eps image to be colorful.