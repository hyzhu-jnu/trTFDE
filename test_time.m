function test_time()

%%%for parallel computation%%%%%%%%%%%
  tic
  num=[40 80 160];
  parfor i=1:length(num)
      [L2e,semiH1e,H1e,Linfty]=FDSTFDE_3(num(i));
  end 
  toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  num=[2,5,10,20,25,50,100,200,250,500,1000];
  %corresponding dt 
  %dt=[0.5,0.2,0.1,0.05,0.04,0.02,0.01,0.005,0.004,0.002,0.001];
  alp=0.5;
  Nx=17;
  T=10;
  %M=5;N=18;n0=7;ns=12;nl=12;%*constants for the reduction
  %[zt,wt]=zwredu(alp,M,N,n0,ns,nl);% compute the original zeros and weights
  %tolybsl=10-log10(9);%* used in reduc_b1.m
  %[nzt,nwt]=reduc_b1(tolybsl,zt,wt);% compute the reducted zeros and weights  
  str=['timeerror-' num2str(alp)  '-' num2str(Nx) '-' num2str(T) '.txt'];
  %fid=fopen(str,'a');
  %s=['M=' num2str(M)  ', N=' num2str(N) ', n0=' num2str(n0) ', ns=' num2str(ns)...
   % ', nl=' num2str(nl)  ', tolybsl=' num2str(tolybsl)];  
  %fprintf(fid,s);
  %fprintf(fid,'\n');
  %fclose(fid);
  for i=1:length(num)
      Nt=num(i)*T;
      tic
      %[L2e,semiH1e,H1e,Linfty]=FDSTFDE_3_fast(alp,Nt,nzt,nwt,Nx);
      [L2e,semiH1e,H1e,Linfty]=FDSTFDE_3(alp,Nt,Nx);
      tim=toc;
      fclose('all');
      fid=fopen(str,'a');
      fprintf(fid,'%17.6f\n',tim);
      fclose(fid);
  end