function test_space_fast()
%%%for parallel computation%%%%%%%%%%%
  tic
  %alp=0.99;
  %bet=alp+1;
  reps=1*10^(-9);%*desired relative error in fast algorithm
  T=1;Nt=10000; 
  a=6:15;
  b=[0.1 0.5 0.99];
  for j=1:length(b)
  alp=b(j); 
  [nzt,nwt,Nact]= SOEappr(alp,reps,T/Nt,T);
  parfor i=1:length(a)
  [L2e,semiH1e,H1e,Linfty]=NIBP_FDSTFDE_3_fast(alp,nzt,nwt,Nact,T,Nt,a(i));
  end
  end 
  toc