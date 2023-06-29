function test_time_fast()

%%%for parallel computation%%%%%%%%%%%
  b=[0.1 0.5 0.99];
  reps=1*10^(-11);%*desired relative error in fast algorithm
  T=1;
  %num=[10 20 25 50 100 200 250 500 1000 2000];
  num=[80];
  dt=T/num(end);
  for j=1:length(b)
  alp=b(j);
  [nzt,nwt,Nact]= SOEappr(alp,reps,dt,T);
  parfor i=1:length(num)
      %dt=T/num(i);
      %[nzt,nwt,Nact]= SOEappr(bet,reps,dt,T);
      [L2e,semiH1e,H1e,Linfty]=NIBP_FDSTFDE_3_fast(alp,nzt,nwt,Nact,T,num(i));
  end 
  end