function test_space()
%%%for parallel computation
  tic
  a=6:15;
  b=[0.1 0.5 0.99];
  for j=1:length(b)
  alp=b(j);
  parfor i=1:length(a)
      [L2e,semiH1e,H1e,Linfty]=FDSTFDE_3(alp,a(i));
  end 
  end
  toc
%%%%%%%%%%%%%%%%%%%%%%%%%% 
  T=1;
  Nt=10000;
  dt=T/Nt;
  alp=0.2;
  
  str=['spaceerror-' num2str(alp)  '-' num2str(dt) '-' num2str(T) '.txt'];
  fid=fopen(str,'a');
  for Nx=6:17
      tic
      [L2e,semiH1e,H1e,Linfty]=FDSTFDE_3(Nx);
      tim=toc;
      fclose('all');
      fid=fopen(str,'a');
      fprintf(fid,'%17.6f\n',tim);
      fclose(fid);
  end

