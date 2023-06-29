function test_time_sub()%this function is test the time error with sub-scheme(FDSTFDE_3_fast.m) or without sub-scheme(FDSTFDE_3_fast_nosub.m) at the first time step
  Nx=24;
  alp=0.5;
  reps=1*10^(-10);%*desired relative error in fast algorithm
  T=1;
  num=[10 20 40 80 160];
  dt=T/num(end);
  [nzt,nwt,Nact]= SOEappr(alp,reps,dt,T);
  parfor i=1:length(num)
      %dt=T/num(i);
      %[nzt,nwt,Nact]= SOEappr(bet,reps,dt,T);
      [L2e,semiH1e,H1e,Linfty]=NIBP_FDSTFDE_3_fast_nosub(Nx,alp,nzt,nwt,Nact,T,num(i));
      [L2e,semiH1e,H1e,Linfty]=NIBP_FDSTFDE_3_fast(Nx,alp,nzt,nwt,Nact,T,num(i));
  end 



%%%%%%%%%%%%%%%%%%%%%%compute ratio with sub-stepping  
  str=['sub-fast-timeerror-' num2str(alp)  '-' num2str(Nx) '-' num2str(T) '.txt'];%*
  fid = fopen(str,'rt');
  A = fscanf(fid,'%f',[7,length(num)])';%*fscanf reads the numbers by line, then store them by column, so the first number in [] means column number of fid, the second number means line number of fid, we have a transform operation after reading, this makes A and fid arranges the elements in the same way.
  fclose(fid);
  [p,q]=size(A);

  %rearrange A by the descending order of time step.
  a=A(:,1);
  [A(:,1),pos]=sort(a,'descend');
  for j=2:q
  A(:,j)=A(pos,j);
  end

  %ratio store the accuracy, it has 4 columns, each means L2e,semiH1e,H1e,Linfty accuracy.
  [p,q]=size(A);
  for j=1:4
      for i=1:p-1
          ratio(i,j)=log2(A(i,j+2))-log2(A(i+1,j+2));%*the codes in yellow label should change according to the order of four norms in A.
      end
  end
  %save ratio
  str1=['subratio-' num2str(alp)  '-' num2str(Nx) '-' num2str(T) '.txt'];%*
  fid = fopen(str1,'wt');
  for i=1:p-1
      t3=ratio(i,:);
      fprintf(fid,'%20.17f  %20.17f  %20.17f  %20.17f\n',t3);
  end
  fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%compute ratio without sub-stepping
  str=['nosub-fast-timeerror-' num2str(alp)  '-' num2str(Nx) '-' num2str(T) '.txt'];%*
  fid = fopen(str,'rt');
  A = fscanf(fid,'%f',[7,length(num)])';%*fscanf reads the numbers by line, then store them by column, so the first number in [] means column number of fid, the second number means line number of fid, we have a transform operation after reading, this makes A and fid arranges the elements in the same way.
  fclose(fid);
  [p,q]=size(A);

  %rearrange A by the descending order of time step.
  a=A(:,1);
  [A(:,1),pos]=sort(a,'descend');
  for j=2:q
  A(:,j)=A(pos,j);
  end

  %ratio store the accuracy, it has 4 columns, each means L2e,semiH1e,H1e,Linfty accuracy.
  [p,q]=size(A);
  for j=1:4
      for i=1:p-1
          ratio(i,j)=log2(A(i,j+2))-log2(A(i+1,j+2));%*the codes in yellow label should change according to the order of four norms in A.
      end
  end
  %save ratio
  str1=['nosubratio-' num2str(alp)  '-' num2str(Nx) '-' num2str(T) '.txt'];%*
  fid = fopen(str1,'wt');
  for i=1:p-1
      t3=ratio(i,:);
      fprintf(fid,'%20.17f  %20.17f  %20.17f  %20.17f\n',t3);
  end
  fclose(fid);