function paint_time(alp,Nx,T)
  str=['timeerror-' num2str(alp)  '-' num2str(Nx) '-' num2str(T) '.txt'];%*
  fid = fopen(str,'rt');
  A = fscanf(fid,'%f',[7,6])';%*fscanf reads the numbers by line, then store them by column, so the first number in [] means column number of fid, the second number means line number of fid, we have a transform operation after reading, this makes A and fid arranges the elements in the same way.
  fclose(fid);
  [p,q]=size(A);

  %rearrange A by the descending order of time step.
  a=A(:,1);
  [A(:,1),pos]=sort(a,'descend');
  for j=2:q
  A(:,j)=A(pos,j);

  %ratio store the accuracy, it has 4 columns, each means L2e,semiH1e,H1e,Linfty accuracy.
  [p,q]=size(A);
  for j=1:4
      for i=1:p-1
          ratio(i,j)=(log10(A(i+1,j+2))-log10(A(i,j+2)))/(log10(A(i+1,1))-log10(A(i,1)));%*the codes in yellow label should change according to the order of four norms in A.
      end
  end
  %save ratio
  str1=['ratio-' num2str(alp)  '-' num2str(Nx) '-' num2str(T) '.txt'];%*
  fid = fopen(str1,'wt');
  for i=1:p-1
      t3=ratio(i,:);
      fprintf(fid,'%20.17f  %20.17f  %20.17f  %20.17f\n',t3);
  end
  fclose(fid);
  
  loglog(A(:,1),A(:,3),'b*-',A(:,1),A(:,4),'go-',A(:,1),A(:,5),'r+-',A(:,1),A(:,6),'yx-',A(:,1),A(:,7),'kd-');
  axis([0.004 0.11 1e-7 3e-1])%* set the range of x,y axis
  %set(gca,'xlim',[0.0008 0.11],'ylimmode','auto'); %*set the range of x axis, y is automatic
  set(gca, 'XTick', [0.005  0.01  0.1]);%* set the position on x axis where I want to mark
  set(gca,'XTickLabel',{'0.005','0.01','0.1'});%* set the corresponding mark on the position defined above
  xlabel('Time step');
  ylabel('Error');
  legend('L2','SEMI','H1','LINF',num2str(3-alp),'Location','SouthEast');
  str2=['timeerror-' num2str(alp)  '-' num2str(Nx) '-' num2str(T) '.eps'];%*
  saveas(gcf,str2,'psc2');%save the picture in eps, 'psc2' makes the eps image to be colorful.