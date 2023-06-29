function paint_space(alp,dt,T)
  str=['spaceerror-' num2str(alp)  '-' num2str(dt) '-' num2str(T) '.txt'];%*
  fid = fopen(str,'rt');
  A = fscanf(fid,'%f',[6,11])';%*fscanf reads the numbers by line, then store them by column, so the first number in [] means column number of fid, the second number means line number of fid, we have a transform operation after reading, this makes A and fid arranges the elements in the same way.
  
  %rearrange A by the ascending order of polynomial degree.
  [p,q]=size(A);
  a=A(:,1);
  [A(:,1),pos]=sort(a,'ascend');
  for j=2:q
  A(:,j)=A(pos,j);
  
  semilogy(A(:,1),A(:,2),'b*-',A(:,1),A(:,3),'go-',A(:,1),A(:,4),'r+-',A(:,1),A(:,5),'yx-');
  axis([5 20 1e-12 8e-1])%* set the range of x,y axis
  %set(gca,'xlim',[6 18],'ylimmode','auto'); %*set the range of x axis, y is automatic
  set(gca, 'XTick', [6 7 8 9 10 11 12 13 15 17 19]);%* set the position on x axis where I want to mark
  %set(gca,'XTickLabel',{'7','9','11','13','15','17'});%* set the corresponding mark on the position defined above, defaultly, it's the real data
  xlabel('Polynomial degree N');
  ylabel('Error');
  legend('L2','SEMI','H1','LINF');
  str2=['spaceerror-' num2str(alp)  '-' num2str(dt) '-' num2str(T) '.eps'];%*
  saveas(gcf,str2,'psc2');%save the picture in eps, 'psc2' makes the eps image to be colorful.