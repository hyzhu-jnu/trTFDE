function [a]=cpu_time()
  num=[10000 20000 40000 80000 160000];
  %num=[1000];
  alp=0.5;  T=1; reps=10^(-10);
  Nx=9;
  [xx,wx]= legslb(Nx+1);%Legendre-Gauss-Lobatto points and weights in [-1,1]
  D=legslbdiff(Nx+1,xx);%the first-order differentiation matrix in [-1,1]
  xx=pi/2*xx+pi/2; %Legendre-Gauss-Lobatto for [a1,a2]
  wx=pi/2*wx;%weights in [a1,a2]

  for i=1:length(num)
    tic
    [nzt,nwt,Nact]= SOEappr(alp,reps,T/num(i),T);
    [L2e,semiH1e,H1e,Linfty]=NIBP_FDSTFDE_3_fast_nosub(alp,Nx,T,num(i),nzt,nwt,Nact,wx,xx,D);
    a(i,1)=toc;
  end
  for i=1:length(num)    
    tic
    [L2e,semiH1e,H1e,Linfty]=FDSTFDE_3_nosub(alp,Nx,T,num(i),wx,xx,D);
    a(i,2)=toc;
  end
  
a=[10000   1.871627283994016e+01     4.084231022204660e+01
   20000   3.073307879357048e+01     1.334713048327487e+02
   40000   5.455079932841189e+01     4.757626207160812e+02
   80000   1.061853004142258e+02     1.645314419472292e+03
   160000  2.031883555756268e+02     6.131028578845317e+03];
a=[10000    3.537922940007715e+01     9.300875067462560e+01
   20000    7.022062887508800e+01     2.800316597201255e+02
   40000    1.395893020415390e+02     9.077502337134046e+02
   80000    3.462930308501694e+02     3.231886932747845e+03
   160000   6.683736621885496e+02     1.214298254868637e+04];
a=log10(a);
loglog(a(:,1),a(:,2),'b*-',a(:,1),a(:,3),'go-');
axis([4 5.21 1.5 4.2]);
set(gca, 'XTick', [4 4.3 4.6 4.9 5.2]);
xlabel('Log_{10}(N_T)');
ylabel('Log_{10}(cputime)');
legend('Fast','Direct','Location','SouthEast');
grid on;
%str2=['cputime.eps'];
%set(gca,'position',[0 0 0 0]);
%saveas(gcf,str2,'psc2','format');%save the picture in eps, 'psc2' makes the eps image to be colorful.
%print('-depsc',str2);
