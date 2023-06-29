function [a,N_varepsilon]=varepsilon_vs_time
  alp=0.5;bet=alp+1;
  T=1;Nt=10000;dt=1/10000;
  reps=[10^(-5),10^(-6),10^(-7),10^(-8),10^(-9),10^(-10)];
  for i=1:length(reps)
     tic
     [nzt,nwt,N_varepsilon(i)]= SOEappr(bet,reps(i),dt,T);
     [L2e,semiH1e,H1e,Linfty]=FDSTFDE_3_fast(alp,dt,T,Nt,nzt,nwt,N_varepsilon(i));
     a(i)=toc;
  end
  return

%%%%%%%%%%%%%%%%%%%%%%%%
%result:
reps=[10^(-5),10^(-6),10^(-7),10^(-8),10^(-9),10^(-10)];
reps=log10(reps);
a=[2.318538073393122e+01     2.257870094969121e+01     2.308706222664892e+01... 
  2.277062478964770e+01     2.315569855354591e+01     2.363953301895063e+01];
N_varepsilon=[35    41    46    51    57    62];
plot(reps,a,'b*-',reps,N_varepsilon,'go-')


reps=[10^(-5),10^(-6),10^(-7),10^(-8),10^(-9),10^(-10)];
reps=log10(reps);
a=[2.717162053919470e+01     2.785340705517891e+01     2.750152907683616e+01...
   2.750936626383086e+01     2.828996351683824e+01     2.738844239436626e+01];
N_varepsilon=[35    41    46    51    57    62];
plot(reps,a,'b*-',reps,N_varepsilon,'go-')
