function [betan,xi]=coefficient_d(alp,gam,n)%direct calculation
  %clear;clc;
  %format long 
%   alp=0.9;
%   N=300;
%   gam=5;%(3-alp)/alp;

  t=1:n;t=t.^gam;xi=0;
  if n==2
      betan=(2*t(2)-2+alp)/((t(2)-1)*t(2)^alp);
      xi(1)=(t(2)-1)*(2-alp-alp*t(2))/(2*t(2)-2+alp);
      xi(2)=alp*t(2)^2/(2*t(2)-2+alp);
  elseif n==3
          betan=(2*(t(3)-t(2))+alp*(t(2)-1))/((t(3)-t(2))*(t(3)-1)^alp);
          xi(1)=(t(3)-t(2))*(t(3)-1)/t(2)-(t(3)-t(2))*(t(3)-1)^alp*t(3)^(1-alp)/t(2)*(2*t(3)-(2-alp)*t(2)-(2-alp))/(2*t(3)-(2-alp)*t(2)-alp);
          xi(2)=(t(3)-t(2))*(t(3)-1)^alp*t(3)^(1-alp)*(2*t(3)-(2-alp)*t(2))/((t(2)-1)*(2*(t(3)-t(2))+alp*(t(2)-1)))-(t(3)-t(2))*t(3)/(t(2)-1);
          xi(3)=(t(3)-1)*t(3)/((t(2)-1)*t(2))-(t(3)-1)^alp*(t(3)-t(2))*(alp-2+2*t(3))*t(3)^(1-alp)/(t(2)*(t(2)-1)*(2*(t(3)-t(2))+alp*(t(2)-1)));
  else
  a=0;b=0;c=0;
  %for n=4:N
     betan=(2*t(n)-(2-alp)*t(n-1)-alp*t(n-2))/((t(n)-t(n-1))*(t(n)-t(n-2))^alp);
     
     for j=2:n-1
         a(n-j)=(2*t(n)+2*(1-alp)*t(j-1)-(2-alp)*(t(j)+t(j+1)))*(t(n)-t(j-1))^(1-alp)...
               +(-2*t(n)+alp*t(j)+(2-alp)*t(j+1))*(t(n)-t(j))^(1-alp);
         a(n-j)=a(n-j)/((t(j)-t(j-1))*(t(j+1)-t(j-1)));
         c(n-j)=(2*t(n)-alp*t(j-1)-(2-alp)*t(j))*(t(n)-t(j-1))^(1-alp)...
               +(-2*t(n)+alp*t(j)+(2-alp)*t(j-1))*(t(n)-t(j))^(1-alp);
         c(n-j)=c(n-j)/((t(j+1)-t(j))*(t(j+1)-t(j-1)));
         b(n-j)=-a(n-j)-c(n-j);
     end
     a(n-1)=(2*t(n)-(2-alp)*(1+t(2)))*t(n)^(1-alp)+(-2*t(n)+alp+(2-alp)*t(2))*(t(n)-1)^(1-alp);
     a(n-1)=a(n-1)/t(2);
     c(n-1)=(2*t(n)-(2-alp))*t(n)^(1-alp)+(-2*t(n)+alp)*(t(n)-1)^(1-alp);
     c(n-1)=c(n-1)/((t(2)-1)*t(2));
     b(n-1)=-a(n-1)-c(n-1);
     
     xi(1)=-a(n-1);
     xi(2)=-a(n-2)-b(n-1);
     for j=3:n-2
         xi(j)=-a(n-j)-b(n-j+1)-c(n-j+2);
     end
     xi(n-1)=-a(1)-b(2)-c(3)-alp*(t(n)-t(n-1))^(2-alp)/((t(n-1)-t(n-2))*(t(n)-t(n-2)));
     xi(n)=-b(1)-c(2)+((2-alp)*(t(n-1)-t(n-2))+alp*(t(n)-t(n-1)))/((t(n-1)-t(n-2))*(t(n)-t(n-1))^alp);
     
     
%      xi(1)=(2*t(n)-alp-(2-alp)*t(2))*(t(n)-1)^(1-alp)-(2*t(n)-(2-alp)*(1+t(2)))*t(n)^(1-alp);
%      xi(1)=xi(1)/t(2);
%      
%      xi(2)=(t(3)-1)*(2*t(n)-(2-alp)*t(2))*t(n)^(1-alp)+(2*t(n)-alp*t(2)-(2-alp)*t(3))*(t(n)-t(2))^(1-alp)...
%           -t(3)*(2*t(n)-(2-alp)*t(2)-alp)*(t(n)-1)^(1-alp);
%      xi(2)=xi(2)/((t(2)-1)*(t(3)-1));
%      
%      for j=3:n-2
%          if(j==3)
%              xi(3)=(2*t(n)-alp*t(3)-(2-alp)*t(4))/((t(3)-t(2))*(t(4)-t(2)))*(t(n)-t(3))^(1-alp)...
%                   -((t(4)-1)*(2*t(n)-alp*t(2)-(2-alp)*t(3)))/((t(2)-1)*(t(3)-t(2))*(t(4)-t(2)))*(t(n)-t(2))^(1-alp)...
%                   +(t(3)*(2*t(n)-alp-(2-alp)*t(2)))/(t(2)*(t(2)-1)*(t(3)-t(2)))*(t(n)-1)^(1-alp)...
%                   -(2*t(n)-(2-alp))/(t(2)*(t(2)-1))*t(n)^(1-alp);
%          else
%              xi(j)=(2*t(n)-alp*t(j)-(2-alp)*t(j+1))/((t(j)-t(j-1))*(t(j+1)-t(j-1)))*(t(n)-t(j))^(1-alp)...
%                   -((t(j+1)-t(j-2))*(2*t(n)-alp*t(j-1)-(2-alp)*t(j)))/((t(j-1)-t(j-2))*(t(j)-t(j-1))*(t(j+1)-t(j-1)))*(t(n)-t(j-1))^(1-alp)...
%                   +((t(j)-t(j-3))*(2*t(n)-alp*t(j-2)-(2-alp)*t(j-1)))/((t(j-1)-t(j-3))*(t(j-1)-t(j-2))*(t(j)-t(j-1)))*(t(n)-t(j-2))^(1-alp)...
%                   -(2*t(n)-alp*t(j-3)-(2-alp)*t(j-2))/((t(j-1)-t(j-3))*(t(j-1)-t(j-2)))*(t(n)-t(j-3))^(1-alp);
%          end
%      end
%      
%      if(n==4)
%          xi(n-1)=(t(3)*(2*t(4)-alp-(2-alp)*t(2)))/(t(2)*(t(2)-1)*(t(3)-t(2)))*(t(4)-1)^(1-alp)...
%             -((t(4)-1)*(2*t(4)-alp*t(2)-(2-alp)*t(3)))/((t(2)-1)*(t(3)-t(2))*(t(4)-t(2)))*(t(4)-t(2))^(1-alp)...
%             -(2*t(4)-(2-alp))/(t(2)*(t(2)-1))*t(4)^(1-alp);
%      else
%          xi(n-1)=((t(n-1)-t(n-4))*(2*t(n)-alp*t(n-3)-(2-alp)*t(n-2)))/((t(n-2)-t(n-4))*(t(n-2)-t(n-3))*(t(n-1)-t(n-2)))*(t(n)-t(n-3))^(1-alp)...
%                 -((t(n)-t(n-3))*(2*t(n)-alp*t(n-2)-(2-alp)*t(n-1)))/((t(n-2)-t(n-3))*(t(n-1)-t(n-2))*(t(n)-t(n-2)))*(t(n)-t(n-2))^(1-alp)...
%                 -(2*t(n)-alp*t(n-4)-(2-alp)*t(n-3))/((t(n-2)-t(n-4))*(t(n-2)-t(n-3)))*(t(n)-t(n-4))^(1-alp);
%      end
%      
%      xi(n)=((t(n)-t(n-3))*(2*t(n)-(2-alp)*t(n-1)-alp*t(n-2)))/((t(n-1)-t(n-2))*(t(n)-t(n-1))*(t(n-1)-t(n-3)))*(t(n)-t(n-2))^(1-alp)...
%           -(2*t(n)-alp*t(n-3)-(2-alp)*t(n-2))/((t(n-1)-t(n-2))*(t(n-1)-t(n-3)))*(t(n)-t(n-3))^(1-alp);
       
     xi=xi./betan;
  end
     return
     bn(n-3)=betan;
     b3(n-3)=xi(n-2);
     b2(n-3)=xi(n-1);
     b1(n-3)=xi(n);
     s1(n-3)=xi(n)^2+4*xi(n-1);
     s2(n-3)=xi(n)+xi(n-1);
  %end
  find(b2<0)
  find(s1<0)
  find(s2<0)
  %for n=2:N-3
      %k1(n-3)=b1(n-1);
      %k2(n-3)=-b2(n-1)/b1(n-2);
      %k3(n-3)=-b3(n-1)/b2(n-2);
      %kk(n-1)=b1(n-1)^2+4*b2(n);
  %end
  %find(kk<0)
  
%   ga=find(s1<0)
%   k1(ga)
%   k2(ga)
%   k3(ga)
      
  
  subplot(2,2,1)
  plot(4:N,b2)
  legend('\xi^n_{n-2}')
  
  subplot(2,2,2)
  plot(4:N,b1)
  legend('\xi^n_{n-1}')
  
  subplot(2,2,3)
  plot(4:N,s1)
  legend('(\xi^n_{n-1})^2+4\xi^n_{n-2}')
  
  subplot(2,2,4)
  plot(4:N,s2)
%   plot(1:N-6,k1,'b')
%   hold on
%   plot(1:N-6,k2,'r')
%   hold on
%   plot(1:N-6,k3,'y')
%   legend('k1','k2','k3')
 
end