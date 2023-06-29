function  [varargout]=coeffi_gk(alp,k,tt,h,gam)%adaptive Gauss-Kronrod quadrature [betak,xi]=coeffi_gk(alp,k,tt,h,gam)
  format long 
  tm=(1-alp)*(2-alp)*h(1)^alp;
  %direct method for computing yi,de be
  t=[k-2 k-1 k];t=t.^gam;t=diff(t);%t(1)=(k-1)^gam-(k-2)^gam, t(2)=k^gam-(k-1)^gam
  yi=(2*t(2)+(2-alp)*t(1))/(t(2)^alp*sum(t));
  de=((2-alp)*t(1)+alp*t(2))/(t(1)*t(2)^alp);
  be=alp*t(2)^(2-alp)/(t(1)*sum(t));
  %Gauss-Kronrod quadrature for computing yi,de,be
  %eps=min([10^(-12),(tt(k+1)-tt(k))/10^4]);%for the singularity at tt(k+1)
  %yi=quadgk(@(s)tm*(tt(k+1)-s).^(-alp).*((2*s-tt(k-1)-tt(k))./(h(k)*(h(k-1)+h(k)))),tt(k),tt(k+1)-eps);
  %de=quadgk(@(s)tm*(tt(k+1)-s).^(-alp).*((2*s-tt(k-1)-tt(k+1))./(h(k-1)*h(k))),tt(k),tt(k+1)-eps);
  %be=quadgk(@(s)tm*(tt(k+1)-s).^(-alp).*((2*s-tt(k)-tt(k+1))./(h(k-1)*(h(k-1)+h(k)))),tt(k),tt(k+1)-eps);
  for i=1:k-1
      a(i)=quadgk(@(s)tm*(tt(k+1)-s).^(-alp).*(2*s-tt(k-i+1)-tt(k-i+2))/h(k-i)/(h(k-i)+h(k-i+1)),tt(k-i),tt(k-i+1));
      b(i)=quadgk(@(s)(-1)*tm*(tt(k+1)-s).^(-alp).*(2*s-tt(k-i)-tt(k-i+2))/h(k-i)/h(k-i+1),tt(k-i),tt(k-i+1));
  end
  c=-a-b;
  betak=c(1)+yi;
  varargout{1}=betak;
 if nargout==1, return; end;
 
  if(k==2)
      xi(1)=-a(1)-be;
      xi(2)=-b(1)+de;
  elseif(k==3)
      xi(1)=-a(2);
      xi(2)=-a(1)-b(2)-be;
      xi(3)=-b(1)-c(2)+de;
  else
      xi(1)=-a(k-1);
      xi(2)=-a(k-2)-b(k-1);
      for i=3:k-2
          xi(i)=-a(k-i)-b(k-i+1)-c(k-i+2);
      end
      xi(k-1)=-a(1)-b(2)-c(3)-be;
      xi(k)=-b(1)-c(2)+de;
  end
  xi=xi/betak;
  varargout{2}=xi;
end