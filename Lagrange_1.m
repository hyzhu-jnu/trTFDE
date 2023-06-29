function y = Lagrange_1(x0,f0,x)
%实现拉格朗日多项式插值公式，x0,f0为已知节点向量，x为插值位置，可以为向量 
  L = length(x0);
  x0 = repmat(reshape(x0,1,L),L,1);
  f0 =  repmat(reshape(f0,1,L),L,1);

  x1 =  diag(x0)*ones(1,L) - x0+ eye(L);
  x1 = prod(x1');

  for i = 1 : length(x)
      %temp = abs(x(i) == x0(1,:))<1e-6%如果插值点与节点距离小于1e-6时，取节点值
      %if sum(temp)
       %   y(i) = f0(1,temp);
      %else
          y(i) = sum(prod(x(i)-x0(1,:))./(x1.*(x(i)-x0(1,:))).*f0(1,:));
      %end
  end