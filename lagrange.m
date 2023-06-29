%{
本函数为lagrange插值基函数，其中共需要输入两个参数：
t为插值结点值，y为插值结点的函数值。
LA是lagrange插值结果，注意该函数自变量为x
%}
function LA=lagrange(t,y,x)
m=length(t);  %m,n分别用来表示结点数和函数值个数
n=length(y);
if m~=n  %节点数和函数值个数不同时处理
    disp('错误：结点个数应该与函数值个数相同');
    return
end
s=0;
for i=1:m %lagrange插值构造开始
    l=1;
    for j=1:m
        if j~=i
            l=l.*(x-t(j))/(t(i)-t(j));
        end
    end
    s=s+l*y(i);
end %lagrange插值构造结束
LA=s;
