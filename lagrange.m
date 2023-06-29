%{
������Ϊlagrange��ֵ�����������й���Ҫ��������������
tΪ��ֵ���ֵ��yΪ��ֵ���ĺ���ֵ��
LA��lagrange��ֵ�����ע��ú����Ա���Ϊx
%}
function LA=lagrange(t,y,x)
m=length(t);  %m,n�ֱ�������ʾ������ͺ���ֵ����
n=length(y);
if m~=n  %�ڵ����ͺ���ֵ������ͬʱ����
    disp('���󣺽�����Ӧ���뺯��ֵ������ͬ');
    return
end
s=0;
for i=1:m %lagrange��ֵ���쿪ʼ
    l=1;
    for j=1:m
        if j~=i
            l=l.*(x-t(j))/(t(i)-t(j));
        end
    end
    s=s+l*y(i);
end %lagrange��ֵ�������
LA=s;
