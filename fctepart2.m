function yp=fctepart2(m,y)
mt=5.9736e24;
G=6.67e-11;
yp=zeros(4,1);
r=sqrt((y(1)^2+y(3)^2));
yp(1)=y(2);
yp(2)=-mt*G/r^3*y(1);
yp(3)=y(4);
yp(4)=-mt*G/r^3*y(3);
end