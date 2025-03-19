function yp=fcte(t,y)
esp=2;
yp=zeros(size(y));
yp(1)=y(2);
yp(2)=-y(1)+esp*y(2)*(1-y(1)^2);
end