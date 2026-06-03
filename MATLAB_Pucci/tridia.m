function  x = tridia(b,a,c,y)

%
% x = tridia(b,a,c,y)
% Input:
% a = elementi delle diagonali
% b = elementi delle sottodiagonali
% c = elementi delle sopradiagonali
% y = termini noti
% 
% Output:
% x = soluzioni del sistema tridiagonale Ax=y
% 
%                                           Rel. 2-11-2025


n=length(a);
if n ~= length(b)+1 ||...
    n ~= length(c)+1 ||...
    n~= length(y) 
    error('dati errati');
end

if a(1)==0
    error('fattorizzazione indefinita');
end

for i = 1:n-1
    b(i) = b(i)/a(i);
    a(i+1) = a(i+1) - b(i)*c(i);
    if a(i+1)==0
        error('fattorizzazione indefinita');
    end
end

x=y(:);
for i=2:n
    x(i)=x(i)-b(i-1)*x(i-1);
end

x(n) = x(n) / a(n);
for i = n-1:-1:1
    x(i) = (x(i) - c(i) * x(i+1)) / a(i);
end

end