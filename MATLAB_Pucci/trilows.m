function x = trilows(L,b)
% x = trilow(L,b)
%
% Input:
% L = la matrice dei coefficienti L è passata 
% come un vettore di lunghezza n(n + 1)/2.
%
% b= vettore dei termini noti
%
% Output: soluzione del sistema
%                                                Rel. 2-11-2025

lung= length(b);
if lung*(lung+1)/2 ~= length(L)
    error ("Dati errati");
end

m=0;
L=L(:);
b=b(:);

for i=1:lung
    if L(i+m)==0
        error ("Matrice singolare");
    end
    b(i) = b(i)- L(m+1:m+i-1)'*b(1:i-1);
    b(i) = b(i) / L(m+i);
    m=m+i;
end

x=b;
return
end




