function x = LUrisolvi(LU, x)

% Risolve il sistema lineare LU * x = b, con pivoting parziale.
%
% Input:
%   LU = matrice quadrata contenente i coefficienti del sistema
%   x  = vettore dei termini noti
%
% Output:
%   x  = soluzione del sistema lineare
%                                           2-11-2025
if nargin < 2 
    error('Errore: numero di argomenti insufficiente. Fornire LU e il vettore b.'); 
end

[rows, cols] = size(LU);
b_len = length(x);
if rows ~= cols || b_len ~= rows
    error('Errore: dimensioni non compatibili. Controllare matrice e vettore.');
end
[LU, perm] = LU_pivoting_parziale(LU);
x = x(:);
x = x(perm);
for i = 1:rows
    x(i+1:rows) = x(i+1:rows) - LU(i+1:rows, i) * x(i);
end

for i = rows:-1:1
    if LU(i,i) == 0
        error('Errore: fattore U singolare, impossibile proseguire.');
    end 
    x(i) = x(i) / LU(i,i);
    x(1:i-1) = x(1:i-1) - LU(1:i-1, i) * x(i);
end

end
