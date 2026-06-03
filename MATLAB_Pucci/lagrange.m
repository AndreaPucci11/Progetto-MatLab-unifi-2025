function  y=lagrange(xi,fi,x)

%
%  y=lagrange(xi,fi,x)
% 
% Input: 
% xi = ascisse di interpolazione
% fi = valore della funzione nelle ascisse di interpolazione
% x = vettore di punti
%
% Output:
% y = vettore contenente le valutazioni del polinomio interpolante
%                                                Rel. 2-11-2025

if nargin < 3
    error('Dati insufficienti');
end


grado_n = length(xi);
if grado_n ~= length(fi)
    error('dati␣inseriti␣errati');
end
y = zeros(size(x));  

    for i = 1:grado_n
        L = ones(size(x));  
        for j = 1:grado_n
            if j ~= i
                L = L .* (x - xi(j)) / (xi(i) - xi(j)); 
            end
        end
        y = y + fi(i) * L; 
    end
end
