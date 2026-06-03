% Esegue la fattorizzazione QR di una matrice A utilizzando
% le trasformazioni di Householder (metodo implicito).
%
% Input:
%   A = matrice di dimensione m x n, con m > n e rango(A) = n
%
% Output:
%   A = matrice contenente la fattorizzazione QR
%                                           2-11-2025

function A = qr_fattorizza(A)
    [m, n] = size(A);
    if m < n
        error('La matrice deve avere più righe che colonne (sistema sovradeterminato).');
    end
    
    for i = 1:n
        alfa = norm(A(i:m, i));
        if alfa == 0
            error('La matrice non è a rango pieno: impossibile proseguire la fattorizzazione.');
        end
        if A(i, i) > 0
            alfa = -alfa;
        end
        v1 = A(i, i) - alfa;
        A(i, i) = alfa;
        A(i+1:m, i) = A(i+1:m, i) / v1;
        beta = -v1 / alfa;
        A(i:m, i+1:n) = A(i:m, i+1:n) - (beta * [1; A(i+1:m, i)]) * [1, A(i+1:m, i)'] * A(i:m, i+1:n);
    end
    
    return;
end
