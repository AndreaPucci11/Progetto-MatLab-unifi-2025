function verificaLDLEs11()

    maxErr = zeros(12,1);
    dim = 1:12;

    for k = 1:12
        [M, v] = creaDatiES11(k);
        try
            v = LDL(M, v);
        catch
            error("Errore durante la risoluzione con LDLsolve per dimensione %d.", k);
        end
        maxErr(k) = norm(abs(1 - v), inf);
    end

    loglog(dim, maxErr, '-o')
    xlabel('n (dimensione della matrice)'), ylabel('Errore ∞-norm')
    title("Andamento dell'errore massimo al crescere di n");
end

% Funzione di supporto: genera la matrice M e il vettore v
% necessari per l'esercizio 11

function [M, v] = creaDatiES11(n)
    if n <= 0
        error("Il parametro n deve essere un intero positivo.");
    end

    M = zeros(n);
    v = zeros(n,1);

    for r = 1:n
        M(r, :) = (r : r+n-1).^(-1);
        v(r) = sum(M(r, :));
    end
end
