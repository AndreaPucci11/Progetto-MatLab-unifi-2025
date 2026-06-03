function verificaLUEs9() 
    maxErr = zeros(20,1);
    dims = 1:20;

    for k = 1:20
        [M, rhs] = generaDatiES9(k);
        try
            rhs = LUrisolvi(M, rhs);
        catch
            error("Errore durante la risoluzione del sistema per n = %d.", k);
        end
        
        maxErr(k) = norm(abs(1 - rhs), inf);
    end

    loglog(dims, maxErr, '-o')
    xlabel('Dimensione matrice (n × n)'), ylabel('Errore massimo (norma inf)')
    title("Andamento dell'errore massimo al variare della dimensione n")
end

% Funzione di supporto che genera la matrice e il termine noto
% per il test 
function [M, rhs] = generaDatiES9(n)
    if n <= 0
        error("Il parametro n deve essere strettamente positivo.");
    end

    M = zeros(n);
    rhs = zeros(n,1);

    for col = 1:n
        M(:, col) = col;
    end

    for row = 1:n
        M(row, :) = M(row, :).^ (row - 1);  
        rhs(row) = sum(M(row, :));         
    end
end
