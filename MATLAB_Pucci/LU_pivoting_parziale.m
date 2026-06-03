function [A, perm] = LU_pivoting_parziale(A)

    % Esegue la fattorizzazione LU di A utilizzando pivoting
    % parziale (scambio di righe).
    %
    % Input:
    %   A = matrice quadrata e non singolare
    %
    % Output:
    %   A    = matrice contenente L e U sovrapposte
    %   perm = vettore di permutazione delle righe
    %                                           2-11-2025
    
    [rows, cols] = size(A);
    
    if rows ~= cols
        error('Errore: la matrice deve essere quadrata.');
    end
    
    perm = (1:rows)';
    
    for k = 1:rows
        [pivotVal, pivotIndex] = max(abs(A(k:rows, k)));
        
        if pivotVal == 0
            error('Errore: la matrice risulta singolare.');
        end
        
        pivotIndex = pivotIndex + k - 1;
        if pivotIndex ~= k
            A([pivotIndex, k], :) = A([k, pivotIndex], :);
            perm([pivotIndex, k]) = perm([k, pivotIndex]);
        end
        
        A(k+1:rows, k) = A(k+1:rows, k) / A(k, k);
        A(k+1:rows, k+1:cols) = A(k+1:rows, k+1:cols) - A(k+1:rows, k) * A(k, k+1:cols);
    end

end
