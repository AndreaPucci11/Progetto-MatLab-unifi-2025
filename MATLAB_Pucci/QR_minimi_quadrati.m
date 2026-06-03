function [x_soluzione, norma_residuo_quadrato] = QR_minimi_quadrati(matrice_A, vettore_b)

% Risolve un sistema sovradeterminato Ax = b nel senso dei minimi quadrati
% utilizzando la fattorizzazione QR tramite trasformazioni di Householder.
%
% Input:
%   matrice_A: matrice dei coefficienti m x n (m > n e rango = n)
%   vettore_b: vettore dei termini noti di dimensione m
%
% Output:
%   x_soluzione: vettore soluzione di dimensione n
%   norma_residuo_quadrato: quadrato della norma 2 del resto (||Ax - b||^2)
%                                           2-11-2025

    A_fattorizzata = qr_fattorizza(matrice_A);
    [m_righe, n_colonne] = size(A_fattorizzata);
    
    if length(vettore_b) ~= m_righe
        error('Il vettore b deve avere lo stesso numero di elementi delle righe di A.');
    end
    
    b_corrente = vettore_b(:);
    for i = 1:n_colonne
        vettore_householder = [1; A_fattorizzata(i+1:m_righe, i)];
        coeff_beta = 2 / (vettore_householder' * vettore_householder);
        b_corrente(i:m_righe) = b_corrente(i:m_righe) - coeff_beta * vettore_householder * (vettore_householder' * b_corrente(i:m_righe));
    end
    
    norma_residuo_quadrato = norm(b_corrente(n_colonne+1:m_righe))^2;
    
    for i = n_colonne:-1:1
        b_corrente(i) = b_corrente(i) / A_fattorizzata(i, i);
        b_corrente(1:i-1) = b_corrente(1:i-1) - A_fattorizzata(1:i-1, i) * b_corrente(i);
    end
    x_soluzione = b_corrente(1:n_colonne);

end
