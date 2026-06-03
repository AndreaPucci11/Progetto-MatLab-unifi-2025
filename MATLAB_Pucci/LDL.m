% Questa funzione utilizza la fattorizzazione LDL per risolvere il sistema
% Ax = b.
% Outuput: 
%   s: vettore soluzione del sistema lineare
% Input:
%   A_matrice: matrice dei coefficienti
%   s_rhs_vettore: Vettore dei termini noti
 %                                          2-11_2025                                       
function s = LDL(A_matrice, s_rhs_vettore)

    if nargin < 2
        error('gli argomenti da inserire sono 2');
    end
    [n, m] = size(A_matrice);

   
    if n ~= m || m ~= length(s_rhs_vettore)
        error(['Dimensioni non valide: la matrice dev\''essere quadrata e' ...
            ' la lunghezza delle righe e delle colonne deve corrispondere alla lunghezza di b']);
    end

    A_fattorizzata = ldl_fattorizzazione(A_matrice);
        vettore_soluzione = s_rhs_vettore(:);
    for i = 1:n
        vettore_soluzione(i+1:n) = vettore_soluzione(i+1:n)- A_fattorizzata(i+1:n, i) * vettore_soluzione(i);
    end
    vettore_soluzione = vettore_soluzione ./ diag(A_fattorizzata);

    for i = n:-1:1
        vettore_soluzione(1:i-1) = vettore_soluzione(1:i-1) - A_fattorizzata(i, 1:i-1)' * vettore_soluzione(i);
    end
    s = vettore_soluzione;
end

function A = ldl_fattorizzazione(A)
    [m,n] = size(A);
    if m ~= n
        error('matrice non quadrata');
    end
    if A(1,1) <= 0
        error('A non sdp');
    end
    
    A(2:n,1) = A(2:n,1) / A(1,1);
     
    for j = 2:n
        v = A(j,1:j-1)' .* diag(A(1:j-1,1:j-1));
        A(j,j) = A(j,j) - A(j,1:j-1) * v;
        if A(j,j) <= 0 
            error('A non semidefinita positiva');
        end
        A(j+1:n,j) = (A(j+1:n,j) - A(j+1:n, 1:j-1)*v) / A(j,j);
    end
    return
end