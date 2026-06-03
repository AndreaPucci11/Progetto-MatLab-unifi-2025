function valori = newton(punti_x, valori_f, x_valori)
  
% Newton polinomio interpolante
%
% Input:    
% punti_x = vettore delle ascisse di interpolazione
% valori_f = valori della funzione corrispondenti alle ascisse
% x_valori = vettore dei punti in cui valutare il polinomio
%    
% Output:
% valori = vettore dei valori del polinomio interpolante
%                                           2-11-2025
    n = length(punti_x);
    if length(valori_f) ~= n
        error('I vettori punti_x e valori_f devono avere la stessa lunghezza.');
    end

    % Costruzione della tabella delle differenze divise
    tabella_divise = zeros(n, n);
    tabella_divise(:, 1) = valori_f(:);
    for col = 2:n
        for riga = 1:n - col + 1
            tabella_divise(riga, col) = ...
                (tabella_divise(riga+1, col-1) - tabella_divise(riga, col-1)) / ...
                (punti_x(riga+col-1) - punti_x(riga));
        end
    end

    % Coefficienti del polinomio di Newton
    coeff_newton = tabella_divise(1, :);

    % Valutazione del polinomio nei punti richiesti
    valori = zeros(size(x_valori));
    for k = 1:length(x_valori)
        x_corrente = x_valori(k);
        termine = coeff_newton(1);
        prodotto = 1;
        for j = 2:n
            prodotto = prodotto * (x_corrente - punti_x(j-1));
            termine = termine + coeff_newton(j) * prodotto;
        end
        valori(k) = termine;
    end
end
