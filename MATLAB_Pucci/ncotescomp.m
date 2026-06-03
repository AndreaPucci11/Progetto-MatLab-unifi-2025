function integrale_approssimato = ncotescomp(funzione_f, limite_a, limite_b, grado_k, numero_intervalli)
% integrale_approssimato = ncotes_composita(funzione_f, limite_a, limite_b, grado_k, numero_intervalli)
% Input:
%    funzione_f        : funzione da integrare (handle della funzione)
%    limite_a, limite_b : estremi di integrazione [a,b]
%    grado_k           : grado della formula di Newton-Cotes (k)
%    numero_intervalli : numero totale di sottointervalli (n)
% Output:
%    integrale_approssimato : approssimazione dell'integrale di funzione_f in [a,b]
%                                   Rel. 2-11-2025

    
    if nargin < 5, error('dati insufficienti'), end
   
    if limite_b <= limite_a || mod(numero_intervalli, grado_k) ~= 0 || numero_intervalli < 1
        error('dati non validi: controllare estremi, n o k');
    end
    
    if grado_k < 1 || grado_k == 8 || grado_k > 9 
        error('problema malcondizionato: k deve essere tra 1 e 9, escluso 8');
    end

    passo_h = (limite_b - limite_a) / numero_intervalli;
    
    
    punti_x = (limite_a : passo_h : limite_b)';
    valori_f = funzione_f(punti_x);
    coefficienti_nc = ncotes(grado_k); 
    integrale_approssimato = 0;
    
    for indice_i = 1 : numero_intervalli / grado_k
        indici_correnti = ((indice_i - 1) * grado_k : indice_i * grado_k) + 1;
        valori_sottointervallo = valori_f(indici_correnti);
           
        integrale_approssimato = integrale_approssimato + sum(coefficienti_nc .* valori_sottointervallo);
    end
        
    integrale_approssimato = integrale_approssimato * passo_h;
    
    return 
end