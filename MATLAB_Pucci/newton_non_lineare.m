function radice_approssimata = newton_non_lineare(funzione_handle, jacobiano_handle, vettore_innesco, tolleranza_relativa, numero_max_iterazioni)
% radice_approssimata = newton_non_lineare_italiano(funzione_handle, jacobiano_handle, vettore_innesco, tolleranza_relativa, numero_max_iterazioni)
% Input:
% funzione_handle : handle a un vettore di n funzioni non lineari F
% jacobiano_handle : handle alla matrice nxn Jacobiana di F
% vettore_innesco : vettore iniziale d'innesco x0
% tolleranza_relativa : tolleranza relativa per il criterio d'arresto
% numero_max_iterazioni : numero massimo di iterazioni consentite
% Output:
%
% radice_approssimata: approssimazione della radice
%                                                          Rel. 2-11-2025

    if nargin < 3
        error('numero argomenti di ingresso errato');
    end

    % I nomi delle variabili interne per le dimensioni rimangono spesso in inglese
    % per consuetudine (rows, cols), ma le rinominiamo per coerenza.
    [righe_jacobiano, colonne_jacobiano] = size(jacobiano_handle); 


    if righe_jacobiano ~= colonne_jacobiano || colonne_jacobiano ~= length(funzione_handle)
        error('dati in input errati: la Jacobiana deve essere quadrata e le sue dimensioni devono corrispondere al numero di funzioni/elementi del vettore d''innesco');
    end


    if nargin < 5
        numero_max_iterazioni = 1000;
    elseif numero_max_iterazioni < 1
        error('numero_max_iterazioni errato');
    end

    if nargin < 4
        tolleranza_relativa = 1e-12;
    elseif tolleranza_relativa < 0
        error('tolleranza_relativa negativa');
    end

    vettore_soluzione_corrente = vettore_innesco(:); % Assicura che sia un vettore colonna

    for contatore_iterazioni = 1:numero_max_iterazioni
        jacobiano_valutata = feval(jacobiano_handle, vettore_soluzione_corrente);

        funzione_valutata = feval(funzione_handle, vettore_soluzione_corrente);
        
        
        delta_soluzione = SoluzioniLU(jacobiano_valutata, -1 * funzione_valutata); % SoluzioniLU non viene rinominata in quanto è una funzione esterna
        
        vettore_soluzione_corrente = vettore_soluzione_corrente + delta_soluzione;
        
        errore_approssimazione_corrente = norm( abs(delta_soluzione) ./ (1 + abs(vettore_soluzione_corrente)));
        
        if errore_approssimazione_corrente <= tolleranza_relativa
            break;
        end

    end

    if errore_approssimazione_corrente > tolleranza_relativa
        error('Non converge entro il numero massimo di iterazioni');
    end
    
    radice_approssimata = vettore_soluzione_corrente; % Assegna il risultato alla variabile di output
    return
end