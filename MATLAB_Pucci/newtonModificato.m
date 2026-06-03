% In input la funzione newtonModificato riceve:
% f = funzione passata
% x0 = punto di partenza scelto in modo arbitrario
% m = molteplicità funzione di partenza
% df = derivata funzione di partenza
% massima_iterazione = iterazioni massime che possono essere fatte
% tol = tolleranza scelta
%
% La funzione restituisce:
% x = approssimazione del risultato della funzione
% iter = numero di iterazioni impiegate 
%                                           2-11-2025
function [x, iterazione] = newtonModificato(f, x0, m, df, massima_iterazione, tol)
    if nargin < 5
        massima_iterazione = 100; 
    end
    if nargin < 6
        tol = 1e-11;
    end

    x = x0;
    controlloTol = 0;

    for iterazione = 1:massima_iterazione
        fx = f(x);
        der = df(x);

        if der == 0 && fx == 0
            break;        
        end
        if der == 0
            error('Non avviene convergenza (derivata nulla)');
        end

        nuovox = x - m * (fx / der);
        controlloTol = abs(nuovox - x);

        if controlloTol <= tol
            x = nuovox;
            return;
        end
        x = nuovox;
    end

    if controlloTol > tol
        warning("Non è stata soddisfatta la tolleranza");
    end
end
