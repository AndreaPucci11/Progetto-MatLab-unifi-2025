function val_spline = spline3(ascisse_interpolazione, valori_funzionali, punti_valutazione, tipo)
% val_spline = spline3(ascisse_interpolazione, valori_funzionali, punti_valutazione, tipo)
% Input:
% ascisse_interpolazione : n+1 ascisse di interpolazione
% valori_funzionali : n+1 valutazioni funzionali nelle ascisse_interpolazione
% punti_valutazione : vettore di punti in cui valutare la spline
% tipo : tipo di spline; 0 naturale, 1 not-a-knot
% Output:
% val_spline : vettore di valori della spline nei punti in punti_valutazione
%
%
%                                                     Rel. 2-11-2025

    if nargin < 3, error('argomenti insufficienti'), end
    if nargin == 3, tipo = 1; end

    n_punti = length(ascisse_interpolazione) - 1;
    if n_punti ~= length(valori_funzionali) - 1 || (tipo ~= 0 && tipo ~= 1), error('dati incorretti'), end

    valori_funzionali = valori_funzionali(:);
    ascisse_interpolazione = ascisse_interpolazione(:);

    delta_h = ascisse_interpolazione(2:n_punti+1) - ascisse_interpolazione(1:n_punti);
    phi = delta_h(1:n_punti-1) ./ (delta_h(1:n_punti-1) + delta_h(2:n_punti));
    Xi = delta_h(2:n_punti) ./ (delta_h(1:n_punti-1) + delta_h(2:n_punti));
    diff_prime = (valori_funzionali(2:n_punti+1) - valori_funzionali(1:n_punti)) ./ delta_h;
    diff_seconde = (diff_prime(2:n_punti) - diff_prime(1:n_punti-1)) ./ (ascisse_interpolazione(3:n_punti+1) - ascisse_interpolazione(1:n_punti-1));

    m = zeros(n_punti + 1, 1);

    % m trovato con tridia
    if tipo == 0
        m(2:n_punti) = tridia(phi(2:n_punti-1),2*ones((n_punti-1),1),Xi(1:n_punti-2),6*diff_seconde);
    else
        % Codice per spline not-a-knot
        if length(m) <= 3, error('ascisse insufficienti'); end
        m(2:n_punti) = tridia([phi(2:n_punti-2); phi(n_punti-1)-Xi(n_punti-1)], ...
                              [2-phi(1); 2*ones(n_punti-3,1); 2-Xi(n_punti-1)], ...
                              [Xi(1)- phi(1); Xi(2:n_punti-2)], ...
                              6*[diff_seconde(1)*Xi(1); diff_seconde(2:n_punti-2); diff_seconde(n_punti-1)*phi(n_punti-1)]);
        m(1) = 6*diff_seconde(1) - m(2) - m(3);
        m(n_punti+1) = 6*diff_seconde(n_punti-1) - m(n_punti) - m(n_punti-1);
    end

    % calcolo q e r
    q = diff_prime + (m(1:n_punti) - m(2:n_punti+1)) .* delta_h / 6;
    r = valori_funzionali(1:n_punti) - (m(1:n_punti) .* (delta_h).^2) / 6;

    % preallocazione per il vettore di output
    val_spline = zeros(size(punti_valutazione));

    % calcolo spline in punti_valutazione
    for j = 1:length(punti_valutazione)
        indice_intervallo = trova_intervallo(punti_valutazione(j), ascisse_interpolazione);

        % Valutazione della spline nel punto j
        val_spline(j) = (m(indice_intervallo+1) * (punti_valutazione(j) - ascisse_interpolazione(indice_intervallo))^3 ...
                        + m(indice_intervallo) * (ascisse_interpolazione(indice_intervallo+1) - punti_valutazione(j))^3) / (6 * delta_h(indice_intervallo)) ...
                        + q(indice_intervallo) * (punti_valutazione(j) - ascisse_interpolazione(indice_intervallo)) ...
                        + r(indice_intervallo);
    end

    return
end

function indice = trova_intervallo(valore_x, ascisse_interp)
% indice = trova_intervallo(valore_x, ascisse_interp)
% Input:
% valore_x : ascissa
% ascisse_interp : ascisse di interpolazione
% Output:
% indice : intervallo in cui si trova valore_x, ritorna -1 se non trovato

    for indice = 1:length(ascisse_interp) - 1
        if valore_x >= ascisse_interp(indice) && valore_x <= ascisse_interp(indice+1)
            return
        end
    end
    indice = -1; % Se non trova l'intervallo, imposta a -1
    error("x non fuori dall'intervallo"); % Lancia un errore se x non è nell'intervallo
end