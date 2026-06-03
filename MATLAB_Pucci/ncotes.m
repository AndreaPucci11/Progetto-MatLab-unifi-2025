function coefficiente_interpolazione = ncotes(grado)

% coefficiente_interpolazione = ncotes(grado)
% Calcola i coefficienti della formula di N-C di grado 'grado'
%                       Rel. 2-11-2025

    if nargin < 1, error('input errato'), end

    coefficiente_interpolazione = zeros(grado + 1, 1);

    for indice = 0 : ceil(grado / 2)
        coefficiente_interpolazione(indice + 1) = calcola_coefficiente_i(indice, grado);
        coefficiente_interpolazione(grado + 1 - indice) = coefficiente_interpolazione(indice + 1);
    end

    return
end

function coeff_i = calcola_coefficiente_i(i, grado)

    vettore_nn = [0 : i - 1, i + 1 : grado];    
    fattore_scala = prod(i - vettore_nn);
    polinomio_a = poly(vettore_nn);
    denominatori_integrazione = (grado + 1) : -1 : 1;
    polinomio_a_integrato = [polinomio_a ./ denominatori_integrazione, 0];
    coeff_i = polyval(polinomio_a_integrato, grado) / fattore_scala;

    return
end