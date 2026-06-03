function grafico_spline_esercizio25()
    % Definisce il numero di punti di interpolazione (n) e il passo (h)
    vettore_n_punti = 100:100:1000;
    passo_h = 10 ./ vettore_n_punti; % Calcolo del passo

    % Funzione da interpolare: f(x) = 1 / (1 + x^2)
    funz_f = @(x) 1 ./ (1 + x.^2);

    % Pre-allocazione degli array per gli errori massimi
    errore_max_naturale = zeros(1, 10);
    errore_max_na_knot = zeros(1, 10);

    % Intervallo di interpolazione: [-5, 5]
    limite_a = -5;
    limite_b = 5;

    % Punti per la valutazione della funzione esatta
    punti_valutazione_x = linspace(limite_a, limite_b, 10001);
    valori_funzione_x = funz_f(punti_valutazione_x);

    % Ciclo per calcolare gli errori per diverse densità di punti (n)
    for i = 1:10
        % Genera i punti di interpolazione e i valori della funzione
        punti_interpolazione_xn = linspace(limite_a, limite_b, vettore_n_punti(i) + 1);
        valori_funzione_xn = funz_f(punti_interpolazione_xn);

        % Calcola l'errore massimo per la spline naturale
        errore_max_naturale(i) = norm(valori_funzione_x - spline3(punti_interpolazione_xn, valori_funzione_xn, punti_valutazione_x, 0), inf);
        % Calcola l'errore massimo per la spline not-a-knot
        errore_max_na_knot(i) = norm(valori_funzione_x - spline3(punti_interpolazione_xn, valori_funzione_xn, punti_valutazione_x, 1), inf);
    end

    % Linee di riferimento per l'ordine di convergenza
    riferimento_h2 = passo_h.^2; % Per spline naturale
    riferimento_h4 = passo_h.^4; % Per spline not-a-knot

    % --- Grafici per l'Intervallo [-5, 5] ---
    % Crea una nuova figura
    figure;
    % Definisce un layout a griglia 1x2 per i sottografici
    tiledlayout(1, 2);

    % Primo sottografico (spline not-a-knot)
    nexttile;
    loglog(passo_h, errore_max_na_knot, 'v--', 'Color', [0.85 0.33 0.10], 'DisplayName', 'Errore spline not-a-knot');
    hold on; % Mantiene il grafico attivo per aggiungere la prossima curva
    loglog(passo_h, riferimento_h4, '-.', 'Color', [0.75 0.75 0.00], 'DisplayName', '$h^4$', 'LineWidth', 1.2);
    hold off; % Rilascia il grafico
    title("Errore massimo spline not-a-knot [-5,5]");
    xlabel('Passo h');
    ylabel('Errore Massimo');
    legend('Location', 'southwest', 'Interpreter', 'latex');

    % Secondo sottografico (spline naturale)
    nexttile;
    loglog(passo_h, errore_max_naturale, 'v--', 'Color', [0.85 0.33 0.10], 'DisplayName', 'Errore spline naturale');
    hold on;
    loglog(passo_h, riferimento_h2, '-.', 'Color', [0.75 0.75 0.00], 'DisplayName', '$h^2$', 'LineWidth', 1.2);
    hold off;
    title("Errore massimo spline naturale [-5,5]");
    xlabel('Passo h');
    ylabel('Errore Massimo');
    legend('Location', 'southwest', 'Interpreter', 'latex');
end