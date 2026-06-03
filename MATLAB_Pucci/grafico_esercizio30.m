function grafico_esercizio30
n_punti = 10:10:50;
errore_trapezi = zeros(1, 5);
errore_simpson = zeros(1, 5);
funzione_da_integrare = @(x) 1./(1+x.^2);
valore_reale_arctan10 = atan(10);

for i=1:5
    integrale_trapezi = ncotescomp(funzione_da_integrare,0,10,1,n_punti(i));
    errore_trapezi(i) = abs(valore_reale_arctan10 - integrale_trapezi);
    integrale_simpson = ncotescomp(funzione_da_integrare,0,10,2,n_punti(i));
    errore_simpson(i) = abs(valore_reale_arctan10 - integrale_simpson);
end
figure;

subplot(1,2,2); 
loglog(n_punti, errore_simpson, 'o-', 'Color', [0.7 0.0 0.0] ,'LineWidth', 2)
xlabel('n')
ylabel('errore')
title('Errore N-C Simpson')
grid on

subplot(1,2,1);
loglog(n_punti, errore_trapezi, 'o-','Color', [0.7 0.0 0.0], 'LineWidth', 2)
xlabel('n')
ylabel('errore')
title('Errore N-C trapezi')
grid on