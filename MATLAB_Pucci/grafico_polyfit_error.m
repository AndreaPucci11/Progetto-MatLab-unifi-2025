function grafico_polyfit_error()
    errore_max = zeros(1, 40); 

    for n = 1:40
        err = calcola_err(n);                
        errore_max(n) = norm((err),inf);
    end

    semilogy(1:40, errore_max, 'o-','Color', [0.8 0.0 0.0], 'LineWidth', 2 )
    xlabel('Grado del polinomio n')
    ylabel('Massimo errore |p(x)-y(x)|')
    title('Massimo errore in funzione del grado del polinomio')

end


function err = calcola_err(n)
x = linspace(-5,5,1001);
y = 1./(1+x.^2);
rng(0)
r = ( rand(1,1001)-rand(1,1001) )/5;
yr = y + r;
%plot(x,yr,'r.',x,y,'b','LineWidth',2)
%xlabel('x'), ylabel('y')
%shg

p = polyfit(x,yr,n);
p = polyval(p,x);
err = p-y;
end