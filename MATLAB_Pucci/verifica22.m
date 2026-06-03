function verifica22()
a = -5;
b = 5;

f = @(x) 1./(1 + x.^2); 
df = @(x) -2*x ./ ((1 + x.^2).^2); 

num_test_points = 10001;
x_vals = linspace(a, b, num_test_points);
f_vals = f(x_vals); 

err_lag = zeros(1,20); 
err_new = zeros(1,20);   
idx = 1;

for n = 5:5:100
    nodes_count = n + 1;
    x_nodes = chebyshev(a, b, nodes_count);
    f_nodes = f(x_nodes);
  
    y_lag = lagrange(x_nodes, f_nodes, x_vals);
    err_lag(idx) = max(abs(f_vals - y_lag));

    y_new = newton(x_nodes, f_nodes, x_vals);
    err_new(idx) = max(abs(f_vals - y_new));

    idx = idx + 1;
end

err_her = zeros(1,10); 
err_her_der = zeros(1,10); 
k = 1;

for n = 5:5:50
    nodes_count = n + 1;
    x_nodes = chebyshev(a, b, nodes_count);
    f_nodes = f(x_nodes);
    df_nodes = df(x_nodes); 
    
    [y_her, dy_her] = interpolazione_hermite(x_nodes, f_nodes, df_nodes, x_vals);
    err_her(k) = max(abs(f_vals - y_her));

    df_vals = df(x_vals); 
    err_her_der(k) = max(abs(df_vals - dy_her));

    k = k + 1;
end

figure; 
semilogy(5:5:100, err_lag, 'o-', 'LineWidth', 2, 'DisplayName', 'Err Lagrange') 
hold on;
semilogy(5:5:100, err_new, 'o-', 'LineWidth', 2, 'DisplayName', 'Err Newton')   
hold off; 

xlim([0 100]);
xticks(0:10:100);
xlabel('n');
ylabel('Errore (log)'); 
title("Confronto Errori Lagrange e Newton"); 
legend('Location', 'best'); 
grid on;

figure;
semilogy(5:5:50, err_her, 'o-', 'LineWidth', 2, 'DisplayName', 'Err Hermite') 
hold on; 
semilogy(5:5:50, err_her_der, 'o-', 'LineWidth', 2, 'DisplayName', 'Err Der Hermite')
hold off; 

xlim([0 60]);
xlabel('n');
ylabel('Errore (log)');
title("Errori Hermite e Derivata"); 
legend('Location', 'best'); 
grid on;

end
