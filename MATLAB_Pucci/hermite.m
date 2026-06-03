function [y, y1] = hermite(xi, fi, fi1, x)
%   [y, y1] = hermite(xi, fi, fi1, x)
%
%   Input:
%     xi  : Ascisse di interpolazione.
%     fi  : Valori della funzione nelle ascisse di interpolazione (f(xi)).
%     fi1 : Valori della derivata della funzione nelle ascisse di interpolazione (f'(xi)).
%     x   : Vettore di punti in cui valutare il polinomio.
%
%   Output:
%     y   : Valutazioni del polinomio di Hermite nei punti x.
%     y1  : Valutazioni del polinomio di Hermite derivato nei punti x.
%
%                                                 Rel. 2-11-2025

% Controllo dei dati messi in input
if nargin < 4
    error('hermite:datiInsufficienti', 'Errore: Dati insufficienti. Servono almeno 4 argomenti.');
end

n = length(xi) - 1;
if n ~= length(fi) - 1
    error('hermite:datiErrati', 'Errore: Le lunghezze di xi e fi errate.');
end


f_estesa = zeros(1, 2 * length(fi));
f_estesa(1:2:end) = fi;   
f_estesa(2:2:end) = fi1;   

% Duplica le ascisse di interpolazione per le differenze divise
xi_estesa = zeros(1, 2 * length(xi));
xi_estesa(1:2:end) = xi;
xi_estesa(2:2:end) = xi;

% Calcola le differenze divise di Hermite
differenza_divise = calcoloHermiteDifferenzaDivise(xi_estesa, f_estesa);


y = differenza_divise(2*n+2) * ones(size(x));
y1 = zeros(size(x));


for i = (2*n+1):-1:1
    y1 = y + (x - xi_estesa(i)) .* y1;
    y = y .* (x - xi_estesa(i)) + differenza_divise(i);
end

end 

% --- Sottofunzione per il calcolo delle differenze divise di Hermite ---
function fi_dd = calcoloHermiteDifferenzaDivise(xi_estese, fi_estese)

%   fi_dd = calcoloHermiteDifferenzeDivise(xi_ext, fi_ext)
%
%   Input:
%     xi_estese : Ascisse estese per le differenze divise.
%     fi_estese : Valori della funzione e della derivata estesi.
%
%   Output:
%     fi_dd  : Vettore delle differenze divise di Hermite.
%
%                                              2-11-2025 

if nargin < 2
    error('dati insufficienti per il calcolo delle differenze divise');
end

numero_punti = length(xi_estese);
fi_dd = fi_estese;


for i = (numero_punti - 1):-2:3
    fi_dd(i) = (fi_dd(i) - fi_dd(i-2)) / (xi_estese(i) - xi_estese(i-1));
end


for j = 2:(numero_punti - 1)
    for i = numero_punti:-1:(j + 1)
        fi_dd(i) = (fi_dd(i) - fi_dd(i-1)) / (xi_estese(i) - xi_estese(i-j));
    end
end

end