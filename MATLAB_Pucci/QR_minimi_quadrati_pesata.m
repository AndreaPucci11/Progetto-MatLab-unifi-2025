function [x_soluzione, norma_residuo_quadrato_pesata] = QR_minimi_quadrati_pesata(matrice_A, vettore_b, vettore_pesi)

% Risolve il problema dei minimi quadrati pesati (WLS) per il sistema sovradeterminato A*x = b.
% L'obiettivo è minimizzare il quadrato della norma euclidea pesata: || D * (A*x - b) ||_2^2, 
% dove D è la matrice diagonale dei fattori di peso (D = diag(sqrt(vettore_pesi))).
% La soluzione è calcolata tramite la fattorizzazione QR (metodo di Householder implicito).
%
% Input:
%   matrice_A: Matrice dei coefficienti m x n (m > n e rank(A) = n).
%   vettore_b: Vettore dei termini noti di dimensione m.
%   vettore_pesi: Vettore dei pesi (w_i) di dimensione m.
% Output:
%   x_soluzione: Vettore soluzione x di dimensione n.
%   norma_residuo_quadrato_pesata: Valore minimo di || D * (A*x - b) ||_2^2.
   %                                           2-11-2025
[num_righe, num_colonne] = size(matrice_A);
if length(vettore_b) ~= num_righe || length(vettore_pesi) ~= num_righe
    error('Le dimensioni del vettore b o del vettore dei pesi non corrispondono al numero di righe (m) della matrice A.');
end

vettore_b_trasformato = vettore_b(:);
fattori_peso_radice = sqrt(vettore_pesi(:));
vettore_b_trasformato = fattori_peso_radice .* vettore_b_trasformato;
matrice_A_trasformata = fattori_peso_radice .* matrice_A;

matrice_R_e_riflessioni = qr_fattorizza(matrice_A_trasformata);

for i = 1:num_colonne
    vettore_householder = [1; matrice_R_e_riflessioni(i+1:num_righe, i)];
    fattore_beta = 2 / (vettore_householder' * vettore_householder);
    vettore_b_trasformato(i:num_righe) = vettore_b_trasformato(i:num_righe) - fattore_beta * vettore_householder * (vettore_householder' * vettore_b_trasformato(i:num_righe));
end

parte_residua_b_tilde = vettore_b_trasformato(num_colonne+1:num_righe);
norma_residuo_quadrato_pesata = norm(parte_residua_b_tilde)^2;

for i = num_colonne:-1:1
    vettore_b_trasformato(i) = vettore_b_trasformato(i) / matrice_R_e_riflessioni(i, i);
    vettore_b_trasformato(1:i-1) = vettore_b_trasformato(1:i-1) - matrice_R_e_riflessioni(1:i-1, i) * vettore_b_trasformato(i);
end
x_soluzione = vettore_b_trasformato(1:num_colonne);

end