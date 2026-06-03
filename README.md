# Progetto-MatLab-unifi-2025
Progetto MatLab per laurea triennale in informatica svolto in gruppo.
# Elaborato di Calcolo Numerico - Università degli Studi di Firenze (UNIFI)

Questo repository contiene l'intera suite di funzioni MATLAB e la relativa relazione tecnica dedicate alla progettazione, implementazione e analisi di stabilità numerica per algoritmi di calcolo scientifico avanzato. Il sistema è stato interamente sviluppato da me in conformità con i requisiti d'esame per il corso di Calcolo Numerico.

L'analisi teorica completa, le tabelle di tolleranza e i grafici di convergenza sono consultabili all'interno della documentazione allegata: `progetto_calcolo_numerico_Pucci.pdf`.

## 👤 Autori (A.A. 2024/2025)
* **Andrea Pucci** - Matricola: 7119077
* **Antonio Mascani** - Matricola: 7118851
---

## 📝 Panoramica del Progetto
Il progetto affronta la risoluzione di problemi matematici e ingegneristici complessi mediante l'ausilio dell'analisi numerica computerizzata, focalizzandosi su cinque aree macro-tematiche distribuite in 30 esercizi guidati:

1. **Teoria degli Errori e Approssimazioni** (Es. 1-2)
2. **Zeri di Funzioni e Sistemi Non Lineari** (Es. 3-6, 15-16)
3. **Algebra Lineare Numerica e Minimi Quadrati** (Es. 7-14, 23)
4. **Approssimazione e Interpolazione Polinomiale / Spline** (Es. 17-22, 24-27)
5. **Quadrature e Integrazione Numerica** (Es. 28-30)

---

## 📐 Struttura dell'Elaborato e Algoritmi Sviluppati

### 1. Teoria degli Errori e Rappresentazione del Calcolatore
* **Approssimazione Numerica di Derivate (Es. 1):** Combinazione lineare di nodi mediante espansione in serie di Taylor per ricavare un'approssimazione della derivata prima $y'(t)$ con un errore di troncamento di ordine $O(h^4)$.
* **Arritmetica in Fluctuating Point (Es. 2):** Studio approfondito dello standard IEEE 754 in quadrupla precisione (128 bit), calcolo della precisione di macchina ($u = 2^{-113}$), limiti di normalizzazione e dell'errore assoluto di rappresentazione per frazioni periodiche ($1/3$).

### 2. Risoluzione di Equazioni e Sistemi Non Lineari
* **Metodo di Newton Standard (Es. 3):** Funzione `newtonStandard` per il calcolo delle radici di una funzione con gestione dinamica dei parametri di tolleranza e iterazioni massime.
* **Metodo di Newton Modificato (Es. 4):** Funzione `newtonModificato` ottimizzata per preservare la convergenza quadratica anche in presenza di radici con molteplicità $m$.
* **Metodo di Bisezione (Es. 5):** Funzione `Bisezione` operante su intervalli con cambio di segno.
* **Analisi di Costo e Convergenza (Es. 6):** Studio e tabulazione comparativa del numero di valutazioni funzionali necessarie per risolvere $f(x) = x - \cos(x)$ e la sua variante con radice tripla $g(x) = f(x)^3$, analizzando la perdita di convergenza quadratica del metodo standard.
* **Sistemi Non Lineari Multivariabili (Es. 15-16):** Funzione `newton_non_lineare` per la ricerca di punti stazionari di funzioni a più variabili mediante l'utilizzo del vettore gradiente e della matrice Jacobiana valutata a ogni iterazione.

### 3. Sistemi Lineari e Problemi di Minimi Quadrati
* **Sistemi Triangolari Compatti (Es. 7):** Risoluzione efficiente tramite la funzione `trilows` di sistemi in cui la matrice dei coefficienti $L$ è memorizzata in forma vettoriale compressa di lunghezza $n(n+1)/2$.
* **Fattorizzazione LU con Pivoting (Es. 8-9):** Algoritmo `LU_pivoting_parziale` per matrici quadrate non singolari e validazione del comportamento asintotico distruttivo dell'errore causato dal malcondizionamento delle matrici di Vandermonde (testato fino a $n=20$).
* **Fattorizzazione $LDL^\top$ (Es. 10-11):** Algoritmo per matrici simmetriche definite positive ed evidenza numerica del limite di stabilità finita sulle matrici di Hilbert (risolvibili solo fino a $n=12$).
* **Sistemi Sovradeterminati e QR (Es. 12):** Risoluzione nel senso dei minimi quadrati lineari tramite decomposizione QR implicita basata sulle trasformazioni di Householder, con calcolo della norma del residuo quadratico.
* **Minimi Quadrati Pesati (Es. 13-14):** Algoritmo `QR_minimi_quadrati_pesata` per la minimizzazione della norma pesata del residuo $\|r\|_\omega^2$ con vettori di peso assegnati.
* **Sistemi Tridiagonali (Es. 23):** Risoluzione a complessità lineare mediante l'algoritmo ottimizzato `tridia`.

### 4. Interpolazione ed Approssimazione di Funzioni
* **Forme di Lagrange e Newton (Es. 17-18):** Funzioni `lagrange` e `newton` per l'interpolazione polinomiale classica (differenze divise).
* **Interpolazione di Hermite (Es. 19):** Calcolo del polinomio interpolante e della sua derivata prima basato su punti e derivate assegnati.
* **Nodi di Chebyshev (Es. 20):** Algoritmo per la generazione delle ascisse ottimali di Chebyshev atte a minimizzare l'errore teorico di interpolazione su intervalli arbitrari.
* **Analisi del Fenomeno di Runge (Es. 21-22):** Studio comparativo dell'errore massimo sulla funzione di Runge $f(x) = (1+x^2)^{-1}$. Evidenza sperimentale del malcondizionamento esponenziale $O(2^n)$ indotto dai nodi equispatiati e dimostrazione della stabilità della forma di Lagrange rispetto a Newton/Hermite in presenza di nodi di Chebyshev in aritmetica finita.
* **Spline Cubiche (Es. 24-26):** Funzione `spline3` per la costruzione di spline cubiche interpolanti con condizioni al bordo di tipo *Naturale* o *Not-a-knot*. Verifica sperimentale dell'ordine di errore asintotico su intervalli multipli ($h^2$ per spline naturali vs $h^4$ per spline not-a-knot).
* **Minimi Quadrati Polinomiali (Es. 27):** Sperimentazione su dati perturbati mediante l'utilizzo combinato di `polyfit` per valutare l'andamento dell'errore al variare del grado polinomiale $n$.

### 5. Integrazione Numerica (Quadratura)
* **Coefficienti di Newton-Cotes (Es. 28):** Sviluppo della funzione `ncotes` per il calcolo esatto in formato razionale dei coefficienti di quadratura per formule di grado $n$.
* **Formule Composite (Es. 29-30):** Implementazione di `ncotes_composita` per l'integrazione numerica su sottointervalli. Analisi dell'efficienza e tabulazione comparativa degli errori residui tra le formule composite dei Trapezi e di Simpson applicate all'approssimazione dell'integrale definito $\int_{0}^{10}\frac{1}{1+x^2}dx \equiv \arctan(10)$.

---

## 🛠️ Requisiti d'Esecuzione
* **Software:** MATLAB (Versione R2020a o successive consigliata) o GNU Octave compatibile.
* **Modalità d'uso:** Ciascun esercizio è mappato su una corrispondente funzione `.m` autonoma. All'interno della relazione in PDF sono riportati gli esempi di inizializzazione e i comandi di test esatti per validare l'output di ogni singolo algoritmo direttamente dalla CLI di MATLAB.
