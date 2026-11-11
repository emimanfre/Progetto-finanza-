# Esercitazione 1 — GJR-GARCH(1,1) (Glosten, Jagannathan e Runkle, 1993)

> **Bozza** del piano di lavoro. Da rivedere e completare man mano.

## Il modello

$$y_t = \sigma_t z_t, \qquad z_t \sim \text{n.i.d.}(0,1), \qquad y_t \mid I_{t-1} \sim N(0, \sigma_t^2)$$

$$\sigma_t^2 = \omega + \alpha\, y_{t-1}^2 + \gamma\, d_{t-1} y_{t-1}^2 + \beta\, \sigma_{t-1}^2, \qquad d_{t-1} = \mathbf{1}(y_{t-1} < 0)$$

- Parametri: $\theta = (\omega, \alpha, \gamma, \beta)'$
- Vincoli: $\omega > 0,\ \alpha \ge 0,\ \alpha + \gamma \ge 0,\ \beta \ge 0$; stazionarietà $\alpha + \gamma/2 + \beta < 1$
- Varianza non condizionata: $V(y_t) = \dfrac{\omega}{1 - \alpha - \gamma/2 - \beta}$

---

## Step 0 — Setup

- [ ] Scegliere il linguaggio (R / MATLAB / Python)
- [ ] Definire la struttura della repository:
  ```
  src/        funzioni (simulazione, stima, diagnostica)
  scripts/    programmi principali (stima su una serie, Monte Carlo)
  results/    tabelle e grafici prodotti
  report/     rapporto scritto
  ```
- [ ] Fissare i parametri "veri" per le simulazioni (proposta: $\omega=0.05,\ \alpha=0.05,\ \gamma=0.10,\ \beta=0.85$, persistenza 0.95)
- [ ] Fissare il seed per la riproducibilità

## Step 1 — Simulazione (punto a, prima parte)

- [ ] Funzione `simula_gjr(n, theta, burn_in)`:
  - genera $z_t \sim N(0,1)$
  - inizializza $\sigma_1^2 = V(y_t)$ (varianza non condizionata)
  - ricorsione per $\sigma_t^2$ e $y_t$
  - scarta un periodo di burn-in (es. 500 oss.)
  - controlla i vincoli su $\theta$ e restituisce un errore se violati
- [ ] Verifiche: grafico della serie, varianza campionaria vs $V(y_t)$ teorica, ACF di $y_t$ e $y_t^2$, curtosi > 3

## Step 2 — Derivate analitiche

- [ ] Log-verosimiglianza: $\ell(\theta) = -\tfrac12 \sum_t \left[\log 2\pi + \log \sigma_t^2 + y_t^2/\sigma_t^2\right]$
- [ ] Derivate prime di $\sigma_t^2$ (ricorsive), con $x_{t-1} = (1,\ y_{t-1}^2,\ d_{t-1}y_{t-1}^2,\ \sigma_{t-1}^2)'$:
  $$\frac{\partial \sigma_t^2}{\partial \theta} = x_{t-1} + \beta \frac{\partial \sigma_{t-1}^2}{\partial \theta}$$
- [ ] Score: $\dfrac{\partial \ell_t}{\partial \theta} = \dfrac{1}{2\sigma_t^2}\left(\dfrac{y_t^2}{\sigma_t^2} - 1\right)\dfrac{\partial \sigma_t^2}{\partial \theta}$
- [ ] Derivate seconde (Hessiano), anch'esse ricorsive
- [ ] Controllo: confrontare score e Hessiano analitici con quelli numerici (differenze finite)

## Step 3 — Stima (punto a, seconda parte)

- [ ] Funzione `varianza_condizionata(theta, y)` (ricorsione, $\sigma_1^2$ = varianza campionaria)
- [ ] Funzione `loglik(theta, y)` (negativa, da minimizzare)
- [ ] Ottimizzazione vincolata (es. riparametrizzazione o bound + vincolo di stazionarietà), più valori iniziali per robustezza
- [ ] Errori standard:
  - **ML**: $\widehat{V} = A^{-1}$ con $A = -\sum_t \partial^2 \ell_t / \partial\theta\partial\theta'$ (in alternativa OPG $B^{-1}$)
  - **QML** (Bollerslev–Wooldridge): $\widehat{V} = A^{-1} B A^{-1}$, con $B = \sum_t s_t s_t'$
- [ ] Test $t$ di significatività $\hat\theta_i / se(\hat\theta_i)$ e p-value (normale standard), per ML e QML
- [ ] Attenzione: per $\omega, \alpha, \gamma, \beta$ sul bordo ($=0$) l'asintotica standard non vale → commentarlo nel report

## Step 4 — Diagnostica (punto a, terza parte)

- [ ] Serie della varianza condizionata stimata $\hat\sigma_t^2$ (grafico)
- [ ] Residui standardizzati $\hat z_t = y_t / \hat\sigma_t$ (grafico)
- [ ] Test di Jarque–Bera su $\hat z_t$ (statistica, p-value)
- [ ] QQ-plot gaussiano di $\hat z_t$
- [ ] (Extra) Ljung–Box su $\hat z_t$ e $\hat z_t^2$, news impact curve stimata
- [ ] Programma principale che stampa una tabella: stima | se ML | t ML | p ML | se QML | t QML | p QML

## Step 5 — Monte Carlo (punto b)

- [ ] Griglia di numerosità: es. $n \in \{250, 500, 1000, 2000, 5000\}$
- [ ] $R$ replicazioni per ogni $n$ (es. $R = 1000$): simula → stima → salva $\hat\theta$, se, convergenza
- [ ] **Correttezza asintotica** $E(\hat\theta) \to \theta$: distorsione media $\bar{\hat\theta} - \theta$ per $n$; test $t$ su $H_0: E(\hat\theta) = \theta$
- [ ] **Consistenza** $V(\hat\theta) \to 0$: varianza e RMSE per $n$; verificare che $n \cdot V(\hat\theta)$ sia ~ costante (tasso $\sqrt{n}$), grafico log-log
- [ ] **Normalità asintotica**: standardizzare $(\hat\theta - \theta)/se$ e verificare su istogrammi, QQ-plot, Jarque–Bera, Kolmogorov–Smirnov / Shapiro–Wilk
- [ ] Confronto se medie (ML e QML) con la deviazione standard Monte Carlo
- [ ] Copertura degli intervalli di confidenza al 95%
- [ ] (Extra) ripetere con $z_t \sim t$ di Student standardizzata per mostrare che servono gli se QML
- [ ] Salvare risultati in `results/`, tabelle e grafici

## Step 6 — Rapporto scritto (punto c)

- [ ] Introduzione: modello, motivazione (effetto leva), proprietà, vincoli
- [ ] Derivazione di $V(y_t)$ e dei momenti
- [ ] Calcolo completo delle derivate analitiche (score e Hessiano)
- [ ] Descrizione delle funzioni implementate
- [ ] Esempio di stima su una serie simulata (tabelle + grafici diagnostici)
- [ ] Risultati Monte Carlo e commento
- [ ] Conclusioni
- [ ] Listato del programma in file separato (come richiesto dalla consegna)

---

## Divisione dei compiti (da definire)

| Step | Chi | Stato |
|---|---|---|
| 0 Setup | | |
| 1 Simulazione | | |
| 2 Derivate | | |
| 3 Stima | | |
| 4 Diagnostica | | |
| 5 Monte Carlo | | |
| 6 Rapporto | | |
