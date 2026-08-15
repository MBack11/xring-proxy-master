# XRing Proxy Master

C++ / Gurobi pipeline that **jointly** searches optical ring topologies and shortcut waveguides for wavelength-routed optical networks-on-chip (WR-ONoCs), then verifies real shortcut geometry and assigns wavelengths.

This is my SHK research code building on the **XRing** synthesis method.

## 1. Starting point — the XRing paper

> T. Huang, Z. Li, et al., *XRing: A Crosstalk-Aware Synthesis Method for Wavelength-Routed Optical Ring Routers*, DATE 2023.  
> DOI: [10.23919/DATE56975.2023.10137181](https://doi.org/10.23919/date56975.2023.10137181)

XRing automatically synthesizes optical ring routers: it places network nodes on concentric ring waveguides, adds shortcuts to shorten long paths, and opens the ring for a crossing-free power-distribution network. The paper shows large gains in insertion loss, crosstalk, and laser power versus earlier WR-ONoC rings.

I do **not** host the IEEE publisher PDF here (copyright). Use the DOI above.

## 2. What I work on — improving the search

XRing (and early prototypes) make it hard to treat **ring shape** and **shortcuts** as one optimisable object with a proof of optimality.

This repository implements a **proxy-master loop**:

```
Ring Proposer (proxy MILP, lower bound LB)
    → Quick Estimator (ranking only)
    → Shortcut Solver (true geometry → W_true)
    → Ring Exclude (no-good cuts)
    → compare W* vs LB  (optima when they meet)
    → Wavelength postprocess (conflict graph → λ)
```

- The **master** is optimistic (abstract shortcuts charged Manhattan δ) → valid **lower bound**.
- The **Shortcut Solver** routes real paths on a grid → achievable **upper bound** `W*`.
- When `W* ≈ LB`, the design is **certified optimal** (within ε), not just “best of the runs I tried”.

Formulations and flow charts (my write-ups):

- [`docs/pdf/master_milp.pdf`](docs/pdf/master_milp.pdf) — proxy master MILP
- [`docs/pdf/shortcut_solver_milp.pdf`](docs/pdf/shortcut_solver_milp.pdf) — fixed-ring shortcut MILP
- [`docs/pdf/flow_diagrams.pdf`](docs/pdf/flow_diagrams.pdf) — full M→E→V→C pipeline

## 3. Current status — layout + results

![Example proxy-master layout](figures/layout_plot_PM_ring_sc.png)

Shared-mode proxy vs baseline (excerpt):

| N | Proxy mean W* | Baseline mean W | Mean Δ (W* − baseline) |
|---|--------------:|----------------:|-----------------------:|
| 6 | 20.83 | 21.03 | −0.20 |
| 8 | 24.38 | 25.05 | −0.67 |
| 12 | 26.40 | 28.33 | −1.93 |
| 14 | 26.19 | 27.11 | −0.92 |
| 16 | 27.99 | 30.26 | −2.27 |

Negative Δ means the proxy-master incumbent is **better** (shorter worst-case path) than the baseline on average.

## Layout

| Path | Contents |
|------|----------|
| [`src/`](src/) | C++ sources — proxy master, Stage E, Method D, wavelength postprocess |
| [`docs/pdf/`](docs/pdf/) | MILP + flow documentation |
| [`figures/`](figures/) | Example ring + shortcut layout |

## Role

Student assistant (SHK) research code at TUM, building on XRing under supervision in the optical NoC / photonic interconnect line of work. The DATE paper is by the XRing authors; the proxy-master formulations, loop, and experiments in this repo are my contribution. XRing (DATE 2023) remains © the original authors / IEEE — cite the DOI, do not redistribute the publisher PDF from this repo.
