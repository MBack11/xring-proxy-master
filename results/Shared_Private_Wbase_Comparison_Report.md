# Shared vs Private Benchmark Report

## Study Setup

- Runner: `snr_crossing_benchmark` (study outputs under `benchmarks/results/snr_study/`)
- Baseline ring: Method B @ 90s (warm-start A); fallback A / `B_timeout`
- Shared W_base: Method D joint shortcuts
- Private W_base: WC-greedy Method A+Shortcuts (no Method-D MILP)
- Proxy Stage V: Method D (both modes)
- Assertion note: **W\* / W_base may differ between Shared and Private** (private changes MIP flow coupling). Wavelength/SNR are postprocess and do not redefine W.
- Seeds: N=6/8 → 20; N=12/14/16 → 10 (report count-totals ×2 for N≥12)
- Wall budget (proxy+base combined): N6=240s, N8=360s, N≥12=600s

## Source Files

- Shared seeds: `benchmarks/results/snr_study/seeds_shared.csv`
- Private seeds: `benchmarks/results/snr_study/seeds_private.csv`
- Shared pairs: `benchmarks/results/snr_study/pairs_shared.csv`
- Private pairs: `benchmarks/results/snr_study/pairs_private.csv`

## Table 1: Shared (Proxy vs W_base)

### By N

| N | Seeds (reported) | Proxy mean W* | W_base mean W | Mean Delta (W* - W_base) | Proxy mean time (s) | W_base mean time (s) | Mean total time (s) | Proxy mean cross_pairs | W_base mean cross_pairs | Proxy mean signals_affected/pair | W_base mean signals_affected/pair | Base fail rate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 20 | 20.825 | 21.028 | -0.203 | 1.2 | 0.7 | 1.9 | 0.20 | 0.35 | 0.20 | 0.50 | 0.0% |
| 8 | 20 | 24.375 | 25.048 | -0.673 | 6.1 | 3.7 | 9.8 | 0.70 | 0.60 | 1.95 | 1.25 | 0.0% |
| 12 | 10 (×2→20) | 26.400 | 28.329 | -1.929 | 18.8 | 27.1 | 45.9 | 1.00 | 0.80 | 3.40 | 0.85 | 0.0% |
| 14 | 10 (×2→20) | 26.190 | 27.112 | -0.922 | 162.1 | 83.2 | 245.2 | 1.50 | 1.00 | 3.07 | 3.00 | 0.0% |
| 16 | 10 (×2→20) | 27.990 | 30.256 | -2.266 | 196.6 | 159.3 | 355.9 | 0.90 | 1.20 | 4.40 | 3.90 | 0.0% |

### Delta Distribution

| N | % (W* < W_base) | % (W* = W_base) | % (W* > W_base) | Median Delta W | P90 Delta W | Min Delta W | Max Delta W |
|---|---:|---:|---:|---:|---:|---:|---:|
| 6 | 20.0 | 80.0 | 0.0 | 0.000 | 0.000 | -2.000 | 0.000 |
| 8 | 30.0 | 70.0 | 0.0 | 0.000 | 0.000 | -6.500 | 0.000 |
| 12 | 30.0 | 70.0 | 0.0 | 0.000 | 0.000 | -9.500 | 0.000 |
| 14 | 40.0 | 50.0 | 10.0 | 0.000 | 0.268 | -7.000 | 2.676 |
| 16 | 40.0 | 60.0 | 0.0 | 0.000 | 0.000 | -9.655 | 0.000 |

## Table 2: Private (Proxy vs W_base)

### By N

| N | Seeds (reported) | Proxy mean W* | W_base mean W | Mean Delta (W* - W_base) | Proxy mean time (s) | W_base mean time (s) | Mean total time (s) | Proxy mean cross_pairs | W_base mean cross_pairs | Proxy mean signals_affected/pair | W_base mean signals_affected/pair | Base fail rate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 20 | 20.864 | 20.125 | 0.739 | 1.2 | 1.2 | 2.4 | 0.05 | 0.10 | 0.15 | 0.20 | 0.0% |
| 8 | 20 | 24.890 | 24.225 | 0.665 | 6.6 | 3.6 | 10.2 | 0.80 | 0.50 | 1.55 | 0.15 | 0.0% |
| 12 | 10 (×2→20) | 28.000 | 28.950 | -0.950 | 143.6 | 21.9 | 165.5 | 0.60 | 0.60 | 1.30 | 0.30 | 0.0% |
| 14 | 10 (×2→20) | 31.140 | 35.400 | -4.260 | 351.9 | 57.8 | 409.7 | 1.30 | 0.50 | 1.65 | 0.40 | 0.0% |
| 16 | 10 (×2→20) | 37.690 | 34.650 | 3.040 | 978.5 | 446.7 | 1425.2 | 0.80 | 0.60 | 1.20 | 0.30 | 0.0% |

### Delta Distribution

| N | % (W* < W_base) | % (W* = W_base) | % (W* > W_base) | Median Delta W | P90 Delta W | Min Delta W | Max Delta W |
|---|---:|---:|---:|---:|---:|---:|---:|
| 6 | 25.0 | 40.0 | 35.0 | 0.000 | 6.650 | -4.500 | 9.000 |
| 8 | 25.0 | 30.0 | 45.0 | 0.000 | 7.650 | -7.100 | 11.000 |
| 12 | 40.0 | 20.0 | 40.0 | 0.000 | 5.250 | -14.000 | 12.000 |
| 14 | 80.0 | 0.0 | 20.0 | -3.800 | 2.000 | -14.500 | 15.500 |
| 16 | 40.0 | 0.0 | 60.0 | 1.750 | 12.070 | -8.500 | 16.300 |

## Table 3: Proxy Shared vs Proxy Private

| N | Common seeds | Mean Delta W* (private - shared) | Median Delta W* | P90 Delta W* | % (Delta > 0) | % (Delta = 0) | % (Delta < 0) | Mean Delta proxy time (s) | Mean Delta cross_pairs | Mean Delta signals_affected/pair |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 20 | 0.039 | 0.000 | 0.038 | 10.0 | 90.0 | 0.0 | 0.0 | -0.15 | -0.05 |
| 8 | 20 | 0.515 | 0.000 | 1.550 | 35.0 | 65.0 | 0.0 | 0.5 | 0.10 | -0.40 |
| 12 | 10 (×2→20) | 1.600 | 0.000 | 4.750 | 40.0 | 60.0 | 0.0 | 124.9 | -0.40 | -2.10 |
| 14 | 10 (×2→20) | 4.950 | 0.700 | 10.000 | 60.0 | 40.0 | 0.0 | 189.8 | -0.20 | -1.42 |
| 16 | 10 (×2→20) | 9.700 | 11.500 | 17.800 | 90.0 | 10.0 | 0.0 | 781.9 | -0.10 | -3.20 |

## Runtime Breakdown

| N | Mode | Proxy time min/med/mean/p90/max (s) | W_base time min/med/mean/p90/max (s) | Total time min/med/mean/p90/max (s) |
|---|---|---|---|---|
| 6 | shared | 0.3/1.1/1.2/1.9/2.3 | 0.2/0.6/0.7/1.2/1.6 | 0.6/1.7/1.9/3.1/3.8 |
| 6 | private | 0.3/1.1/1.2/1.6/3.6 | 0.2/1.0/1.2/2.1/2.2 | 0.8/1.9/2.4/3.6/5.0 |
| 8 | shared | 2.3/3.8/6.1/11.1/22.2 | 2.2/3.7/3.7/4.9/5.0 | 5.8/8.0/9.8/15.0/27.2 |
| 8 | private | 2.7/3.7/6.6/14.3/27.3 | 1.3/3.4/3.6/4.6/5.5 | 4.5/7.8/10.2/18.9/30.4 |
| 12 | shared | 15.0/18.6/18.8/22.3/23.7 | 18.9/23.3/27.1/34.5/55.1 | 36.5/41.6/45.9/54.0/73.7 |
| 12 | private | 12.8/60.6/143.6/421.7/492.5 | 7.6/20.4/21.9/27.6/52.6 | 36.0/79.3/165.5/433.6/500.1 |
| 14 | shared | 34.3/50.4/162.1/191.2/1138.3 | 42.9/75.5/83.2/133.0/133.6 | 77.1/127.8/245.2/319.1/1214.6 |
| 14 | private | 28.1/216.5/351.9/973.5/1065.1 | 10.9/38.9/57.8/106.2/116.4 | 54.2/291.7/409.7/1011.5/1163.5 |
| 16 | shared | 60.8/131.8/196.6/350.2/464.7 | 66.2/165.7/159.3/183.3/215.0 | 144.4/329.2/355.9/513.6/620.2 |
| 16 | private | 65.9/1026.3/978.5/1790.3/1997.2 | 15.4/283.6/446.7/959.3/973.4 | 171.6/1434.5/1425.2/2188.5/2970.6 |

## Notes

- Private W_base uses WC-greedy single-shortcut insertion; it is a heuristic (order may be suboptimal) but matches private/single-owner semantics without a second MILP.
- For N=12/14/16, means/medians/percentiles use the measured 10 seeds; only count-like totals are scaled ×2 when labeled.

