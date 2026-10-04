# Relay-BP tuning notes

* Full grid (gamma0 in {0.05, 0.125, 0.2, 0.3} x intervals {(-0.24,0.66), (-0.12,0.5), (-0.36,0.8)}),
  300 legs, stop_nconv=1, scored by logical failures (censored decodes count as failures):
  [[72,12,6]] (p=0.004, 20k shots), [[90,8,10]] (p=0.004, 12k), [[108,8,10]] (p=0.004, 12k).
  The interval (-0.12, 0.5) won for every code.
* [[144,12,12]] (p=0.005, 8k shots): the full grid was stopped after 5 configurations for cost; then
  gamma0 in {0.2, 0.3} at interval (-0.12, 0.5). Combined results at interval (-0.12, 0.5):
  gamma0=0.05: 345 failures, 0.125: 346, 0.2: 366, 0.3: 393. gamma0=0.05 and 0.125 are tied within
  noise; 0.125 (the published gross-code value) was selected, overriding the reduced run's pick of 0.2.
* [[288,12,18]] (p=0.005, 3k shots), interval (-0.12, 0.5): gamma0=0.125: 160, 0.2: 167, 0.3: 159
  failures (all failures were non-converged decodes at the 300-leg limit). Tied within noise; 0.3 kept.
* Raw per-configuration rows: relay_tuning_grid.jsonl.txt (full-grid run) and relay_tuning.csv.
