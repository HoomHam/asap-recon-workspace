# In-lung error vs fully-sampled truth — 2024-10-22_021CH

rms err and texture err are inside the lung, relative to the truth lung mean. "none" = bin hole pattern only (pure aliasing); "matched" = plus complex noise on the sampled cells scaled so P1 far-bg σ equals the real bin; "real data" = the actual bin recon (motion + noise + aliasing) vs the all-data truth.

| bin | noise | variant | rms err | texture err | corr | bias | lung mean ratio |
|---|---|---|---|---|---|---|---|
| 0 | none | P1 | 0.0100 | 0.0080 | 0.9998 | -0.0030 | 0.9970 |
| 0 | none | PFs | 0.0984 | 0.0822 | 0.9862 | +0.0125 | 1.0125 |
| 0 | none | PFs+k | 0.0960 | 0.0800 | 0.9868 | +0.0118 | 1.0118 |
| 0 | none | PFc | 0.0666 | 0.0581 | 0.9927 | +0.0170 | 1.0170 |
| 0 | matched | P1 | 0.0494 | 0.0471 | 0.9941 | -0.0025 | 0.9975 |
| 0 | matched | PFs | 0.1245 | 0.1106 | 0.9753 | +0.0125 | 1.0125 |
| 0 | matched | PFs+k | 0.1227 | 0.1090 | 0.9757 | +0.0118 | 1.0118 |
| 0 | matched | PFc | 0.1377 | 0.1326 | 0.9633 | +0.0175 | 1.0175 |
| 0 | real data | P1 | 0.4113 | 0.2190 | 0.7121 | -0.2083 | 0.7917 |
| 0 | real data | PFs | 0.4364 | 0.2503 | 0.6976 | -0.1974 | 0.8026 |
| 0 | real data | PFs+k | 0.4300 | 0.2403 | 0.7021 | -0.1995 | 0.8005 |
| 0 | real data | PFc | 0.4284 | 0.2463 | 0.7013 | -0.1923 | 0.8077 |
