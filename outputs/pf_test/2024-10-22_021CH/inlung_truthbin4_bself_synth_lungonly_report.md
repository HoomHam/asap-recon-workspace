# In-lung error vs fully-sampled truth — 2024-10-22_021CH

rms err and texture err are inside the lung, relative to the truth lung mean. "none" = bin hole pattern only (pure aliasing); "matched" = plus complex noise on the sampled cells scaled so P1 far-bg σ equals the real bin; "real data" = the actual bin recon (motion + noise + aliasing) vs the all-data truth.

| bin | noise | variant | rms err | texture err | corr | bias | lung mean ratio |
|---|---|---|---|---|---|---|---|
| 0 | none | P1 | 0.0803 | 0.0667 | 0.9884 | -0.0259 | 0.9741 |
| 0 | none | PFs | 0.0089 | 0.0071 | 0.9999 | -0.0054 | 0.9946 |
| 0 | none | PFs+k | 0.0263 | 0.0230 | 0.9987 | -0.0097 | 0.9903 |
| 0 | none | PFc | 0.0071 | 0.0071 | 0.9999 | -0.0002 | 0.9998 |
| 0 | matched | P1 | 0.0943 | 0.0817 | 0.9815 | -0.0255 | 0.9745 |
| 0 | matched | PFs | 0.0792 | 0.0766 | 0.9852 | -0.0054 | 0.9946 |
| 0 | matched | PFs+k | 0.0832 | 0.0799 | 0.9833 | -0.0097 | 0.9903 |
| 0 | matched | PFc | 0.1064 | 0.1048 | 0.9736 | +0.0002 | 1.0002 |
| 0 | real data | P1 | 0.4113 | 0.2190 | 0.7121 | -0.2083 | 0.7917 |
| 0 | real data | PFs | 0.4364 | 0.2503 | 0.6976 | -0.1974 | 0.8026 |
| 0 | real data | PFs+k | 0.4300 | 0.2403 | 0.7021 | -0.1995 | 0.8005 |
| 0 | real data | PFc | 0.4284 | 0.2463 | 0.7013 | -0.1923 | 0.8077 |
