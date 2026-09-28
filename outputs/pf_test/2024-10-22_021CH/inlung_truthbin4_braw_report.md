# In-lung error vs fully-sampled truth — 2024-10-22_021CH

rms err and texture err are inside the lung, relative to the truth lung mean. "none" = bin hole pattern only (pure aliasing); "matched" = plus complex noise on the sampled cells scaled so P1 far-bg σ equals the real bin; "real data" = the actual bin recon (motion + noise + aliasing) vs the all-data truth.

| bin | noise | variant | rms err | texture err | corr | bias | lung mean ratio |
|---|---|---|---|---|---|---|---|
| 0 | none | P1 | 0.0100 | 0.0082 | 0.9998 | -0.0029 | 0.9971 |
| 0 | none | PFs | 0.0971 | 0.0789 | 0.9868 | +0.0188 | 1.0188 |
| 0 | none | PFs+k | 0.0947 | 0.0766 | 0.9874 | +0.0180 | 1.0180 |
| 0 | none | PFc | 0.1083 | 0.0968 | 0.9824 | +0.0259 | 1.0259 |
| 0 | matched | P1 | 0.0496 | 0.0473 | 0.9941 | -0.0025 | 0.9975 |
| 0 | matched | PFs | 0.1240 | 0.1087 | 0.9757 | +0.0188 | 1.0188 |
| 0 | matched | PFs+k | 0.1221 | 0.1071 | 0.9761 | +0.0181 | 1.0181 |
| 0 | matched | PFc | 0.1571 | 0.1485 | 0.9588 | +0.0263 | 1.0263 |
| 0 | real data | P1 | 0.4112 | 0.2235 | 0.7124 | -0.2058 | 0.7942 |
| 0 | real data | PFs | 0.4365 | 0.2538 | 0.6979 | -0.1949 | 0.8051 |
| 0 | real data | PFs+k | 0.4300 | 0.2440 | 0.7024 | -0.1969 | 0.8031 |
| 0 | real data | PFc | 0.4285 | 0.2500 | 0.7016 | -0.1897 | 0.8103 |
