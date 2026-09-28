# In-lung error vs fully-sampled truth — 2024-10-22_021CH

rms err and texture err are inside the lung, relative to the truth lung mean. "none" = bin hole pattern only (pure aliasing); "matched" = plus complex noise on the sampled cells scaled so P1 far-bg σ equals the real bin; "real data" = the actual bin recon (motion + noise + aliasing) vs the all-data truth.

| bin | noise | variant | rms err | texture err | corr | bias | lung mean ratio |
|---|---|---|---|---|---|---|---|
| 0 | none | P1 | 0.0100 | 0.0081 | 0.9998 | -0.0028 | 0.9972 |
| 0 | none | PFs | 0.0996 | 0.0801 | 0.9871 | +0.0267 | 1.0267 |
| 0 | none | PFs+k | 0.0970 | 0.0777 | 0.9877 | +0.0260 | 1.0260 |
| 0 | none | PFc | 0.1048 | 0.0917 | 0.9845 | +0.0271 | 1.0271 |
| 0 | matched | P1 | 0.0499 | 0.0476 | 0.9941 | -0.0028 | 0.9972 |
| 0 | matched | PFs | 0.1266 | 0.1103 | 0.9756 | +0.0264 | 1.0264 |
| 0 | matched | PFs+k | 0.1247 | 0.1086 | 0.9761 | +0.0257 | 1.0257 |
| 0 | matched | PFc | 0.1503 | 0.1404 | 0.9633 | +0.0271 | 1.0271 |
| 0 | real data | P1 | 0.4104 | 0.2234 | 0.7136 | -0.2035 | 0.7965 |
| 0 | real data | PFs | 0.4358 | 0.2539 | 0.6990 | -0.1926 | 0.8074 |
| 0 | real data | PFs+k | 0.4293 | 0.2440 | 0.7035 | -0.1946 | 0.8054 |
| 0 | real data | PFc | 0.4279 | 0.2500 | 0.7026 | -0.1874 | 0.8126 |
| 8 | none | P1 | 0.0077 | 0.0061 | 0.9999 | -0.0022 | 0.9978 |
| 8 | none | PFs | 0.0945 | 0.0768 | 0.9876 | +0.0209 | 1.0209 |
| 8 | none | PFs+k | 0.0927 | 0.0750 | 0.9880 | +0.0204 | 1.0204 |
| 8 | none | PFc | 0.1051 | 0.0919 | 0.9844 | +0.0272 | 1.0272 |
| 8 | matched | P1 | 0.0871 | 0.0833 | 0.9821 | -0.0012 | 0.9988 |
| 8 | matched | PFs | 0.1704 | 0.1571 | 0.9515 | +0.0215 | 1.0215 |
| 8 | matched | PFs+k | 0.1694 | 0.1562 | 0.9517 | +0.0211 | 1.0211 |
| 8 | matched | PFc | 0.2185 | 0.2097 | 0.9220 | +0.0282 | 1.0282 |
| 8 | real data | P1 | 1.2012 | 0.4828 | 0.8286 | +0.9275 | 1.9275 |
| 8 | real data | PFs | 1.3694 | 0.6891 | 0.7958 | +0.9764 | 1.9764 |
| 8 | real data | PFs+k | 1.3356 | 0.6404 | 0.8044 | +0.9692 | 1.9692 |
| 8 | real data | PFc | 1.3568 | 0.6746 | 0.7986 | +0.9790 | 1.9790 |
