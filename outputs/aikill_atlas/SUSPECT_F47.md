# ⚠ SUSPECT [F47] — dissolved (dp) pages and dp videos

Generated 2026-08-27 from `d/recon.mat` → `dissolved_phase_magnitude` = |aRBC + i·aTP|.
For sessions where Steve's RBC/TP split ran (86/97 dissolved sessions), that magnitude is NOT
the dissolved image magnitude: the split is ill-conditioned (noise ×1/|sinΔφ|, corr −cosΔφ) and
its phase stop is off target in 60/86 sessions. Gas pages are unaffected.
Regenerate dp pages from the covariance-whitened magnitude (helpers/snr_calc.py, F48) before trusting.
See workspace/canon/facts.md F47/F48, workspace/reference/DP_RBCTP_Split_Noise_Diagnostics.md.
