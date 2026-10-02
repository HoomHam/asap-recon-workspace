"""Scratch: temporal-coupling sweep of the 4D least-squares arm on the corrected model (kit v2 scores, thermal units)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c_dyn as DY, _h1c_score as SC
d = DY.Dyn()
for lt in (0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0):
    SC.show(SC.run_arm(f'rx_lt{lt:g}_it15', lt, 1e-3, 15, d=d))
d0 = DY.Dyn(filt=False)
SC.show(SC.run_arm('nofilt_lt0.1_it15', 0.1, 1e-3, 15, d=d0))
