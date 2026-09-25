"""Scratch: rerun Steve's spectral fit (raw.load_from_arr, CPU) and compare the RBC/TP split's
assumed phase separation (2*pi*df*TEeff, i.e. the t=0 intercept) with the phase difference
measured at the first spectral sample (t = TE + killpts*dtspec). Evidence for 2steve/05.
usage: python _f59_check.py <session> [<session> ...]   (sessions under /Volumes/HoomHamExt/AIkill_Dynamic)
writes workspace/outputs/f59_te_term_2026-09-24/{f59.json, 05_rbctp_split_te_term.png}
"""
import sys, json, numpy as np, mrd, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon')
from gtypes import gvar
from raw import traj, raw
OUT='/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/outputs/f59_te_term_2026-09-24'
def wrap(d): return (d + 180) % 360 - 180
def run(p):
    g=gvar(); gt=traj(); gr=raw()
    ref=dyn=pne=gtr=dtr=None
    with mrd.BinaryMrdReader(p) as r:
        h=r.read_header()
        for it in r.read_data():
            if isinstance(it,(mrd.StreamItem.NdArrayDouble,mrd.StreamItem.NdArrayFloat)):
                m=it.value.meta
                if m.get('gas_phase_trajectory'): gtr=it.value.data
                elif m.get('dissolved_phase_trajectory'): dtr=it.value.data
                elif m.get('pneumotach'): pne=it.value.data
            elif isinstance(it,mrd.StreamItem.NdArrayComplexFloat):
                m=it.value.meta
                if m.get('reference_acquisition'): ref=it.value.data
                elif m.get('dynamic_acquisition'): dyn=it.value.data
    ul={q.name:q.value for q in h.user_parameters.user_parameter_long}
    ud={q.name:q.value for q in h.user_parameters.user_parameter_double}
    kp=int(ul.get('killpts',2)); gt.killpts=kp
    gt.load_traj_from_array(gtr,dtr,int(ul.get('nusimg',32)))
    dyn=dyn[:,kp:,:]; ref=ref[:,kp:,:] if ref is not None else None
    meta={'TR':ud['TR'],'TE':ud['TE'],'DPoff':ud['DPoff'],'dtdyn':ud['dtdyn'],'dtspec':ud['dtspec'],'numspec':int(ul['numspec'])}
    gr.load_from_arr(gt,ref,dyn,pne,'mrd_siemens',meta)
    if not len(gr.fRBC): return None
    df=gr.fRBC[0]-gr.fTP[0]
    dphi=np.degrees(np.array(gr.RBCphase)-np.array(gr.TPphase))
    tt=(np.array(gr.TEphase)+kp/gt.spectBW)*1e3   # true time of each pseudo-TE sample, ms (killpts already trimmed)
    d=dict(session=p.split('/')[-3], TE_ms=gr.TE*1e3, df_Hz=df, ratio=gr.RBCTPratio[0],
           intercept_deg=np.degrees(gr.deltaphase), TEeff_us=gr.TEeff*1e6, TEeff2_us=gr.TEeff2*1e6,
           model_deg=wrap(np.degrees(2*np.pi*df*gr.TEeff)), measured_first_deg=wrap(dphi[0]),
           TEterm_deg=np.degrees(2*np.pi*df*gr.TE), dphi_deg=list(map(float,dphi)), t_ms=list(map(float,tt)),
           slope_deg_per_sample=float(np.degrees(2*np.pi*df/gt.spectBW)))
    d['gain_model']=1/abs(np.sin(np.radians(d['model_deg']))); d['gain_measured']=1/abs(np.sin(np.radians(d['measured_first_deg'])))
    d['fRBC_Hz']=float(gr.fRBC[0]); d['fTP_Hz']=float(gr.fTP[0]); d['DPoff_ppm']=ud['DPoff']; d['sspect']=[float(v) for v in np.abs(gr.sspect[0])]; d['sfreq']=[float(v) for v in gr.sspectfreq[0]]
    print(f"{d['session']}: DPoff={ud['DPoff']:.0f}ppm fRBC={gr.fRBC[0]:.0f} fTP={gr.fTP[0]:.0f} TE={d['TE_ms']:.2f} df={df:.0f}Hz ratio={d['ratio']:.3f} intercept={d['intercept_deg']:.1f} "
          f"model={d['model_deg']:.1f} (gain {d['gain_model']:.1f}) measured={d['measured_first_deg']:.1f} (gain {d['gain_measured']:.2f}) "
          f"gap={wrap(d['measured_first_deg']-d['model_deg']):.1f} 2pi*df*TE={d['TEterm_deg']:.1f}")
    return d
res=[r for r in (run(f'/Volumes/HoomHamExt/AIkill_Dynamic/{s}/d/input.mrd') for s in sys.argv[1:]) if r]
json.dump(res, open(f'{OUT}/f59.json','w'), indent=1)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,2,figsize=(12,4.6))
for d in res:
    ph=np.array(d['dphi_deg']); ph=np.degrees(np.unwrap(np.radians(ph))); ph=ph-360*np.round((ph[0]-d['measured_first_deg'])/360)
    t=np.array(d['t_ms'])
    l,=ax[0].plot(t,ph,'o-',ms=4,label=f"{d['session']} (Δf={d['df_Hz']:.0f} Hz)")
    ax[0].plot([0,t[0]],[d['intercept_deg']+360*np.round((ph[0]-d['intercept_deg']-d['TEterm_deg'])/360), ph[0]],'--',color=l.get_color(),alpha=.6)
ax[0].axvline(res[0]['TE_ms'],color='k',lw=.8,ls=':'); ax[0].text(res[0]['TE_ms']+.01,ax[0].get_ylim()[0]+5,'TE (k0 of image)',fontsize=8)
ax[0].axhline(90,color='gray',lw=.8,ls='--'); ax[0].set_xlabel('time after RF centre (ms)'); ax[0].set_ylabel('φRBC − φTP (deg)')
ax[0].set_title('fitted RBC−TP phase vs pseudo-TE; dashed = linear extrapolation to t=0 (Steve\'s intercept)'); ax[0].legend(fontsize=7); ax[0].set_xlim(0,None)
names=[d['session'][-5:] for d in res]; x=np.arange(len(res))
ax[1].bar(x-0.2,[d['model_deg'] for d in res],0.4,label='Δφ used by split = 2πΔf·TEeff (t=0 intercept)')
ax[1].bar(x+0.2,[d['measured_first_deg'] for d in res],0.4,label='Δφ measured at first sample (t≈TE)')
ax[1].plot(x,[d['TEterm_deg'] for d in res],'k_',ms=18,label='2πΔf·TE (missing term)')
ax[1].axhline(90,color='gray',lw=.8,ls='--'); ax[1].set_xticks(x); ax[1].set_xticklabels(names); ax[1].set_ylabel('deg'); ax[1].legend(fontsize=7)
ax[1].set_title('noise gain 1/|sin Δφ|: ' + ', '.join(f"{d['gain_model']:.1f}→{d['gain_measured']:.2f}" for d in res),fontsize=9)
plt.tight_layout(); plt.savefig(f'{OUT}/05_rbctp_split_te_term.png',dpi=130); print('fig saved')
