"""Scratch (2steve/06): compute Steve's calcb b for an 8-channel session with the XeCS numpy replica
at navigator resolution (MS=104, IS=100) and measure the static combine weight W = sum_ch |b_ch|
after his b /= max_ch|b| normalisation. usage: python _calcb_8ch_weight.py <session>"""
import sys, os, json, numpy as np
ASAP=os.path.expanduser('~/Hooman/Work/Codes/2026_ASAP_Recon'); XE=os.path.expanduser('~/Hooman/Work/Codes/2026_XeCS_Recon')
sys.path.insert(0,ASAP); sys.path.insert(0,os.path.join(ASAP,'workspace/helpers/recon')); sys.path.insert(0,os.path.join(XE,'recon'))
from selftest_steve_tyger import read_input_mrd, exclude_mask
from steve_kernel_numpy import steve_recon
from scipy.optimize import curve_fit
from gtypes import gvar, imgtype, graddir
from raw import traj, raw as R
s=sys.argv[1]; run=f'/Volumes/HoomHamExt/AIkill_Dynamic/{s}/d'
OUT=os.path.join(ASAP,'workspace/outputs/calcb_imprint_2026-09-24'); os.makedirs(OUT,exist_ok=True)
header,arrs=read_input_mrd(os.path.join(run,'input.mrd'))
ul={p.name:p.value for p in header.user_parameters.user_parameter_long}; ud={p.name:p.value for p in header.user_parameters.user_parameter_double}
g=gvar(); g.MS=104; g.IS=100; g.gplb=int(ul.get('gplb',300)); kp=int(ul.get('killpts',2))
meta={'TR':ud['TR'],'TE':ud['TE'],'DPoff':ud['DPoff'],'dtdyn':ud['dtdyn'],'dtspec':ud['dtspec'],'numspec':int(ul['numspec'])}
gt=traj(); gt.killpts=kp; gt.load_traj_from_array(arrs['gp'],arrs['dp'],int(ul.get('nusimg',32)))
dyn=arrs['dyn'][:,kp:,:]; ref=arrs['ref'][:,kp:,:] if arrs['ref'] is not None else None
gr=R(); gr.load_from_arr(gt,ref,dyn,arrs['pneumo'],'mrd_siemens',meta); gt.rescale_to_MS(g.MS,g.IS)
tg=np.stack([np.asarray(gt.gettraj(imgtype.GPDYN,d)) for d in (graddir.X,graddir.Y,graddir.Z)],axis=1)
npts=gr.npts; nch=gr.nch; em=exclude_mask(gr); print('nch',nch,'npts',npts,'nilv',gr.ntotalilvs,flush=True)
MS,IS=g.MS,g.IS; MSc=int(MS/2+.1); ll=int(MSc-IS/2+.1)
b=np.zeros((nch,IS,IS,IS),complex)
for ich in range(nch):
    data=np.reshape(np.asarray(gr.getimg(imgtype.GPDYN))[:,ich,:],npts*gr.ntotalilvs,order='F')*em
    F=steve_recon(tg,data,npts,MS=MS,IS=MS,smoothing=g.gplb)
    F=F/np.mean(np.abs(F[0:10,0:10,0:10])); b[ich]=np.conj(F[ll:ll+IS,ll:ll+IS,ll:ll+IS]); print('ch',ich,'done',flush=True)
# calcb phase fix per channel (results.py:126-151)
x,y,z=np.meshgrid(np.linspace(-1,1,IS),np.linspace(-1,1,IS),np.linspace(-1,1,IS)); u=np.ones((IS,IS,IS))
basis_full=[u,x,y,z,x*x,y*y,z*z,x*y,y*z,x*z]
for ich in range(nch):
    bm=(np.abs(b[ich])>np.mean(np.abs(b[ich][0:10,0:10,0:10])*10)).astype(int); sel=bm.reshape(-1)==1
    fbp=b[ich].reshape(-1)[sel]; ma=np.mean(fbp); ma/=abs(ma); fbp=fbp*np.conj(ma)
    bas=[bb.reshape(-1)[sel] for bb in basis_full]
    def ff(_, *c): return sum(ci*bi for ci,bi in zip(c,bas))
    p,_=curve_fit(ff,bas[0],np.angle(fbp),np.zeros(10))
    ph=sum(pi*bb for pi,bb in zip(p,basis_full))
    b[ich]=bm*b[ich]+(1-bm)*ma*np.abs(b[ich])*np.exp(1j*ph)
bmax=np.max(np.abs(b),0); b=b/bmax                     # results.py:154
W=np.sum(np.abs(b),0)                                  # effective static weight of sum_ch real(F_ch b_ch)
Fmag=np.sqrt(np.sum(np.abs(b*bmax)**2,0)); noise=np.mean(Fmag[0:10,0:10,0:10]); lung=Fmag>5*noise
res=dict(session=s,nch=nch,MS=MS,IS=IS,lung_vox=int(lung.sum()),W_lung_min=float(W[lung].min()),W_lung_p5=float(np.percentile(W[lung],5)),
         W_lung_median=float(np.median(W[lung])),W_lung_p95=float(np.percentile(W[lung],95)),W_lung_max=float(W[lung].max()),
         W_bg_median=float(np.median(W[~lung])))
print(json.dumps(res,indent=1)); json.dump(res,open(f'{OUT}/W_{s}.json','w'),indent=1); np.save(f'{OUT}/W_{s}.npy',W.astype('float32'))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
yc=np.argmax(lung.sum((0,2))); fig,ax=plt.subplots(1,3,figsize=(13,4.2))
im=ax[0].imshow(Fmag[:,yc,:],cmap='gray'); ax[0].set_title(f'{s}: cycle-average |F| (SoS over 8 ch)'); 
im1=ax[1].imshow(W[:,yc,:],cmap='viridis',vmin=1,vmax=nch); plt.colorbar(im1,ax=ax[1]); ax[1].set_title('W = Σ_ch |b_ch| (static gain on every bin)')
im2=ax[2].imshow(np.where(lung[:,yc,:],W[:,yc,:]/np.median(W[lung]),np.nan),cmap='coolwarm',vmin=0.5,vmax=1.5); plt.colorbar(im2,ax=ax[2]); ax[2].set_title('W / median(W in lung), lung only')
for a in ax: a.axis('off')
plt.tight_layout(); plt.savefig(f'{OUT}/06_W_map_{s}.png',dpi=110); print('fig saved',flush=True)
