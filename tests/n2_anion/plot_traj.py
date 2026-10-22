import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
EV=27.211386245988; INK='#0b0b0b'; INK2='#52514e'; GRID='#e4e3df'; BLUE='#2a78d6'; GRAY='#b9b8b3'
plt.rcParams.update({'font.size':11,'axes.edgecolor':INK2,'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2,
                     'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb'})
d=np.load('traj_big_r3.0_500_1202.npz'); eta,t,h=d['etas'],d['traj'],d['h']
res=json.load(open('res_r3.json')); sel=[c['eta'] for c in res['stationary_candidates']]
low=[m for m in range(len(h)) if h[m]<0.4]; R=4
fig,ax=plt.subplots(1,3,figsize=(16,4.8),constrained_layout=True)
# (a) complex plane
for m in low:
    if m!=R: ax[0].plot(t[:,m].real*EV,-2*t[:,m].imag*EV,color=GRAY,lw=1.5)
ax[0].plot(t[:,R].real*EV,-2*t[:,R].imag*EV,color=BLUE,lw=2.5)
for c in res['stationary_candidates']:
    ax[0].plot(c['raw']['energy_ev'],c['raw']['width_ev'],'o',ms=8,color=BLUE,mec='#fcfcfb',mew=2)
ax[0].set(xlim=(0,10),ylim=(0,8),xlabel='Re E  (eV)',ylabel='Γ = −2 Im E  (eV)',title='(a) πg trajectories, η = 1e-5 → 1e-2')
ax[0].annotate('resonance (stays near real axis)',xy=(3.75,1.3),xytext=(4.6,3.2),color=INK,arrowprops=dict(arrowstyle='-',color=INK2))
k=np.argmin(abs(-2*t[:,1].imag*EV-5)); ax[0].annotate('discretized continuum',xy=(t[k,1].real*EV,5),xytext=(5.6,6.3),color=INK2,arrowprops=dict(arrowstyle='-',color=INK2))
# (b) speed
sp=np.abs(np.gradient(t,np.log(eta),axis=0))*EV
for m in low:
    if m!=R: ax[1].loglog(eta,sp[:,m],color=GRAY,lw=1.5)
ax[1].loglog(eta,sp[:,R],color=BLUE,lw=2.5)
for c in res['stationary_candidates']:
    ax[1].plot(c['eta'],c['logarithmic_speed']*EV,'o',ms=8,color=BLUE,mec='#fcfcfb',mew=2)
ax[1].set(xlabel='η',ylabel='|dE/d ln η|  (eV)',title='(b) trajectory speed; dots = stationary points')
# (c) onset dependence, resonance branch only (ordinal blue ramp)
for r,col in zip([2,3,4,5],['#86b6ef','#3987e5','#1c5cab','#0d366b']):
    tr=np.load(f'traj_big_r{r}.0_500_1202.npz')['traj'][:,R]
    ax[2].plot(tr.real*EV,-2*tr.imag*EV,color=col,lw=2,label=f'r_cut = {r} bohr')
    c=json.load(open(f'res_r{r}.json'))['selected_candidate']
    ax[2].plot(c['raw']['energy_ev'],c['raw']['width_ev'],'o',ms=8,color=col,mec='#fcfcfb',mew=2)
ax[2].set(xlabel='Re E  (eV)',ylabel='Γ  (eV)',title='(c) resonance branch vs CAP onset')
ax[2].legend(frameon=False,labelcolor=INK)
for a in ax: a.grid(color=GRID,lw=0.8); a.title.set_color(INK)
fig.suptitle('N2⁻ ²Πg, CAP/Koopmans, aug-cc-pVTZ+4p5d, smooth Voronoi CAP (r_cut = 3 bohr in a, b)',color=INK)
fig.savefig('n2m_trajectories.png',dpi=160)
