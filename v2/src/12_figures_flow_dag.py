"""Figure 1 (selection flow chart) and Figure 2 (conceptual DAG) of the manuscript.
Counts: mortality source file (260,219 neonatal deaths; 196 with a state-only residence code), v2/data/dataset_report.txt and
the analysis tables. The font is Calibri when available (e.g. from an Office installation), otherwise a fallback sans-serif."""
import glob
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib import font_manager as fm
for _f in glob.glob('/Applications/Microsoft Excel.app/Contents/Resources/DFonts/Calibri*.ttf'):
    fm.fontManager.addfont(_f)
names = {f.name for f in fm.fontManager.ttflist}
FONT = next((n for n in ('Calibri','Carlito','Arial','Helvetica','DejaVu Sans') if n in names), 'DejaVu Sans')
plt.rcParams['font.family'] = FONT
print('fonte usada:', FONT)
OUTFIG = Path(__file__).resolve().parents[1] / 'outputs' / 'figures'
OUTFIG.mkdir(parents=True, exist_ok=True)
FS = 12   # PPE: Calibri 12 pt at design size (figures are reduced for print)
GREY = '#e6e6e6'

def box(ax, x, y, w, h, text, dashed=False, fc='white', bold=False):
    p = FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle='round,pad=0.0,rounding_size=0.012',
                       linewidth=1.0, edgecolor='black', facecolor=fc, linestyle='--' if dashed else '-', zorder=2)
    ax.add_patch(p)
    ax.text(x, y, text, ha='center', va='center', fontsize=FS, fontweight='bold' if bold else 'normal', zorder=3, linespacing=1.25)
    return (p, (x, y))

def arrow(ax, a, b, rad=0.0):
    (pa, ca), (pb, cb) = a, b
    ax.annotate('', xy=cb, xytext=ca, zorder=1,
                arrowprops=dict(arrowstyle='-|>', lw=1.0, color='black', patchA=pa, patchB=pb,
                                shrinkA=0, shrinkB=0, connectionstyle=f'arc3,rad={rad}'))

def side_arrow(ax, x_from, x_to, y):
    ax.annotate('', xy=(x_to, y), xytext=(x_from, y), zorder=1,
                arrowprops=dict(arrowstyle='-|>', lw=1.0, color='black', shrinkA=0, shrinkB=0))

# ---------------- Figure 1: flow chart ----------------
fig, ax = plt.subplots(figsize=(10.0, 11.0)); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
src = box(ax, 0.30, 0.93, 0.50, 0.11, 'Source records, Brazil, 2014-2024\nNeonatal deaths (SIM): 260,219\nLive births (Sinasc): 30,463,522')
exA = box(ax, 0.79, 0.81, 0.40, 0.15, 'Without an identifiable municipality\n196 deaths (18 state-only codes)\n980 births (23 state-only codes)', fc=GREY)
idn = box(ax, 0.30, 0.69, 0.50, 0.11, 'With an identified municipality\n260,023 deaths\n30,462,542 live births')
exB = box(ax, 0.79, 0.575, 0.40, 0.15, 'Five municipalities\nwithout MHDI data\n58 deaths; 7,084 births', fc=GREY)
ana = box(ax, 0.30, 0.45, 0.50, 0.14, 'Analysed set\n5,565 municipalities\n259,965 neonatal deaths\n30,455,458 live births')
prim = box(ax, 0.17, 0.25, 0.30, 0.13, 'Primary model\nBYM2 with restricted\nspatial regression\n2014-2024')
per = box(ax, 0.58, 0.25, 0.32, 0.13, 'Period-specific models\n2014-2019: 152,418 deaths\n2020-2024: 107,547 deaths')
sim = box(ax, 0.17, 0.06, 0.32, 0.11, 'Simulation with known residual\nrisk, from the primary model\n40 datasets with residual risk\n12 datasets without')
arrow(ax, src, idn); side_arrow(ax, 0.30, 0.59, 0.81)
arrow(ax, idn, ana); side_arrow(ax, 0.30, 0.59, 0.575)
arrow(ax, ana, prim); arrow(ax, ana, per); arrow(ax, prim, sim)
fig.savefig(str(OUTFIG/'figure1_selection_flow.png'), dpi=300, bbox_inches='tight'); fig.savefig(str(OUTFIG/'figure1_selection_flow.pdf'), bbox_inches='tight'); plt.close(fig)

# ---------------- Figure 2: DAG ----------------
fig, ax = plt.subplots(figsize=(10.0, 9.4)); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
ctx = box(ax, 0.50, 0.94, 0.38, 0.08, 'Geographical and historical\ncontext (neighbourhood)')
un  = box(ax, 0.14, 0.77, 0.22, 0.12, 'Unmeasured\narea-level risk\nfactors', dashed=True)
st  = box(ax, 0.38, 0.77, 0.23, 0.12, 'Socioeconomic\nstructure (MHDI income,\neducation, population)')
ac  = box(ax, 0.62, 0.77, 0.22, 0.12, 'Access to and\nquality of care (proxied\nin sensitivity models)', dashed=True)
rq  = box(ax, 0.87, 0.77, 0.22, 0.12, 'Registration quality\n(completeness,\nclassification)', dashed=True)
mc  = box(ax, 0.11, 0.55, 0.22, 0.12, 'Maternal composition\n(age, schooling,\nrace/skin colour)')
tm  = box(ax, 0.38, 0.40, 0.28, 0.08, 'True neonatal\nmortality', dashed=True, fc=GREY, bold=True)
rd  = box(ax, 0.52, 0.21, 0.24, 0.09, 'Registered\nneonatal deaths', fc=GREY)
rb  = box(ax, 0.82, 0.21, 0.24, 0.09, 'Registered\nlive births', fc=GREY)
for t in (un, st, ac, rq): arrow(ax, ctx, t)
arrow(ax, st, mc); arrow(ax, st, tm); arrow(ax, st, ac)
arrow(ax, mc, tm); arrow(ax, un, tm); arrow(ax, ac, tm)
arrow(ax, tm, rd); arrow(ax, rq, rd); arrow(ax, rq, rb)
ax.text(0.50, 0.055,
        'Solid white: measured covariate groups, or the partly observed context node (geographical neighbourhood defines the adjacency graph; historical context is not measured). '
        'Dashed white:\nconstructs not directly measured in the primary model. Grey: outcome processes (dashed, latent; solid, registered counts). '
        'The diagram simplifies the birth and death\nobservation processes and does not identify a sufficient adjustment set for causal effects.',
        ha='center', va='center', fontsize=11, linespacing=1.3)
fig.savefig(str(OUTFIG/'figure2_dag.png'), dpi=300, bbox_inches='tight'); fig.savefig(str(OUTFIG/'figure2_dag.pdf'), bbox_inches='tight'); plt.close(fig)
