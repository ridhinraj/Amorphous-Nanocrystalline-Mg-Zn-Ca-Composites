"""Shared figure style for the Mg-Zn-Ca AN-composite manuscript.  One palette, fixed meanings:
   Zn content / Zn atoms ........ blue        Ca content / Ca atoms ... vermilion
   nanocrystal volume fraction .. green        seed radius ............. reddish purple
   temperature .................. amber        strain rate ............. grey
   pure-glass control / reference  black       (Okabe-Ito based, colour-blind safe; lightness steps carry the ordering)"""
import os
import matplotlib as mpl, matplotlib.pyplot as plt
ZN  = {"light":"#9CC7E4","mid":"#0072B2","dark":"#003F66"}
CA  = {"light":"#F2B38A","mid":"#D55E00","dark":"#7F3700"}
VF5 = ["#66C2A0","#2CA25F","#009E73","#00755A","#004D3A"]       # 5, 10, 15, 18.5, 18.9 %
SEED= {"12":"#D98CB9","15":"#C2187F","20":"#7A2159"}
TEMP= {250:"#F0D060",300:"#E69F00",350:"#B87700",400:"#6E4600"}
RATE= {"1e8":"#BDBDBD","5e8":"#8C8C8C","1e9":"#4D4D4D","2e9":"#111111"}
GLASS="#000000"; MG="#7F7F7F"; GRID="#E6E6E6"; TXT="#222222"
CA_LEVEL={0.03:CA["light"],0.05:CA["mid"],0.07:CA["dark"]}
def apply():
    mpl.rcParams.update({"font.family":"Liberation Sans","font.size":8,"axes.labelsize":8.5,"axes.titlesize":9,
      "axes.titleweight":"bold","axes.titlelocation":"left","xtick.labelsize":7.5,"ytick.labelsize":7.5,"legend.fontsize":7,
      "legend.frameon":False,"axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.8,"xtick.major.width":0.8,
      "ytick.major.width":0.8,"xtick.major.size":3,"ytick.major.size":3,"lines.linewidth":1.4,"lines.markersize":5,
      "figure.dpi":150,"savefig.dpi":600,"pdf.fonttype":42,"ps.fonttype":42,"axes.edgecolor":TXT,"text.color":TXT,
      "axes.labelcolor":TXT,"xtick.color":TXT,"ytick.color":TXT,"mathtext.fontset":"custom","mathtext.rm":"Liberation Sans",
      "mathtext.it":"Liberation Sans:italic","mathtext.bf":"Liberation Sans:bold"})
def save(fig,name):
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rendered")
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, name + ".png"), bbox_inches="tight", pad_inches=0.03)
    fig.savefig(os.path.join(out, name + ".pdf"), bbox_inches="tight", pad_inches=0.03)
def panel(ax,letter,title): ax.set_title(f"({letter}) {title}",fontsize=8.5,fontweight="bold",loc="left",pad=4)
