import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from scipy import stats
import numpy as np
import time
from tools import Table, Manifest

def get_color(x):
    if x < 0.5:
        return (1-2*x,1-2*x,1)
    else:
        return (0,0,1.5-x)
        
def get_text_color(x):
    if x>0.5:
        return 'white'
    else:
                return 'black'

class PvalMatrix:
    def __init__(self,manifest,pvals):
        groups = self.read_manifest(manifest)
        self.header,self.ylabels,self.group_indices,self.group_counts,self.group_props = self.read_pvals(pvals,groups)
        self.xlabels = []
        for label,index in self.group_indices.items():
            a=5
            if len(label) > a+1:
                new = label[:a]
                for i in range(a,len(label),a):
                    new = new+'\n'+label[i:i+a]
                label = new
            self.xlabels.append(f"{label}\n(n={len(index)})")

    def read_manifest(self,filename):
        with open(filename) as manifest:
            groups = {}
            for line in manifest:
                row = line.rstrip().split("\t")
                if row[1] == "Tumor":
                    continue
                groups[row[0]] = row[1]
        return groups
    
    def read_pvals(self,pvals,groups,threshold=0.05):
        with open(pvals) as tsv:
            header = tsv.readline().split('\t')
            group_indices = {}
            for i,name in enumerate(header):
                if name in groups:
                    if groups[name] not in group_indices:
                        group_indices[groups[name]] = [i]
                    else:
                        group_indices[groups[name]].append(i)
            group_counts = []
            group_props = []
            ylabels = []
            for line in tsv:
                row = line.rstrip().split('\t')
                if "Tumor" in row[0]:
                    continue
                counts = []
                props = []
                for index in group_indices.values():
                    counts.append(len([i for i in index if float(row[i])<threshold]))
                    props.append(counts[-1]/len(index))
                ylabels.append(row[0])
                group_counts.append(counts)
                group_props.append(props)
        return header[1:],ylabels[::-1],group_indices,group_counts[::-1],group_props[::-1]

    def plot_table(self,filename,color=get_color,text=get_text_color):
        pw = .6 * len(self.xlabels)
        ph = .32 * len(self.ylabels)
        w = pw + 2.5
        h = ph + 1
        plt.figure(figsize=(w,h))
        panel = plt.axes([1.75/w,.9/h,pw/w,ph/h])
        colorbar = plt.axes([(1.85+pw)/w,.9/h,.15/w,.5/h]) 
        panel.set_yticks(range(1,len(self.ylabels)+1),self.ylabels)
        panel.set_xticks(range(1,len(self.xlabels)+1),self.xlabels)   
        colorbar.tick_params(left=False,right=True,bottom=False,top=False,
                             labelleft=False,labelright=True,labelbottom=False,labeltop=False)
        colorbar.set_yticks([0,0.5,1],['0%','50%','100%'])
        for i in range(100):
            p = i/100
            rectangle = patches.Rectangle([0,p],1,.01,lw=0,fc=color(p))
            colorbar.add_patch(rectangle)
        for i in range(len(self.ylabels)):
            for j in range(len(self.xlabels)):
                p = self.group_props[i][j]
                rectangle = patches.Rectangle([j+.5,i+.5],1,1,lw=0.5,ec='black',fc=color(p))
                panel.add_patch(rectangle)
                panel.text(j+1,i+1,self.group_counts[i][j],color=text(p),va='center',ha='center')
        panel.set_xlim(0.5,len(self.xlabels)+.5)
        panel.set_ylim(0.5,len(self.ylabels)+.5)
        plt.savefig(f"{filename}_matchtable.png")

class ColorLoop:
    def __init__(self,colors=["darkblue","darkorange","green","gray","purple","black"]):
        self.colors = colors
        self.i = -1

    def next(self):
        self.i = (self.i + 1) % len(self.colors)
        return self.colors[self.i]


def read_betas(filename,interval_set):
    with open(filename) as tsv:
        header = tsv.readline().rstrip().split('\t')[1:]
        indices = {}
        for i,column in enumerate(header):
            stat,name = column.split("_",1)
            if stat == "median":
                indices[name] = [i]
            elif stat == "alpha":
                indices[name].append(i)
            elif stat == "beta":
                indices[name].append(i)
        betas = {}
        for line in tsv:
            interval,row = line.split("\t",1)
            if interval in interval_set:
                row = [float(x) for x in row.split('\t')]
                mabs = {}
                for name,index in indices.items():
                    mabs[name] = [row[i] for i in index]
                betas[interval] = mabs
        return betas


class ColorBox:

    pal = {'doughton_green' : ("#155b51", "#216f63", "#2d8277", "#3a9387", "#45a395"),
           'doughton_purple' : ("#c468b2","#af509c", "#803777", "#5d2155", "#45113f"),
           "fritsch" : ("#0f8d7b", "#8942bd", "#eadd17", "#1e1a1a") }
    

    mld = ["#0f8d7b", "#8942bd", "#eadd17", "#1e1a1a",  # momacolors:Fritsch
           "orange","blue","red","gray"]

    def __init__(self,colors=None):
        if colors:
            self.colors = colors
        else:
            self.colors = ColorBox.mld
        self.i = 0
        self.labels = {}
        
    def add_label(self,label,color=None):
        if color:
            self.labels[label] = color
        else:
            self.labels[label] = self.colors[self.i]
            self.i = (self.i + 1) % len(self.colors)

    def get_color(self,label,default="darkgray"):
        if label in self.labels:
            return self.labels[label]
        else:
            return default

    def get_light(self,label,default="lightgray"):
        if label in self.labels:
            color = self.labels[label]
            return [x+(1-x)/2 for x in mcolors.to_rgb(color)]
        else:
            return default
        
    def get_dark(self,label,default="black"):
        if label in self.labels:
            color = self.labels[label]
            return [x*0.75 for x in mcolors.to_rgb(color)]
        else:
            return default
class PS_distribution:
    def __init__(self,interval,row,group_indices=None,betas={},width=0.05):
        self.interval = interval
        self.ymax = 0
        self.width = width
        self.bins = np.arange(0,1+width,width)
        self.bars = [[] for bin in self.bins[:-1]]
        self.colors = ColorBox()
        if group_indices:
            self.group_indices = group_indices
            for group in group_indices:
                self.colors.add_label(group)
        else:
            self.group_indices = {"Values":[i for i in range(len(row))]}
            self.colors.add_label(values)
        fw,fh = 6,3
        pw,ph = 4,2
        self.pw = pw
        self.ph = ph
        self.fig = plt.figure(figsize=(fw,fh))
        self.panel = self.fig.add_axes([0.5/fw,0.5/fh,pw/fw,ph/fh])
        self.hist_labels = []
        self.beta_labels = []
        for group,indices in self.group_indices.items():
            values = [row[i] for i in indices]
            self.add_hist(values,label=group)
        self.plot_hists()
        for name,mab in betas.items():
            m,a,b = mab
            self.add_beta(a,b,label=name)
        legend_size = len(self.hist_labels) + len(self.beta_labels)
        lw = 1 + (legend_size//6)
        self.legend = self.fig.add_axes([4.6/fw,0.5/fh,lw/fw,ph/fh])


    def add_hist(self,values,label,density=True):
        counts,bins = np.histogram(values,bins=self.bins,density=density)
        for i,count in enumerate(counts):
            self.bars[i].append((count,label))
        self.hist_labels.append(label)
        self.ymax = max(self.ymax,max(counts))

    def plot_hists(self):
        thick = (self.pw/self.ph) * (0.005*(1.1*self.ymax))
        for i,stack in enumerate(self.bars):
            left,right = self.bins[i],self.bins[i+1]
            alpha = 1
            for count,label in sorted(stack,reverse=True):
                r = patches.Rectangle((self.bins[i],0),self.width,count,
                #r = patches.Rectangle((left,0),right-left,count,
                                      linewidth=0.08,edgecolor="black",
                                      facecolor=self.colors.get_light(label),
                                      alpha=alpha,zorder=1)
                self.panel.add_patch(r)
                alpha = 0.5
                top_edge = patches.Rectangle((left,count),right-left,thick,
                                      linewidth=0.08,edgecolor=self.colors.get_dark(label),
                                      facecolor=self.colors.get_color(label),
                                      alpha=1,zorder=3)
                self.panel.add_patch(top_edge)
                left = left+0.005
                right = right-0.005
        ymax = self.ymax*1.1        
        self.panel.set_ylim(0,ymax)
        self.panel.set_xlim(0,1)
        self.panel.add_patch(patches.Rectangle((0,0),1,ymax,fill=None,lw=1,edgecolor="black",zorder=999))

    def add_beta(self,a,b,label):
        xs,ys = self.beta_points(a,b)
        self.panel.plot(xs,ys,color=self.colors.get_dark(label))
        self.beta_labels.append(label)

    def beta_points(self,a,b,xdist=0.005):
        xs,ys = [],[]
        for x in np.arange(xdist/2,1,xdist):
            xs.append(x)
            ys.append(stats.beta.pdf(x,a,b))
        return xs,ys
        
    def fill_legend(self):
        y = -1
        x = -3
        for label in self.hist_labels:
            y = (y+1) % 6
            if y == 0:
                x += 3
            r = patches.Rectangle((x+.1,y+.2),.3,.6,edgecolor="black",
                                  linewidth=.1,facecolor=self.colors.get_light(label))
            self.legend.add_patch(r)
            r = patches.Rectangle((x+.1,y+.2),.3,.1,edgecolor=self.colors.get_dark(label),
                                  linewidth=.1,facecolor=self.colors.get_color(label))
            self.legend.add_patch(r)
            self.legend.text(x+.5,y+.5,f"{label} PS values",ha='left',va='center')
        for label in self.beta_labels:
            y = (y+1) % 6
            if y == 0:
                x += 3
            self.legend.plot([x+.1,x+.4],[y+.5,y+.5],color=self.colors.get_dark(label))
            self.legend.text(x+.5,y+.5,f"{label} beta dist.",ha='left',va='center')
        self.legend.set_xticks([])
        self.legend.set_yticks([])  
        self.legend.set_xlim(0,x+3)
        self.legend.set_ylim(6,0)
           
    def save_fig(self,out_prefix,dpi=600):
        self.fill_legend()
        self.fig.savefig(f"{out_prefix}_{self.interval}.{time.time()}.png",dpi=dpi,bbox_inches="tight") 
        
class PCA_plot:
    def __init__(self,xs,ys,xy_pairs=[]):
        fw,fh = 6,4
        pw,ph = 4,2
        self.pw = pw
        self.ph = ph
        self.fig = plt.figure(figsize=(fw,fh))
        self.panel = self.fig.add_axes([0.5/fw,0.5/fh,pw/fw,ph/fh])
        self.hist_labels = []

    def save_fig(self,out_prefix,dpi=600):
        self.fill_legend()
        self.fig.savefig(f"{out_prefix}_pca.{time.time()}.png",dpi=dpi,bbox_inches="tight") 

def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--intervals",default=None,
                        help="List of intervals to plot (comma-separated).")
    parser.add_argument("-q","--query",default=None,
                        help="Table of p-values from signature query (pvals.tsv).")
    parser.add_argument("-m","--manifest",default=None,
                        help="Manifest of samples with groups for combining or labeling.")
    parser.add_argument("-p","--ps_table",default=None,
                        help="ps.tsv file with PS values for plotting.")
    parser.add_argument("-b","--beta",default=None,
                        help="beta.tsv file with parameters for fit beta distributions.")
    parser.add_argument("-o","--out_prefix",default="splicedice",
                        help="Output path and filename before extensions [Default: 'splicedice']")
    return parser.parse_args()




def main():
    args = get_args()
    if args.query and args.manifest:
        pmat = PvalMatrix(args.manifest,args.query)
        pmat.plot_table(args.out_prefix)

    if args.ps_table and args.intervals and args.beta and args.manifest:
        manifest = Manifest(args.manifest)
        intervals = set(args.intervals.split(","))
        betas = read_betas(args.beta,intervals)
        ps_table = Table(args.ps_table)
        group_indices = manifest.get_group_indices(ps_table.get_samples())
        for interval,row in ps_table.get_rows(interval_set=intervals):
            ps_plot = PS_distribution(interval,row,group_indices,betas[interval])
            ps_plot.save_fig(args.out_prefix)




if __name__=="__main__":
    main()

