"""Render manuscript cell atlases and annotation evidence from archived coordinates."""
from figure_support import *
from PIL import Image, ImageDraw
from reportlab.lib.utils import ImageReader
import io,base64,colorsys

INP=ROOT/'revision_inputs'
ANN=INP/'annotation_sensitivity'
emb={t:pd.read_csv(INP/'embeddings'/(t+'_embedding.tsv.gz'),sep='\t') for t in ['SGN','Lateral_wall']}
sc=pd.read_csv(ANN/'figure_inputs/CELL_TSNE_ANNOTATION_EXPRESSION_UCELL.tsv.gz',sep='\t')
MAP={'Endothelial_vascular':'Endothelial','Fibrocyte_mesenchymal':'Fibrocyte / mesenchymal','Immune_myeloid':'Immune / myeloid','Glial':'Glial / Schwann-like','Supporting_epithelial':'Epithelial / supporting-associated','Hair_sensory':'Hair / sensory','UNRESOLVED':'UNRESOLVED'}
PAL={**col,'Hair / sensory':'#46A5AE'}
SHORT={**short,'Hair / sensory':'Hair / sensory'}
sc['broad_annotation']=sc.original_broad.map(MAP)
CONDITION={'Naive':BLUE,'Noise':RED,'Control':BLUE,'HHL':'#C5A044','NIHL':RED}

def cloud(f,d,xy,field,palette,x,y,w,h,axes=('UMAP 1','UMAP 2'),cluster_labels=False,continuous=False):
    a=d[list(xy)].to_numpy(float)
    lo=a.min(axis=0);hi=a.max(axis=0);span=hi-lo
    lo-=span*.035;hi+=span*.035
    # Equal coordinate scaling preserves the original embedding geometry.
    scale=min(w/(hi[0]-lo[0]),h/(hi[1]-lo[1]))
    center=(hi+lo)/2
    px=(a[:,0]-center[0])*scale+w/2
    py=-(a[:,1]-center[1])*scale+h/2
    ss=600/72
    im=Image.new('RGB',(round(w*ss),round(h*ss)),'white'); dr=ImageDraw.Draw(im)
    vals=d[field].astype(str).to_numpy()
    if continuous:
        vv=d[field].to_numpy(float);vmax=max(float(vv.max()),1e-9)
        order=np.argsort(vv,kind='stable')
        r=.28*ss
        for i in order:
            c='#E3E7EA' if vv[i]==0 else seq(float(vv[i]),vmax)
            xx=px[i]*ss;yy=py[i]*ss;dr.ellipse((xx-r,yy-r,xx+r,yy+r),fill=c)
    else:
        groups=list(d[field].astype(str).value_counts().index)
        if 'Other cells' in groups:groups.remove('Other cells');groups.insert(0,'Other cells')
        if 'Cluster 19' in groups:groups.remove('Cluster 19');groups.append('Cluster 19')
        for g in groups:
            inds=np.flatnonzero(vals==g);r=(.47 if g=='Cluster 19' else .28)*ss
            c=palette[g]
            for i in inds:
                xx=px[i]*ss;yy=py[i]*ss;dr.ellipse((xx-r,yy-r,xx+r,yy+r),fill=c)
    bio=io.BytesIO();im.save(bio,format='PNG');raw=bio.getvalue()
    f.pdf.drawImage(ImageReader(io.BytesIO(raw)),x,f.h-y-h,width=w,height=h)
    f.svg.append(f'<image x="{x}" y="{y}" width="{w}" height="{h}" href="data:image/png;base64,{base64.b64encode(raw).decode()}"/>')
    ax=x+4;ay=y+h-4
    f.line(ax,ay,ax+25,ay,BLACK,.65);f.line(ax,ay,ax,ay-25,BLACK,.65)
    f.text(ax+14,ay+10,axes[0],6.8,align='center');f.text(ax-6,ay-13,axes[1],6.8,align='center',rotate=-90)
    if cluster_labels:
        rows=[]
        for g in sorted(set(vals),key=lambda z:int(z)):
            ii=vals==g; rows.append([g,float(np.median(px[ii])),float(np.median(py[ii]))])
        anchors=np.array([r[1:] for r in rows]);pos=anchors.copy()
        for _ in range(90):
            for i in range(len(pos)):
                for j in range(i):
                    delta=pos[i]-pos[j];dist=np.linalg.norm(delta)
                    if dist<10:
                        if dist<.01:delta=np.array([1.,.5]);dist=np.linalg.norm(delta)
                        step=(10-dist)*.12*delta/dist;pos[i]+=step;pos[j]-=step
            pos+=.015*(anchors-pos);pos[:,0]=np.clip(pos[:,0],6,w-6);pos[:,1]=np.clip(pos[:,1],6,h-6)
        for (g,_,_),p,a0 in zip(rows,pos,anchors):
            xx=x+p[0];yy=y+p[1]
            if np.linalg.norm(p-a0)>4:f.line(xx,yy,x+a0[0],y+a0[1],'#666666',.3)
            f.circle(xx,yy,4.8,'#FFFFFF',None);f.text(xx,yy+2,g,6.4,bold=True,align='center')
    return (lo,hi)

def caption_key(f):
    order=types+['Hair / sensory']
    legend(f,[(SHORT[t],PAL[t]) for t in order],18,508,cols=2,dx=250,dy=14,size=8)


def fig1():
    f=Fig('Figure1_cell_atlases',h=580)
    f.panel('A','Spiral ganglion region | 31,556 cells',8,16)
    cloud(f,emb['SGN'],['dim1','dim2'],'broad_annotation',PAL,25,35,215,172)
    f.panel('B','Lateral wall | 33,149 cells',263,16)
    cloud(f,emb['Lateral_wall'],['dim1','dim2'],'broad_annotation',PAL,280,35,215,172)
    f.panel('C','GSE252285 | 36,111 cells',8,255)
    cloud(f,sc,['tSNE_1','tSNE_2'],'broad_annotation',PAL,25,275,215,172,('t-SNE 1','t-SNE 2'))
    f.panel('D','Cellular composition',263,255)
    frames=[emb['SGN'],emb['Lateral_wall'],sc]
    counts=[];x0=298;y0=289;y1=439
    y_axis(f,x0,y0,y1,[0,50,100],0,100,'Cells (%)')
    for j,(name,d) in enumerate(zip(['SGN','Lateral wall','GSE252285'],frames)):
        x=x0+19+j*62;top=y1
        for ct in types+['Hair / sensory']:
            n=int((d.broad_annotation==ct).sum());height=n/len(d)*(y1-y0)
            if height:f.rect(x,top-height,31,height,PAL[ct]);top-=height
            counts.append({'dataset':name,'cell_type':ct,'cells':n,'fraction':n/len(d)})
        f.text(x+15.5,455,name,7.2,align='center')
    caption_key(f);save_table(pd.DataFrame(counts),'Figure1D_cell_type_composition.tsv');f.finish()


def fig6():
    f=Fig('Figure6_cluster19_expression',h=590)
    highlight=sc.copy();highlight['highlight']='Other cells'
    highlight.loc[highlight.original_broad=='Fibrocyte_mesenchymal','highlight']='Fibrocyte'
    highlight.loc[highlight.original_broad=='Endothelial_vascular','highlight']='Endothelial'
    highlight.loc[highlight.cluster==19,'highlight']='Cluster 19'
    pp={'Other cells':'#E2E6E9','Fibrocyte':PAL['Fibrocyte / mesenchymal'],'Endothelial':PAL['Endothelial'],'Cluster 19':'#8E4C99'}
    f.panel('A','Cluster 19 in the cochlear cell atlas',8,16)
    cloud(f,highlight,['tSNE_1','tSNE_2'],'highlight',pp,25,35,208,170,('t-SNE 1','t-SNE 2'))
    legend(f,[(g,pp[g]) for g in ['Endothelial','Fibrocyte','Cluster 19']],18,229,cols=3,dx=81,size=6.8)
    f.panel('B','Mesenchymal markers in cluster 19',263,16)
    genes=['Coch','Otor','Pou3f4','Col1a2','Pdgfra','Pecam1','Kdr','Vwf','Apold1']
    groups=[('Endothelial',sc.original_broad=='Endothelial_vascular'),('Fibrocyte',sc.original_broad=='Fibrocyte_mesenchymal'),('Cluster 19',sc.cluster==19)]
    rec=[]
    for name,mask in groups:
        for g in genes:
            v=sc.loc[mask,'expr__'+g];rec.append({'population':name,'gene':g,'n_cells':len(v),'mean_log_normalized':v.mean(),'pct_expressing':(v>0).mean()*100})
    dt=pd.DataFrame(rec);vmax=float(dt.mean_log_normalized.max())
    for i,(name,mask) in enumerate(groups):
        yy=60+i*43;f.text(323,yy+3,name,7.5,align='right');f.text(323,yy+13,'n = '+str(int(mask.sum())),6.5,GRAY,align='right')
        for j,g in enumerate(genes):
            r=dt[(dt.population==name)&(dt.gene==g)].iloc[0];xx=338+j*18
            f.circle(xx,yy,7*math.sqrt(r.pct_expressing/100),seq(r.mean_log_normalized,vmax))
    for j,g in enumerate(genes):f.text(338+j*18,173,g,7,italic=True,rotate=-55)
    f.text(285,218,'Mean expression',6.8)
    for i in range(70):f.rect(358+i,212,1.1,7,seq(vmax*i/69,vmax))
    f.text(356,231,'0',6.5);f.text(431,231,f'{vmax:.1f}',6.5)
    f.circle(459,216,3.5,GRAY);f.circle(485,216,7,GRAY);f.text(459,231,'25%',6.5,align='center');f.text(485,231,'100%',6.5,align='center')
    save_table(dt,'Figure6B_cluster19_marker_dotplot.tsv')
    summary=pd.read_csv(ANN/'TOP20_EXPRESSION_BOTH_ANNOTATIONS.tsv',sep='\t')
    c19=pd.read_csv(ANN/'CLUSTER19_MARKER_EXPRESSION.tsv',sep='\t')
    values=[]
    for condition in ['ALL','Control','HHL','NIHL']:
        for label,var,ct in [('Endothelial','original','Endothelial_vascular'),('Cluster 19',None,None),('Fibrocyte: original','original','Fibrocyte_mesenchymal'),('Fibrocyte: reassigned','c19_fibrocyte','Fibrocyte_mesenchymal')]:
            if var is None:r=c19[(c19.condition==condition)&(c19.gene=='Apold1')].iloc[0]
            else:r=summary[(summary.condition==condition)&(summary.gene=='Apold1')&(summary.variant==var)&(summary.broad_cell_type==ct)].iloc[0]
            values.append({'condition':condition,'population':label,'n_cells':int(r.n_cells),'mean_log_normalized':r.mean_log_normalized,'pct_expressing':r.pct_expressing})
    v=pd.DataFrame(values);save_table(v,'Figure6CD_Apold1_expression.tsv')
    fills=[pp['Endothelial'],pp['Cluster 19'],'#F2E8A9',pp['Fibrocyte']]
    for lab,title,key,yy,limit,ticks in [('C','Apold1 expression across populations','mean_log_normalized',280,max(2.,math.ceil(v.mean_log_normalized.max()*2)/2),None),('D','Apold1 detection across populations','pct_expressing',438,100,[0,50,100])]:
        f.panel(lab,title,8,yy-18);x0=59;y0=yy;y1=yy+90
        if ticks is None:ticks=list(np.arange(0,limit+.01,.5))
        y_axis(f,x0,y0,y1,ticks,0,limit,'Mean log-normalized' if key.startswith('mean') else 'Cells expressing (%)')
        for tick in ticks:
            y=y1-tick/limit*90;f.line(x0,y,493,y,'#E5E8EB',.35)
        for j,condition in enumerate(['ALL','Control','HHL','NIHL']):
            sub=v[v.condition==condition]
            for k,r in enumerate(sub.itertuples()):
                xx=79+j*109+k*18;val=getattr(r,key);h=val/limit*90
                f.rect(xx,y1-h,13,h,fills[k],pp['Fibrocyte'] if k==2 else None)
                f.text(xx+6.5,y1-h-4,f'{val:.2f}' if key.startswith('mean') else f'{val:.0f}',6.3,align='center')
            f.text(112+j*109,y1+15,'Pooled' if condition=='ALL' else condition,8,align='center')
    legend(f,list(zip(['Endothelial','Cluster 19','Fibrocytes','Fibrocytes + cluster 19'],fills)),18,574,cols=4,dx=125,size=7)
    f.finish()


def supp5():
    f=Fig('FigureS5_cluster_and_condition_maps',h=625)
    data=[('SGN',emb['SGN'],['dim1','dim2'],'leiden',('UMAP 1','UMAP 2')),('Lateral wall',emb['Lateral_wall'],['dim1','dim2'],'leiden',('UMAP 1','UMAP 2')),('GSE252285',sc,['tSNE_1','tSNE_2'],'cluster',('t-SNE 1','t-SNE 2'))]
    for row,(name,d,xy,cl,axes) in enumerate(data):
        y=16+row*202
        f.panel(chr(65+row*2),name+' | clusters',8,y)
        dd=d.copy();dd[cl]=dd[cl].astype(str)
        gs=sorted(dd[cl].unique(),key=int)
        pal={g:'#'+''.join(f'{round(c*255):02x}' for c in colorsys.hsv_to_rgb((i*.6180339)%1,.48,.77)) for i,g in enumerate(gs)}
        cloud(f,dd,xy,cl,pal,25,y+16,213,154,axes,True)
        f.panel(chr(66+row*2),name+' | condition',263,y)
        cloud(f,d,xy,'condition',CONDITION,280,y+16,213,154,axes)
        conds=['Naive','Noise'] if row<2 else ['Control','HHL','NIHL']
        legend(f,[(c,CONDITION[c]) for c in conds],305,y+183,cols=3,dx=67,size=7)
    f.finish()


def supp8():
    genes=['Apold1','Coch','Otor','Pecam1','Vwf','Errfi1','Nfkbia','Egr1']
    f=Fig('FigureS8_gene_feature_maps',h=720)
    stats=[]
    for i,g in enumerate(genes):
        coln=i%2;row=i//2;x=8+coln*255;y=16+row*176
        f.text(x,y,chr(65+i).lower(),11,bold=True);f.text(x+16,y,g,9,italic=True)
        cloud(f,sc,['tSNE_1','tSNE_2'],'expr__'+g,{},x+17,y+12,190,135,('t-SNE 1','t-SNE 2'),continuous=True)
        vmax=float(sc['expr__'+g].max())
        for k in range(65):f.rect(x+163+k*.8,y+155,.85,5,seq(vmax*k/64,vmax))
        f.text(x+163,y+169,'0',6);f.text(x+215,y+169,f'{vmax:.2f}',6,align='right')
        stats.append({'gene':g,'cells':len(sc),'scale_min':0,'scale_max':vmax,'n_expressing':int((sc['expr__'+g]>0).sum())})
    save_table(pd.DataFrame(stats),'FigureS8_feature_scales.tsv');f.finish()

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--only',choices=['fig1','fig6','supp5','supp8','all'],default='all');a=ap.parse_args()
    for n,fun in [('fig1',fig1),('fig6',fig6),('supp5',supp5),('supp8',supp8)]:
        if a.only in ['all',n]:fun()
