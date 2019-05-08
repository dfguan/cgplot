import argparse, sys,os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from wiggle import Reader 

def get_abspath(p):
    return os.path.abspath(p)

def cgplot_func2(hits, covs, chrg_t, qrybd, fn, tit):
    [chrn, xmin, xmax] = chrg_t
    title = tit 
    # title2 = "t2" 
           
    height = 10
    width = 10
    dpi = 300
    

    colors = ["#000000", "#ef2928", "#ad7fa8", "#8ae234", "#729fcf", "#f2c27e", "#fcaf3e", "#fce94f"]
    t = 0
    color_idx = {}
    for k in covs:
        color_idx[k] = t
        t += 1

    maxcov = 0
    for k in covs:
        for z in covs[k]:
            if z[2] > maxcov:
                maxcov = z[2]
    
    # xmin = 0
    # xmax = 0
    # xlabel = []
    # xpos = []
    # start = 0
    # for k, q in sorted(qrybd.items(), key=lambda x:(x[1][2] + x[1][3]) / 2): 
        # xmax += q[1] - q[0] 
        # start -= q[0]
        # xpos.append(start + q[0])
        # xlabel.append("{0:.2f}".format(q[0]/1000000))
        # xpos.append(start + int((q[0] + q[1])/2))
        # xlabel.append("{}".format(k))
        # xpos.append(start + q[1] - 1000)
        # xlabel.append("{0:.2f}".format(q[1]/1000000))
        # start += q[1] - q[0]

    ymin2 = 0 
    ymax2 = int ((maxcov + 9) / 10 * 10) 
    

    fig = plt.figure(num=None, figsize=(width, height)) 
    
    rect1 = [0.07, 0.3, 0.90, 0.65]
    ax1 = fig.add_axes(rect1)
    ax1.set_ylim(xmin, xmax)
    ax1.set_xlim(xmin, xmax)
    # ax1.set_xticks(xpos)
    # ax1.set_xticklabels(xlabel, fontsize='xx-small')
    # ax1.tick_params(axis='x', direction='inout')
    # ax1.set_xticks(xpos[4:6])
    # ax1.set_xticklabels(xlabel[4:6], fontsize='xx-small')
    # ax1.tick_params(axis='x', direction='out')
    # i = 0
    # for txt in ax1.get_xticklabels():
        # if int(i / 3) % 2 == 1:
            # txt.set_y(0.04) 
        # i += 1
    # plt.setp(ax1.get_xticklabels(), visible=True)
    # ax1.set_yticks([])
    start = xmin 
    ttl = 0
    for k, q in sorted(qrybd.items(), key=lambda x:(x[1][2] + x[1][3]) / 2): 
        start = ttl - q[0]
        ttl += q[1] - q[0]
        for t in hits[k]:
            ax1.plot([t[2], t[3]], [start + t[0], start + t[1]], color=colors[color_idx[k]], label=k)  
    # ax1.set_title(title)
    ax1.set_xlabel("Contig (Mb)")
    ax1.set_ylabel("{} (Mb)".format(chrn))
    # zip legend
    handle, labels = ax1.get_legend_handles_labels()
    nhandles = []
    nlabels = []
    for h, l in zip(handle, labels):
        if l not in nlabels:
            nhandles.append(h)
            nlabels.append(l)
    ax1.legend(nhandles, nlabels, loc=0)
     
    rect2 = [0.07, 0.02, 0.90, 0.23]
    ax2 = fig.add_axes(rect2 ) 
    ax2.set_ylim([ymin2, ymax2])
    ax2.invert_yaxis()
    ax2.set_xlim([xmin, xmax])
    ax2.set_xticks([])
    # start = xmin
    for k, q in sorted(qrybd.items(), key=lambda x:(x[1][2] + x[1][3]) / 2): 
        start = q[2]
        for t in covs[k]:
            # print (t)
            ax2.fill([start + t[0], start + t[1], start + t[1], start + t[0]],[0, 0, t[2], t[2]],facecolor = colors[color_idx[k]], label=k)
        # start += q[3] - q[2]
    # for k in covs:
        # _st = qrybd[k][2]
    # ax2.set_title("pp")
    ax2.set_ylabel("read-depth")
    ax2.legend().set_visible(False)
    # handle, labels = ax2.get_legend_handles_labels()
    # nhandles = []
    # nlabels = []
    # for h, l in zip(handle, labels):
        # if l not in nlabels:
            # nhandles.append(h)
            # nlabels.append(l)
    # ax2.legend(nhandles, nlabels, loc=1)
    # ax2.legend(loc=1)
    # fig.tight_layout()
    fig.savefig(fn, dpi = dpi)

def cgplot_func(hits, covs, chrg_t, qrybd, fn, tit):
    [chrn, ymin, ymax] = chrg_t
    title = tit 
    title2 = "t2" 
           
    height = 10
    width = 10
    dpi = 300
    

    colors = ["#000000", "#ef2928", "#ad7fa8", "#8ae234", "#729fcf", "#f2c27e", "#fcaf3e", "#fce94f"]
    t = 0
    color_idx = {}
    for k in covs:
        color_idx[k] = t
        t += 1

    totcov = 0
    totlen = 0 
    for k in covs:
        for z in covs[k]:
            totcov += (z[1] - z[0] + 1) * z[2]
            totlen += (z[1] - z[0] + 1)
    avgcov = int ( totcov / totlen) 
    xmin = 0
    xmax = 0
    xlabel = []
    xpos = []
    start = 0
    beg = True
    for k, q in sorted(qrybd.items(), key=lambda x:(x[1][2] + x[1][3]) / 2): 
        start = xmax - q[0]
        xmax += q[1] - q[0] 
        # print ("{0} {1} {2}".format(k, q[0], q[1]))
        if not beg:
            xlabel[-1] += " | {0:.2f}".format(q[0]/1000000) 
        else:
            xpos.append(start + q[0])
            xlabel.append("{0:.2f}".format(q[0]/1000000))
            beg = False
        xpos.append(start + int((q[0] + q[1])/2))
        xlabel.append("{}".format(k))
        xpos.append(start + q[1] )
        xlabel.append("{0:.2f}".format(q[1]/1000000))

    # print (xpos)
    ymin2 = 0 
    ymax2 = int (2.5 * avgcov) 
    

    fig = plt.figure(num=None, figsize=(width, height)) 
    
    rect1 = [0.07, 0.3, 0.90, 0.65]
    ax1 = fig.add_axes(rect1)
    ax1.set_ylim(ymin/1000000, ymax/1000000)
    ax1.set_xlim(xmin, xmax)
    ax1.set_xticks(xpos)
    ax1.set_xticklabels(xlabel, fontsize='xx-small')
    # ax1.tick_params(axis='x', direction='inout')
    # ax1.set_xticks(xpos[4:6])
    # ax1.set_xticklabels(xlabel[4:6], fontsize='xx-small')
    # ax1.tick_params(axis='x', direction='out')
    # i = 0
    # for txt in ax1.get_xticklabels():
        # if int(i / 3) % 2 == 1:
            # txt.set_y(0.04) 
        # i += 1
    # plt.setp(ax1.get_xticklabels(), visible=True)
    # ax1.set_yticks([])
   
    for k in hits:
        for t in hits[k]:
            ax1.plot([t[0], t[1]], [t[2]/1000000, t[3]/1000000], color=colors[color_idx[k]], label=k)  
    # ax1.set_title(title)
    ax1.set_xlabel("Contig (Mb)")
    ax1.set_ylabel("{} (Mb)".format(chrn))
    # zip legend
    handle, labels = ax1.get_legend_handles_labels()
    nhandles = []
    nlabels = []
    for h, l in zip(handle, labels):
        if l not in nlabels:
            nhandles.append(h)
            nlabels.append(l)
    ax1.legend(nhandles, nlabels, loc=0)
     
    rect2 = [0.07, 0.02, 0.90, 0.23]
    ax2 = fig.add_axes(rect2 ) 
    ax2.set_ylim([ymin2, ymax2])
    ax2.invert_yaxis()
    ax2.set_xlim([xmin, xmax])
    ax2.set_xticks([])
    for k in covs:
        for t in covs[k]:
            # print (t)
            ax2.fill([t[0], t[1], t[1], t[0]],[0, 0, t[2], t[2]],facecolor = colors[color_idx[k]], label=k)
    # ax2.set_title("pp")
    ax2.set_ylabel("read-depth")
    ax2.legend().set_visible(False)
    # handle, labels = ax2.get_legend_handles_labels()
    # nhandles = []
    # nlabels = []
    # for h, l in zip(handle, labels):
        # if l not in nlabels:
            # nhandles.append(h)
            # nlabels.append(l)
    # ax2.legend(nhandles, nlabels, loc=1)
    # ax2.legend(loc=1)
    # fig.tight_layout()
    fig.savefig(fn, dpi = dpi)

def chk_fl(p):
    return os.path.isfile(p)

def update_coords(hits, covs, qrybd):
    start = 0 
    # for t in qrybd.items():
        # print (t)
    ttl = 0
    for k, q in sorted(qrybd.items(), key=lambda x:(x[1][2] + x[1][3]) / 2): 
        qe = q[1]
        qs = q[0]
        # print (q)
        start = ttl - qs
        ttl += qe - qs
        for i in range(len(hits[k])):
            
            hits[k][i][0] += start 
            hits[k][i][1] += start 
        for i in range(len(covs[k])):
            covs[k][i][0] += start
            covs[k][i][1] += start
def parse_chrg(chrg):
    chrglst = chrg.strip().split(':')

    chrn = chrglst[0]
    st = 0
    ed = 1000000000
    if len(chrglst) > 1:
        chrn = chrglst[0]
        cordlst = chrglst[1].split('-')
        if len(cordlst) < 2:
            print ("chromsome region format error")
            return [] 
        st = int(cordlst[0])
        ed = int(cordlst[1])
    return [chrn, st, ed]



def get_select_hits(paf_fn, chrg_t, qnslist, ml, qrybd):
    cnt = {}
    hitsd = {}
    sel_hitsd = {}
    for k in qnslist:
        cnt[k] = [0, 0]
        sel_hitsd[k] = []
        hitsd[k] = []
    [chrn, st, ed] = chrg_t
    with open(paf_fn, 'r') as f:
        for ln in f:
            hit = ln.split('\t')
            qn = hit[0]
            qs = int(hit[2])
            qe = int(hit[3])
            _chrn  = hit[5]
            _st = int(hit[7])
            _ed = int(hit[8])
            if qn in qnslist and  _chrn == chrn and _st < ed and _ed > st and qe - qs >= ml :
                if hit[4] == '-': 
                    cnt[qn][1] += qe - qs 
                    # sel_hitsd[qn].append([qs, qe, _ed, _st])
                else:
                    cnt[qn][0] += qe - qs 
                    # sel_hitsd[qn].append([qs, qe, _st, _ed])
                
                hitsd[qn].append([qs, qe, hit[4], _st, _ed])
                # if qs < qrybd[qn][0]: 
                    # qrybd[qn][0] = qs
                # if qe > qrybd[qn][1]:
                    # qrybd[qn][1] = qe
        f.close()
    for qn in hitsd:
        sel = '-'
        if cnt[qn][0] > cnt[qn][1]:
            sel = '+'  
        elif cnt[qn][0] == cnt[qn][1]: 
           continue  
        for hit in hitsd[qn]:
            if sel == hit[2]:
                if hit[2] == '-':
                    sel_hitsd[qn].append([hit[0], hit[1], hit[4], hit[3]])
                else:
                    sel_hitsd[qn].append([hit[0], hit[1], hit[3], hit[4]])
                
                if hit[0] < qrybd[qn][0]: 
                    qrybd[qn][0] = hit[0]
                if hit[1] > qrybd[qn][1]:
                    qrybd[qn][1] = hit[1]
                if hit[3] < qrybd[qn][2]:
                    qrybd[qn][2] = hit[3]
                if hit[4] > qrybd[qn][3]:
                    qrybd[qn][3] = hit[4]

    return sel_hitsd

def get_select_coverage(wigfn, qrybd, c, n):   
    wigd = {}
    sel_wigd = {}
    for k in qrybd:
        wigd[k] = []
        sel_wigd[k] = [] 
    with open(wigfn) as f:
        rd = Reader(f)
        for chrn, st, ed, v in rd:
            if chrn in qrybd:
                if st <= qrybd[chrn][1] and ed >= qrybd[chrn][0]:
                    _st = max(st, qrybd[chrn][0])
                    _ed = min(ed, qrybd[chrn][1])
                    wigd[chrn].append([_st, _ed, v])
    f.close()
    

    for k in wigd:
        fst = True
        _st = 0
        _ed = 0
        length = 0
        tot = 0
        wigdl = len(wigd[k])
        # print (wigdl)
        for i in range(wigdl + 1):
            # _st = wigd[k][i][0]
            if i == wigdl or i % n == 0:
                if not fst:
                    sel_wigd[k].append([_st, _ed, tot / length])
                if i != wigdl:
                    _st = wigd[k][i][0]
                    _ed = wigd[k][i][1]
                    tot = (_ed - _st) * wigd[k][i][2] 
                    length = _ed - _st
                fst = False
            else:
                _ed = wigd[k][i][1]
                tot += (_ed - wigd[k][i][0]) * wigd[k][i][2] 
                length += _ed - wigd[k][i][0]
    return sel_wigd

def worker(opts):
    pafn = (opts.paf_file)
    wigfn =(opts.wig_file)
    # print (wigfn)
    if not chk_fl(pafn) or not chk_fl(wigfn):
        print ("paf file or wig file is not found")
        return 1
    chrg_t = parse_chrg(opts.chrg)
    if not len(chrg_t):
        return 1
    qnslst = opts.qns.split(',')
    qrybd = {}
    for qry in qnslst:
        qrybd[qry] = [10000000000, 0, 10000000000, 0]
    hitsd = get_select_hits(pafn, chrg_t, qnslst, opts.mmapl, qrybd) 
    for k in hitsd:
        if len(hitsd[k]) == 0:
            print ("Error: a contig has no hits, please check query or chromsome name")
            return 1
    cov_lim = [0]
    # print (qrybd)
    covs = get_select_coverage(wigfn, qrybd, cov_lim, 5)
    for k in covs:
        if len(covs[k]) == 0:
            print ("Error: a contig has no coverage")
            return 1

    # cgplot_func2(hitsd, covs, chrg_t, qrybd)
    update_coords(hitsd, covs, qrybd)
    
    cgplot_func(hitsd, covs, chrg_t, qrybd, opts.out, opts.title) 
    return 0 





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Genome Comparison plot')

    parser.add_argument('-c', '--chrom', type=str, action="store", dest = "chrg", help ='chromsome region in chr:start-end or chr format', required=True)
    parser.add_argument('-q', '--query', type=str, action = "store", dest = "qns", help = 'query name(s) that fall(s) into the chromsome region, add comma to join multiple query names, support a maximum of 5 contigs', required = True)
    parser.add_argument('-l', '--mmapl', type = int, action = "store", dest = "mmapl", help = 'minimum mapped length [4000]', default= 4000)
    parser.add_argument('-o', '--out', type = str, action = "store", dest = "out", help = 'output file [plot.png]', default= "plot.png")
    parser.add_argument('-t', '--title', type = str, action = "store", dest = "title", help = 'figure title [NULL]')
    parser.add_argument('--version', action='version', version='%(prog)s 0.0.0')
    parser.add_argument('paf_file', type=str, action="store", help = "a paf file")
    parser.add_argument('wig_file', type=str, action="store", help = "a wig file")
    opts = parser.parse_args()
    sys.exit(worker(opts))
