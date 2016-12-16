__version__ = '1.0.0'

import argparse
import os
import sys
import logging
import time
import random
import colorsys
import subprocess
import pysam
from Kleat import *
from Matrix import *

def genRHeatmap(gene, path, clust='TRUE'):
    out = os.path.join(path, 'plot_{}.r'.format(gene))
    with open(out, 'w') as f:
        f.write(
            'setwd("{}")\n' \
            'library(pheatmap)\n' \
            'library(RColorBrewer)\n' \
            'mtx = read.table("{}")\n' \
            'mybreaks = seq(min(mtx), max(mtx), by = 1)\n' \
            'colors = rev(brewer.pal(max(min(length(mybreaks),11),3), "RdYlBu"))\n' \
            'png("{}.png", width=1000, height=1000, res=150)\n' \
            'pheatmap(mtx, main="{}", legend_breaks = mybreaks, cluster_rows = {}, cluster_cols = {})\n' \
            'dev.off()'.format(path, gene, gene, gene, clust, clust)
        )
    return out

def mapSamples(samplemap):
    mapping = {}
    with open(samplemap, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            mapping[line[0]] = line[1]
    return mapping

def centroid(_list):
    return float(sum(_list))/len(_list)

# List of lists must be sorted in ascending order
def getListIndexMinMetric(_list, linkage='centroid', metric='euclidian'):
    if metric == 'euclidian':
        index = None
        _min = float('inf')
        for i in xrange(len(_list) - 1):
            if linkage == 'centroid':
                dist = centroid([x.cleavage_site for x in _list[i+1]]) - centroid([x.cleavage_site for x in _list[i]])
            if dist < _min:
                index = i
                _min = dist
        return index, _min

# Given a sorted list, a linkage criteria, 
# and a threshold, return a
# list of lists where each internal list
# represents a cluster
def iterAHC(_list, linkage='centroid', threshold=20):
    clusters = [[x] for x in _list]
    index, _min = getListIndexMinMetric(clusters)
    while (_min <= threshold) and (len(clusters) > 1):
        clusters[index] += clusters[index+1]
        del(clusters[index+1])
        index, _min = getListIndexMinMetric(clusters)
    return clusters

#def outputClusters(clusters, tissues):
#    lclusters = len(clusters)
#    for i,cluster in enumerate(clusters):
#        values = [x.cleavage_site for x in cluster]
#        centroid = centroid(values)
#        tissue = tissue.keys()[i]

def groupUtr3s(utr3s):
    d = {}
    for utr in utr3s:
        if utr.name not in d:
            d[utr.name] = [utr]
            continue
        d[utr.name].append(utr)
    return d

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Given KLEAT files identify tissue-specific cleavage sites as well as differentially utilized sites')
    parser.add_argument('kleats', nargs='+', help='One or more KLEAT files to analyze')
    parser.add_argument('-u', '--utr3s', default='/home/dmacmillan/annotations/ensembl/ensembl.fixed.sorted.utr3s_only.genes.sorted.gz', help='Bgzipped tabix-indexed BED file containing all viable 3\'UTRs to use')
    parser.add_argument('-e', '--exclude', default=[], nargs='+', help='Tissue types to exclude')
    parser.add_argument('-s', '--samplemap', default=os.path.dirname(os.path.realpath(__file__))+'/samplemap', help='A file mapping sample names to tissue site')
    parser.add_argument('-n', '--normalize', action='store_true', help='Enable this to randomly subsample all tissues to equal the minimum number of samples')
    parser.add_argument("-l", "--log", dest="logLevel", default='WARNING', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the logging level. Default = "WARNING"')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))
    
    args = parser.parse_args()
    
    utr3s = pysam.tabix_iterator(open(args.utr3s), parser=pysam.asBed())
    gutr3s = groupUtr3s(utr3s)

    if not os.path.isdir(args.outdir):
        try:
            os.makedirs(args.outdir)
        except OSError:
            pass

    # Logging
    logging.basicConfig(filename=os.path.join(args.outdir, 'log_{}'.format(os.path.basename(__file__))), level=getattr(logging, args.logLevel), filemode='w')

    kleats = []
    mapping = mapSamples(args.samplemap)
    logging.debug('mapping: {}'.format(mapping))
    #logging.debug('mapping: {}'.format(mapping))
    tissues = {}
    #tissues = {x[1]:{'file': open(os.path.join(args.outdir, '{}.bg'.format(x[1])), 'w'), 'samples': set()} for x in mapping.items()}
    #for tissue in args.exclude:
    #    del(tissues[tissue])
    n = len(args.kleats)
    for i,k in enumerate(args.kleats):
        sample = os.path.basename(k).split('.')[0]
        logging.debug('sample: {}'.format(sample))
        tissue = mapping[sample]
        if tissue in args.exclude:
            continue
        if tissue not in tissues:
            tissues[tissue] = {
                'file': None,#open(os.path.join(args.outdir, '{}.bg'.format(tissue)), 'w'),
                'samples': set()
            }
        tissues[tissue]['samples'].add(sample)
        kleats += Kleat.parseKleatFast(k, sample)
        sys.stdout.write('Processed {}/{}\r'.format(i+1,n))
        sys.stdout.flush()
    print

    logging.debug('tissues: {}'.format(tissues))
    tkeys = tissues.keys()
    N = len(tissues)
    hsv = [(x*1.0/N, 1, 1) for x in range(N)]
    rgb = map(lambda x: colorsys.hsv_to_rgb(*x), hsv)
    #for i,tissue in enumerate(tissues):
    #    tissues[tissue]['file'].write('track type=\'bedGraph\' visibility=\'2\' name=\'{}\' color=\'{}\'\n'.format(tissue, (',').join([str(x*255) for x in rgb[i]])))
    
    kleats = Kleat.groupKleat(kleats)

    for chrom in kleats:
        for gene in kleats[chrom]:
            kleats[chrom][gene] = iterAHC(kleats[chrom][gene])
            clusters = kleats[chrom][gene]
            lclusters = len(clusters)
            mtx = Matrix(N, lclusters)
            centroids = []
            myutr3 = None
            for cluster in clusters:
                values = [x.cleavage_site for x in cluster]
                centroids.append(centroid(values))
            try:
                utrs = gutr3s[gene]
            except KeyError:
                continue
            if len(utrs) > 1:
                for u in utrs:
                    if all([(u.start-20 <= x <= u.end+20) for x in centroids]):
                        myutr3 = u
            else:
                myutr3 = utrs[0]
                if not all([(myutr3.start-20 <= x <= myutr3.end+20) for x in centroids]):
                    continue
            try:
                lmyutr3 = myutr3.end - myutr3.start
            except AttributeError:
                continue
            for i in xrange(mtx.m):
                tissue = tkeys[i]
                for j in xrange(mtx.n):
                    values = [x.cleavage_site for x in clusters[j] if mapping[x.name] == tissue]
                    if not values:
                        continue
                    cs = centroid(values)
                    score = len(values)
                    factor = len(tissues[tissue]['samples'])
                    if args.normalize:
                        score /= float(factor)
                    #tissues[tissue]['file'].write('{}\t{}\t{}\t{}\n'.format(chrom, cs-1, cs, score))
                    mtx.mtx[i][j] = score
            mtx.rownames = [t.replace(' ', '_') for t in tissues.keys()]
            #mtx.colnames = ['c{}_{}'.format(i, int(x)) for i,x in enumerate(centroids)]
            if myutr3.strand == '+':
                mtx.colnames = ['c_{}'.format(int(100*(x - myutr3.start)/(lmyutr3))) for x in centroids]
            else:
                mtx.colnames = ['c_{}'.format(int(100*(myutr3.end - x)/(lmyutr3))) for x in centroids]
            for i in xrange(mtx.n):
                col = mtx.col(i)
                if len(col) - col.count(0) == 1:
                    with open(os.path.join(args.outdir, gene), 'w') as f:
                        f.write(str(mtx))
                    path = os.path.abspath(args.outdir)
                    if mtx.n == 1:
                        rscript = genRHeatmap(gene, path, 'FALSE')
                    else:
                        rscript = genRHeatmap(gene, path)
                    #print rscript
                    #subprocess.Popen(['/gsc/btl/linuxbrew/bin/Rscript', rscript])
#            deltas = Matrix(N, lclusters-1)
#            for i in xrange(deltas.m):
#                for j in xrange(deltas.n):
#                    value = mtx.mtx[i][j+1] - mtx.mtx[i][j]
#                    deltas.mtx[i][j] = value
