__version__ = '1.0.1'

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

def samplemapToTissueCount(samplemap):
    tissues = {}
    for sample in samplemap:
        tissue = samplemap[sample]
        if tissue not in tissues:
            tissues[tissue] = 1
            continue
        tissues[tissue] += 1
    return tissues

def centroid(_list):
    return float(sum(_list))/len(_list)

# List of lists must be sorted in ascending order
def getListIndexMinMetric(_list, linkage='centroid', metric='euclidian'):
    if metric == 'euclidian':
        index = None
        _min = float('inf')
        for i in xrange(len(_list) - 1):
            if linkage == 'centroid':
                dist = centroid([x.end for x in _list[i+1]]) - centroid([x.end for x in _list[i]])
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Given KLEAT files identify tissue-specific cleavage sites as well as differentially utilized sites')
    parser.add_argument('kleat_bed', help='A tabix-indexed BED file containing all KLEAT calls you wish to analyze from all tissues. The name column should be the same as the first column in the samplemap argument')
    parser.add_argument('samplemap', help='A file mapping sample names to tissue site')
    parser.add_argument('-u', '--utr3s', default='/home/dmacmillan/annotations/ensembl/ensembl.fixed.sorted.utr3s_only.genes.sorted.gz', help='Bgzipped tabix-indexed BED file containing all viable 3\'UTRs to use')
    parser.add_argument('-n', '--normalize', action='store_true', help='Normalize frequencies by number of samples in tissue')
    parser.add_argument('-e', '--exclude', default=[], nargs='+', help='Tissue types to exclude')
    parser.add_argument("-l", "--log", dest="logLevel", default='WARNING', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the logging level. Default = "WARNING"')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))
    
    args = parser.parse_args()

    utr3s = pysam.tabix_iterator(open(args.utr3s), parser=pysam.asBed())

    kleats = pysam.TabixFile(args.kleat_bed, parser=pysam.asBed())

    if not os.path.isdir(args.outdir):
        try:
            os.makedirs(args.outdir)
        except OSError:
            pass

    # Logging
    logging.basicConfig(filename=os.path.join(args.outdir, 'log_{}'.format(os.path.basename(__file__))), level=getattr(logging, args.logLevel), filemode='w')

    mapping = mapSamples(args.samplemap)
    logging.debug('mapping: {}'.format(mapping))

    if args.normalize:
        normalize = samplemapToTissueCount(mapping)
    
    for utr3 in utr3s:
        logging.debug('utr3: {}'.format(utr3))
        gene = utr3.name
        logging.debug('gene: {}'.format(gene))
        lutr3 = utr3.end - utr3.start
        calls = [x for x in kleats.fetch(utr3.contig, utr3.start-20, utr3.end+20)]
        tissues = {mapping[x.name]: None for x in calls}
        tkeys = tissues.keys()
        ltissues = len(tissues)
        logging.debug('calls: {}'.format(calls))
        clusters = iterAHC([x for x in calls])
        lclusters = len(clusters)
        mtx = Matrix(ltissues, lclusters)
        centroids = []
        for cluster in clusters:
            values = [x.end for x in cluster]
            centroids.append(centroid(values))
        for i in xrange(mtx.m):
            tissue = tkeys[i]
            for j in xrange(mtx.n):
                values = [x.end for x in clusters[j] if mapping[x.name] == tissue]
                if not values:
                    continue
                cs = centroid(values)
                score = len(values)
                if args.normalize:
                    score /= normalize[tissue]
                mtx.mtx[i][j] = score
        mtx.rownames = [t.replace(' ', '_') for t in tissues.keys()]
        if utr3.strand == '+':
            mtx.colnames = ['c_{}'.format(int(100*(x - utr3.start)/(lutr3))) for x in centroids]
        else:
            mtx.colnames = ['c_{}'.format(int(100*(utr3.end - x)/(lutr3))) for x in centroids]
        for i in xrange(mtx.n):
            col = mtx.col(i)
            if len(col) - col.count(0) == 1:
                gname = '{}_{}_{}'.format(gene, utr3.start, utr3.end)
                with open(os.path.join(args.outdir, gname), 'w') as f:
                    f.write(str(mtx))
                path = os.path.abspath(args.outdir)
                if mtx.n == 1:
                    rscript = genRHeatmap(gname, path, 'FALSE')
                else:
                    rscript = genRHeatmap(gname, path)
