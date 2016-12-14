import argparse
import os
import sys
import logging
import time
import random
from Kleat import *

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
def iterAHC(_list, linkage='centroid', threshold=15):
    clusters = [[x] for x in _list]
    index, _min = getListIndexMinMetric(clusters)
    while (_min <= threshold) and (len(clusters) > 1):
        clusters[index] += clusters[index+1]
        del(clusters[index+1])
        index, _min = getListIndexMinMetric(clusters)
    return clusters

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Given KLEAT files identify tissue-specific cleavage sites as well as differentially utilized sites')
    parser.add_argument('-k1', nargs='+', help='One or more KLEAT files to analyze')
    parser.add_argument('-k2', nargs='+', help='One or more KLEAT files to compare to files in k1')
    parser.add_argument('-n', '--normalize', action='store_true', help='Enable this to randomly subsample k1 or k2 to both equal the minimum number of samples between the two')
    parser.add_argument('-a', '--names', nargs=2, default=('k1', 'k2'), help='Name k1 and k2 sets')
    parser.add_argument("-l", "--log", dest="logLevel", default='WARNING', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the logging level. Default = "WARNING"')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))
    
    args = parser.parse_args()
    if not args.k1 or not args.k2:
        sys.exit('-k1 and -k2 required')
    
    if not os.path.isdir(args.outdir):
        try:
            os.makedirs(args.outdir)
        except OSError:
            pass

    colours = ('255,0,0', '0,0,255')
    datasets = dict.fromkeys(args.names)
    keys = datasets.keys()
    for i,data in enumerate(datasets):
        datasets[data] = {'out': open(os.path.join(args.outdir, '{}.bg'.format(data)), 'w')}
        datasets[data]['out'].write('track type=\'bedGraph\' visibility=\'2\' name=\'{}\' color=\'{}\'\n'.format(data, colours[i]))
    
    # Logging
    logging.basicConfig(filename=os.path.join(args.outdir, 'log_{}'.format(os.path.basename(__file__))), level=getattr(logging, args.logLevel))

    if args.normalize:
        if len(args.k1) < len(args.k2):
            args.k2 = random.sample(args.k2, len(args.k1))
        else:
            args.k1 = random.sample(args.k1, len(args.k2))
    
    kleats = []
    for k in args.k1:
        kleats += Kleat.parseKleatFast(k, keys[0])
    for k in args.k2:
        kleats += Kleat.parseKleatFast(k, keys[1])

    kleats = Kleat.groupKleat(kleats)

    for chrom in kleats:
        for gene in kleats[chrom]:
            kleats[chrom][gene] = iterAHC(kleats[chrom][gene])
            diffs = []
            for cluster in kleats[chrom][gene]:
                cdiff = []
                for data in datasets:
                    cdiff.append(len([x for x in cluster if x.name == data]))
                diffs.append(cdiff)
            diffs = [x[0] - x[1] for x in diffs]
            _min = min(diffs)
            _max = max(diffs)
            if abs(_max - _min) > 10:
                print gene
                for data in datasets:
                    for i,cluster in enumerate(kleats[chrom][gene]):
                        css = [x.cleavage_site for x in cluster if x.name == data]
                        if not css:
                            continue
                        score = len(css)
                        cs = int(centroid(css))
                        datasets[data]['out'].write('{}\t{}\t{}\t{}\n'.format(chrom, cs-1, cs, score))

    for data in datasets:
        datasets[data]['out'].close()
