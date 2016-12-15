import argparse
import os
import sys
import logging
import time
import random
import colorsys
from Kleat import *
from Matrix import *

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
    parser.add_argument('kleats', nargs='+', help='One or more KLEAT files to analyze')
    parser.add_argument('-e', '--exclude', default=[], nargs='+', help='Tissue types to exclude')
    parser.add_argument('-s', '--samplemap', default=os.path.dirname(os.path.realpath(__file__))+'/samplemap', help='A file mapping sample names to tissue site')
    parser.add_argument('-n', '--normalize', action='store_true', help='Enable this to randomly subsample all tissues to equal the minimum number of samples')
    parser.add_argument("-l", "--log", dest="logLevel", default='WARNING', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the logging level. Default = "WARNING"')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))
    
    args = parser.parse_args()
    
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
                'file': open(os.path.join(args.outdir, '{}.bg'.format(tissue)), 'w'),
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
    for i,tissue in enumerate(tissues):
        tissues[tissue]['file'].write('track type=\'bedGraph\' visibility=\'2\' name=\'{}\' color=\'{}\'\n'.format(tissue, (',').join([str(x*255) for x in rgb[i]])))

    # Normalize
    #if args.normalize:
    #    smallest = min([len(tissues[tissue]['samples']) for tissue in tissues])
    #    logging.debug('smallest: {}'.format(smallest))
    #    keep = set([x for y in [random.sample(tissues[tissue]['samples'], smallest) for tissue in tissues] for x in y])
        #keep = set([x for x in random.sample(tissues[tissue]['samples'], smallest) for tissue in tissues])
    #    kleats = [x for x in kleats if x.name in keep]
    
    kleats = Kleat.groupKleat(kleats)

    for chrom in kleats:
        for gene in kleats[chrom]:
            kleats[chrom][gene] = iterAHC(kleats[chrom][gene])
            clusters = kleats[chrom][gene]
            lclusters = len(clusters)
            mtx = Matrix(N, lclusters)
            centroids = []
            for cluster in clusters:
                values = [x.cleavage_site for x in cluster]
                centroids.append(centroid(values))
            for i in xrange(mtx.m):
                tissue = tissues.keys()[i]
                for j in xrange(mtx.n):
                    values = [x.cleavage_site for x in clusters[j] if mapping[x.name] == tissue]
                    if not values:
                        continue
                    cs = centroid(values)
                    score = len(values)
                    factor = len(tissues[tissue]['samples'])
                    if args.normalize:
                        score /= float(factor)
                    tissues[tissue]['file'].write('{}\t{}\t{}\t{}\n'.format(chrom, cs-1, cs, score))
                    mtx.mtx[i][j] = score
            mtx.rownames = [t.replace(' ', '_') for t in tissues.keys()]
            mtx.colnames = ['c{}_{}'.format(i, int(x)) for i,x in enumerate(centroids)]
            with open(os.path.join(args.outdir, gene), 'w') as f:
                f.write(str(mtx))
#            deltas = Matrix(N, lclusters-1)
#            for i in xrange(deltas.m):
#                for j in xrange(deltas.n):
#                    value = mtx.mtx[i][j+1] - mtx.mtx[i][j]
#                    deltas.mtx[i][j] = value
