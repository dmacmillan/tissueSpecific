import os

class Kleat:

    columns = ('gene', 'transcript', 'transcript_strand', 'coding', 'contig', 'chromosome', 'cleavage_site', 'within_UTR', 'distance_from_annotated_site', 'ESTs', 'length_of_tail_in_contig', 'number_of_tail_reads', 'number_of_bridge_reads', 'max_bridge_read_tail_length', 'bridge_read_identities', 'tail+bridge_reads', 'number_of_link_pairs', 'max_link_pair_length', 'link_pair_identities', 'hexamer_loc+id', '3UTR_start_end', 'flag')
    
    @staticmethod
    def strToInt(string):
        try:
            return int(string)
        except ValueError:
            return '-'

    def __init__(self, gene, transcript, transcript_strand, coding, contig, chromosome, cleavage_site, within_UTR, distance_from_annotated_site, ESTs, length_of_tail_in_contig, number_of_tail_reads, number_of_bridge_reads, max_bridge_read_tail_length, bridge_read_identities, tail_and_bridge_reads, number_of_link_pairs, max_link_pair_length, link_pair_identities, pas, utr3, flag=None):
        self.gene = gene
        self.transcript = transcript
        self.transcript_strand = transcript_strand
        self.coding = coding
        self.contig = contig
        self.chromosome = chromosome
        self.cleavage_site = int(cleavage_site)
        if (within_UTR == 'no'):
            self.within_UTR = False
        else:
            self.within_UTR = True
        self.distance_from_annotated_site = Kleat.strToInt(distance_from_annotated_site)
        self.ESTs = ESTs
        self.length_of_tail_in_contig = Kleat.strToInt(length_of_tail_in_contig)
        self.number_of_tail_reads = Kleat.strToInt(number_of_tail_reads)
        self.number_of_bridge_reads = Kleat.strToInt(number_of_bridge_reads)
        self.max_bridge_read_tail_length = Kleat.strToInt(max_bridge_read_tail_length)
        self.bridge_read_identities = bridge_read_identities.split(',')
        self.tail_and_bridge_reads = Kleat.strToInt(tail_and_bridge_reads)
        self.number_of_link_pairs = Kleat.strToInt(number_of_link_pairs)
        self.max_link_pair_length = Kleat.strToInt(max_link_pair_length)
        self.link_pair_identities = link_pair_identities
        self.pas = self.utr3 = None
        self.flag = flag
        try:
            self.pas = [int(x) for x in pas.split(':')]
        except ValueError:
            pass
        try:
            self.utr3 = [int(x) for x in utr3.split('-')]
        except ValueError:
            pass

    def __str__(self):
        pas = '-'
        if self.pas:
            pas = (':').join([str(x) for x in self.pas])
        if self.within_UTR:
            within_UTR = 'yes'
        else:
            within_UTR = 'no'
        self.bridge_read_identities = (',').join(self.bridge_read_identities)
        atts = [self.gene, self.transcript, self.transcript_strand, self.coding, self.contig, self.chromosome, self.cleavage_site, within_UTR, self.distance_from_annotated_site, self.ESTs, self.length_of_tail_in_contig, self.number_of_tail_reads, self.number_of_bridge_reads, self.max_bridge_read_tail_length, self.bridge_read_identities, self.tail_and_bridge_reads, self.number_of_link_pairs, self.max_link_pair_length, self.link_pair_identities, pas, self.utr3, self.flag]
        atts = [str(x) if x is not None else '-' for x in atts]
        return ('\t').join(atts)

    def __repr__(self):
        pas = '-'
        if self.pas:
            pas = (':').join([str(x) for x in self.pas])
        atts = [self.gene, self.transcript, self.transcript_strand, self.coding, self.contig, self.chromosome, self.cleavage_site, self.within_UTR, self.distance_from_annotated_site, self.ESTs, self.length_of_tail_in_contig, self.number_of_tail_reads, self.number_of_bridge_reads, self.max_bridge_read_tail_length, self.bridge_read_identities, self.tail_and_bridge_reads, self.number_of_link_pairs, self.max_link_pair_length, self.link_pair_identities, pas, self.utr3]
        atts = [str(x) if x is not None else '-' for x in atts]
        return ('\t').join(atts)

    @staticmethod
    def getCleavageSites(kleat):
        with open(kleat, 'r') as f:
            header = f.readline()
            results = set([int(x[6]) for x in f])
        return results

    @staticmethod
    def KleatFile(kleat):
        with open(kleat, 'r') as f:
            header = f.readline()
            for line in f:
                yield Kleat(*line.strip().split('\t'))
        
    @staticmethod
    def parseKleatFast(kleat, name=None):
        results = []
        with open(kleat, 'r') as f:
            header = f.readline()
            for line in f:
                result = Kleat(*line.strip().split('\t'))
                if name:
                    result.name = name
                results.append(result)
        return results

    @staticmethod
    def parseKleat(kleat, 
                   max_dist_ann = float('inf'), 
                   min_len_tail_contig = None,
                   min_num_tail_reads = None,
                   min_num_bridge_reads = None,
                   #min_bridge_read_tail_len = None,
                   has_pas = None,
                   discard_added = False):
        results = []
        with open(kleat, 'r') as f:
            header = f.readline()
            for line in f:
                result = Kleat(*line.strip().split('\t'))
                if discard_added and result.flag == '2':
                    continue
                if has_pas:
                    if not result.pas:
                        continue
                if max_dist_ann and (result.distance_from_annotated_site != '-'):
                    if max_dist_ann < result.distance_from_annotated_site:
                        continue
                if min_len_tail_contig and (result.length_of_tail_in_contig != '-'):
                    if min_len_tail_contig > result.length_of_tail_in_contig:
                        continue
                if min_num_tail_reads and (result.number_of_tail_reads != '-'):
                    if min_num_tail_reads > result.number_of_tail_reads:
                        continue
                if min_num_bridge_reads and (result.number_of_bridge_reads != '-'):
                    if min_num_bridge_reads > result.number_of_bridge_reads:
                        continue
                results.append(result)
        return results

    @staticmethod
    def groupKleat(parsed, results={}):
        for r in parsed:
            if r.chromosome not in results:
                results[r.chromosome] = {r.gene: [r]}
            elif r.gene not in results[r.chromosome]:
                results[r.chromosome][r.gene] = [r]
            else:
                results[r.chromosome][r.gene].append(r)
        return results
