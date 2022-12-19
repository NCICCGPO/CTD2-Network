#!/usr/bin/env python

import argparse
import os
import fnmatch
import math

import numpy as np
# from scipy import interpolate
from scipy import stats as spstats


import logging
logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logger = logging.getLogger('ibis2d')
logger.setLevel(logging.INFO)

def float_safe(strval):
    if (strval == ''):
        strval = '0.0'
    try:
      val = float(strval)
      return val
    except:
      # logger.info('converting %s to 0.0', strval)
      return 0.0
        
def read_rows(filename):
    fp = open(filename, 'r')
    header_toks = fp.readline().rstrip('\n').split('\t')
    ntoks = len(header_toks)
    logger.info('header: %d toks, %s', len(header_toks), str(header_toks))
    data_rows = [ ]
    i = 0
    for line in fp:
        i += 1
        toks = line.rstrip('\n').split('\t')
        assert(len(toks) == ntoks),'bad token count: line %d has %d toks %s' % (i, len(toks), line)
        data_rows.append(toks)
    fp.close()
    return(header_toks, data_rows)

def print_rows(filename, header, rows):
    logger.info('writing %d rows to %s', len(rows), filename)
    fp = open(filename, 'w')
    fp.write('\t'.join(header) + '\n')
    for row in rows:
        fp.write('\t'.join(row) + '\n')
    fp.close()

def write_list(filename, mylist, header=None):
    fp = open(filename, 'w')
    if header is not None:
        fp.write('\t'.join(header) + '\n' )
    for rec in mylist:
        fp.write('\t'.join([str(x) for x in rec]) + '\n')
    fp.close()
    
# gmtfile has name, comment, list of genes, all tab-delimited, ragged length
def read_gmt(gmtfile):
    fp = open(gmtfile, 'r')
    gmt_rows = [ ]
    for line in fp:
        toks = line.rstrip('\n').split('\t')
        assert(len(toks) >= 3)
        symbols = [ s.upper() for s in toks[2:] ]
        gmt_rows.append([ toks[0], toks[1], symbols ])
    return(gmt_rows)
    
def analyze_geneset(geneset, gene_list, gene_zscore):
    nset = len(geneset)
    tmpdict = dict(zip(geneset, [True]*nset))
    inset = [ abs(gene_zscore[g]) for g in gene_list if g in tmpdict ]
    outset = [ abs(gene_zscore[g]) for g in gene_list if g not in tmpdict ]
    nin = len(inset)
    nout = len(outset)
    (tstat, pvalue) = (0, 0.5)
    if ((nin > 2) and (nout > 2)):
        (tstat, pvalue) = spstats.ttest_ind(inset, outset)
    if (tstat > 0):
        pvalue = 0.5 * pvalue
    else:
        pvalue = 1 - 0.5 * pvalue
    z_list = [ gene_zscore[g] for g in gene_list if g in tmpdict ]
    if nin>0:
      z_mean = sum(z_list)/float(nin)
    else:
      z_mean = 0
    g_list = [ g for g in gene_list if g in tmpdict ]
    zg = sorted( zip(z_list, g_list) )
    top5 = zg[0:min(5, nin)]
    bottom5 = zg[max(0,nin-6):max(0,(nin-1))]
    return(pvalue, tstat, nset, nin, nout, z_mean, top5, bottom5)

def main():

    parser = argparse.ArgumentParser(description='Run gene set enrichment analysis using t-test',
                                     epilog='Sample call: see run.sh',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('zscorefile', help='file with pvalues')
    parser.add_argument('gmtfile', help='file with gene sets in gmt format')
    parser.add_argument('outputfile', help='file to write output')
    parser.add_argument('--method', help='[raw|rank] for raw z-scores or convert to rank', default='raw')
    args = parser.parse_args()
    
    # read the files
    (gene_header, gene_rows) = read_rows(args.zscorefile)
    
    SYMBOL = 1
    ZSCORE = 2
    assert(gene_header[ZSCORE - 1] == 'zscore')
    gene_list_mixedcase = [ row[0] for row in gene_rows ]
    gene_list = [ g.upper() for g in gene_list_mixedcase ]
    zscore_orig = [ float(row[ZSCORE - 1 ]) for row in gene_rows ]
    zscore_list = [ z for z in zscore_orig ]
    if args.method == 'rank':
        logger.info('converting zscores to ranks')
        zscore_abs = [ abs(z) for z in zscore_orig ]
        rank_list = spstats.rankdata(zscore_abs)
        sgn_list = np.sign(zscore_orig)
        zscore_list = [ x[0]*x[1] for x in zip(sgn_list, rank_list) ]
        
    
    gene_zscore = dict( zip(gene_list, zscore_list) )
    
    gmt_rows = read_gmt(args.gmtfile)
    
    fp = open(args.outputfile, 'w')
    fp.write('\t'.join(['geneset','comment','pvalue','tstat','nset','nin','nout','zmean','top5','bottom5'])+'\n')
    for (name, comment, geneset) in gmt_rows:
        (pvalue, tstat, nset, nin, nout, zmean, top5, bot5) = analyze_geneset(geneset, gene_list, gene_zscore)
        fp.write('\t'.join([ name, comment, str(pvalue), str(tstat),
                            str(nset), str(nin), str(nout),
                            str(zmean),
                            str(top5), str(bot5) ]) + '\n')
    fp.close()
        

if __name__ == "__main__":
    main()
