import Levenshtein
import argparse, gzip, os

max_Ldist = 1 # can't be changed to larger than 1 without requiring code changes elsewhere

def find_startpoint(read, query, def_query_from_5prim):
    for query_from_5prim in (def_query_from_5prim, def_query_from_5prim-1, def_query_from_5prim+1):
        read_part = read[query_from_5prim:query_from_5prim+len(query)]
        Ldist = Levenshtein.distance(read_part, query)
        if Ldist <= max_Ldist or Ldist > max_Ldist+2: break
    if Ldist <= max_Ldist:
        return Ldist, query_from_5prim, query_from_5prim + len(query)
    else:
        return Ldist, None, 0

def reverseDNA(seq_in):
	""" returns nucleotide sequence string """
	sequencetools_reverseDNAdict = {"A":"T", "C":"G", "G":"C", "T":"A", "R":"Y","Y":"R","K":"M","M":"K","S":"S","W":"W","B":"V","D":"H","H":"D","V":"B","N":"N", "a":"t","c":"g","g":"c","t":"a","n":"n","\n":"\n"}
	seq_out = ""
	for bi in range(len(seq_in)):
		seq_out = sequencetools_reverseDNAdict[seq_in[bi]] + seq_out
	return seq_out

def match_barcode(barcode_raw, barcodes):
    # returns the barcode in barcodes most similar to barcode_raw, or throw a ValueError if none is similar enough
    dists = [(Levenshtein.distance(barcode, barcode_raw), barcode) for barcode in barcodes]
    dists.sort()
    if dists[0][0] <= max_Ldist and dists[1][0] - dists[0][0] >= 2:
        return dists[0][1]
    raise ValueError

def cut_read(read, queries):
    # returns sequence until it matches either of the queries
    
    # check for full matches
    for query in queries:
        try:
            pos = read.index(query)
        except ValueError:
            continue
        return read[pos]
    
    # check for match to first half
    for query in queries:
        try:
            pos = read.index(query[:len(query)//2])
        except ValueError:
            continue
        other_half = query[len(query)//2:]
        dist_other = Levenshtein.distance(read[pos+len(other_half):pos+len(query)], other_half)
        if dist_other <= max_Ldist:
            return read[pos]
    
     # check for match to second half
    for query in queries:
        this_half = query[len(query)//2:]
        try:
            pos = read.index(this_half)
        except ValueError:
            continue
        other_half = query[:len(query)//2]
        dist_other = Levenshtein.distance(read[pos-len(other_half):pos+len(this_half)], other_half)
        if dist_other <= max_Ldist:
            return read[pos-len(other_half)]
        
    return read

def write_fastq_read(outfh, seqname, read, qual):
    print('@'+seqname, file=outfh)
    print(read, file=outfh)
    print('+', file=outfh)
    print(qual, file=outfh)

def parse_pasted_step1(barcodes_pasted):
    ret = []
    for ri, line in enumerate(barcodes_pasted.splitlines()):
        if ri < 2: continue # skip the 2 header lines
        p = line.strip().split('\t')
        if len(p) != 2: continue
        ret.append(p[1])
    return ret

row_barcodes_pasted = '''
rowname	barcode
PlateNN1A	GCAGAGATAAGC
PlateNN1B	GCAGAGATGCAC
PlateNN1C	GCAGAGACTCAG
PlateNN1D	GCAGAGAGGAAT
PlateNN1E	GCAGAGACGAGG
PlateNN1F	GCAGAGAAGGAG
PlateNN1G	GCAGAGATGTTG
PlateNN1H	GCAGAGACAACT
PlateNN2A	TCGAAGATAAGC
PlateNN2B	TCGAAGATGCAC
PlateNN2C	TCGAAGACTCAG
PlateNN2D	TCGAAGAGGAAT
PlateNN2E	TCGAAGACGAGG
PlateNN2F	TCGAAGAAGGAG
PlateNN2G	TCGAAGATGTTG
PlateNN2H	TCGAAGACAACT
PlateNN3A	AACAAGATAAGC
PlateNN3B	AACAAGATGCAC
PlateNN3C	AACAAGACTCAG
PlateNN3D	AACAAGAGGAAT
PlateNN3E	AACAAGACGAGG
PlateNN3F	AACAAGAAGGAG
PlateNN3G	AACAAGATGTTG
PlateNN3H	AACAAGACAACT
PlateNN4A	GGTGCGATAAGC
PlateNN4B	GGTGCGATGCAC
PlateNN4C	GGTGCGACTCAG
PlateNN4D	GGTGCGAGGAAT
PlateNN4E	GGTGCGACGAGG
PlateNN4F	GGTGCGAAGGAG
PlateNN4G	GGTGCGATGTTG
PlateNN4H	GGTGCGACAACT
PlateNN5A	TTGGTGATAAGC
PlateNN5B	TTGGTGATGCAC
PlateNN5C	TTGGTGACTCAG
PlateNN5D	TTGGTGAGGAAT
PlateNN5E	TTGGTGACGAGG
PlateNN5F	TTGGTGAAGGAG
PlateNN5G	TTGGTGATGTTG
PlateNN5H	TTGGTGACAACT
PlateNN6A	CATTCGATAAGC
PlateNN6B	CATTCGATGCAC
PlateNN6C	CATTCGACTCAG
PlateNN6D	CATTCGAGGAAT
PlateNN6E	CATTCGACGAGG
PlateNN6F	CATTCGAAGGAG
PlateNN6G	CATTCGATGTTG
PlateNN6H	CATTCGACAACT
PlateNN7A	ATTGGGATAAGC
PlateNN7B	ATTGGGATGCAC
PlateNN7C	ATTGGGACTCAG
PlateNN7D	ATTGGGAGGAAT
PlateNN7E	ATTGGGACGAGG
PlateNN7F	ATTGGGAAGGAG
PlateNN7G	ATTGGGATGTTG
PlateNN7H	ATTGGGACAACT
PlateNN8A	CGGTTGATAAGC
PlateNN8B	CGGTTGATGCAC
PlateNN8C	CGGTTGACTCAG
PlateNN8D	CGGTTGAGGAAT
PlateNN8E	CGGTTGACGAGG
PlateNN8F	CGGTTGAAGGAG
PlateNN8G	CGGTTGATGTTG
PlateNN8H	CGGTTGACAACT
PlateNN9A	ATCCTGATAAGC
PlateNN9B	ATCCTGATGCAC
PlateNN9C	ATCCTGACTCAG
PlateNN9D	ATCCTGAGGAAT
PlateNN9E	ATCCTGACGAGG
PlateNN9F	ATCCTGAAGGAG
PlateNN9G	ATCCTGATGTTG
PlateNN9H	ATCCTGACAACT
PlateNN10A	ATGTCGATAAGC
PlateNN10B	ATGTCGATGCAC
PlateNN10C	ATGTCGACTCAG
PlateNN10D	ATGTCGAGGAAT
PlateNN10E	ATGTCGACGAGG
PlateNN10F	ATGTCGAAGGAG
PlateNN10G	ATGTCGATGTTG
PlateNN10H	ATGTCGACAACT
PlateNN11A	TCACGGATAAGC
PlateNN11B	TCACGGATGCAC
PlateNN11C	TCACGGACTCAG
PlateNN11D	TCACGGAGGAAT
PlateNN11E	TCACGGACGAGG
PlateNN11F	TCACGGAAGGAG
PlateNN11G	TCACGGATGTTG
PlateNN11H	TCACGGACAACT
PlateNN12A	AGACCGATAAGC
PlateNN12B	AGACCGATGCAC
PlateNN12C	AGACCGACTCAG
PlateNN12D	AGACCGAGGAAT
PlateNN12E	AGACCGACGAGG
PlateNN12F	AGACCGAAGGAG
PlateNN12G	AGACCGATGTTG
PlateNN12H	AGACCGACAACT
PlateNN13A	CCCCAGATAAGC
PlateNN13B	CCCCAGATGCAC
PlateNN13C	CCCCAGACTCAG
PlateNN13D	CCCCAGAGGAAT
PlateNN13E	CCCCAGACGAGG
PlateNN13F	CCCCAGAAGGAG
PlateNN13G	CCCCAGATGTTG
PlateNN13H	CCCCAGACAACT
PlateNN14A	GCGCTGATAAGC
PlateNN14B	GCGCTGATGCAC
PlateNN14C	GCGCTGACTCAG
PlateNN14D	GCGCTGAGGAAT
PlateNN14E	GCGCTGACGAGG
PlateNN14F	GCGCTGAAGGAG
PlateNN14G	GCGCTGATGTTG
PlateNN14H	GCGCTGACAACT
PlateNN15A	TCCTTGATAAGC
PlateNN15B	TCCTTGATGCAC
PlateNN15C	TCCTTGACTCAG
PlateNN15D	TCCTTGAGGAAT
PlateNN15E	TCCTTGACGAGG
PlateNN15F	TCCTTGAAGGAG
PlateNN15G	TCCTTGATGTTG
PlateNN15H	TCCTTGACAACT
PlateNN16A	TATATGATAAGC
PlateNN16B	TATATGATGCAC
PlateNN16C	TATATGACTCAG
PlateNN16D	TATATGAGGAAT
PlateNN16E	TATATGACGAGG
PlateNN16F	TATATGAAGGAG
PlateNN16G	TATATGATGTTG
PlateNN16H	TATATGACAACT
PlateNN17A	CGTAAGATAAGC
PlateNN17B	CGTAAGATGCAC
PlateNN17C	CGTAAGACTCAG
PlateNN17D	CGTAAGAGGAAT
PlateNN17E	CGTAAGACGAGG
PlateNN17F	CGTAAGAAGGAG
PlateNN17G	CGTAAGATGTTG
PlateNN17H	CGTAAGACAACT
PlateNN18A	AAGGTGATAAGC
PlateNN18B	AAGGTGATGCAC
PlateNN18C	AAGGTGACTCAG
PlateNN18D	AAGGTGAGGAAT
PlateNN18E	AAGGTGACGAGG
PlateNN18F	AAGGTGAAGGAG
PlateNN18G	AAGGTGATGTTG
PlateNN18H	AAGGTGACAACT
PlateNN19A	AGCTCGATAAGC
PlateNN19B	AGCTCGATGCAC
PlateNN19C	AGCTCGACTCAG
PlateNN19D	AGCTCGAGGAAT
PlateNN19E	AGCTCGACGAGG
PlateNN19F	AGCTCGAAGGAG
PlateNN19G	AGCTCGATGTTG
PlateNN19H	AGCTCGACAACT
PlateNN20A	CCTGCGATAAGC
PlateNN20B	CCTGCGATGCAC
PlateNN20C	CCTGCGACTCAG
PlateNN20D	CCTGCGAGGAAT
PlateNN20E	CCTGCGACGAGG
PlateNN20F	CCTGCGAAGGAG
PlateNN20G	CCTGCGATGTTG
PlateNN20H	CCTGCGACAACT
PlateNN21A	GTATCGATAAGC
PlateNN21B	GTATCGATGCAC
PlateNN21C	GTATCGACTCAG
PlateNN21D	GTATCGAGGAAT
PlateNN21E	GTATCGACGAGG
PlateNN21F	GTATCGAAGGAG
PlateNN21G	GTATCGATGTTG
PlateNN21H	GTATCGACAACT
PlateNN22A	TATGAGATAAGC
PlateNN22B	TATGAGATGCAC
PlateNN22C	TATGAGACTCAG
PlateNN22D	TATGAGAGGAAT
PlateNN22E	TATGAGACGAGG
PlateNN22F	TATGAGAAGGAG
PlateNN22G	TATGAGATGTTG
PlateNN22H	TATGAGACAACT
PlateNN23A	CACACGATAAGC
PlateNN23B	CACACGATGCAC
PlateNN23C	CACACGACTCAG
PlateNN23D	CACACGAGGAAT
PlateNN23E	CACACGACGAGG
PlateNN23F	CACACGAAGGAG
PlateNN23G	CACACGATGTTG
PlateNN23H	CACACGACAACT
PlateNN24A	ACACTGATAAGC
PlateNN24B	ACACTGATGCAC
PlateNN24C	ACACTGACTCAG
PlateNN24D	ACACTGAGGAAT
PlateNN24E	ACACTGACGAGG
PlateNN24F	ACACTGAAGGAG
PlateNN24G	ACACTGATGTTG
PlateNN24H	ACACTGACAACT
PlateNN25A	ACTACGATAAGC
PlateNN25B	ACTACGATGCAC
PlateNN25C	ACTACGACTCAG
PlateNN25D	ACTACGAGGAAT
PlateNN25E	ACTACGACGAGG
PlateNN25F	ACTACGAAGGAG
PlateNN25G	ACTACGATGTTG
PlateNN25H	ACTACGACAACT
PlateNN26A	GTTACGATAAGC
PlateNN26B	GTTACGATGCAC
PlateNN26C	GTTACGACTCAG
PlateNN26D	GTTACGAGGAAT
PlateNN26E	GTTACGACGAGG
PlateNN26F	GTTACGAAGGAG
PlateNN26G	GTTACGATGTTG
PlateNN26H	GTTACGACAACT
'''

col_barcodes_pasted = '''
Column1	Column2
BC1	GTTCA
BC2	CAGGA
BC3	TTATA
BC4	CCTGT
BC5	ACCGC
BC6	ACTTA
BC7	GCTAG
BC8	GACGT
BC9	GGCTA
BC10	GAATG
BC11	CCAAC
BC12	GAGAC
'''

def parse_pasted_step2(barcodes_pasted):
    ret = []
    for ri, line in enumerate(barcodes_pasted.splitlines()):
        if ri < 2: continue # skip the 2 header lines
        p = line.strip().split('\t')
        if len(p) != 2: continue
        ret.append(p)
    return ret

def load4lines(infh):
    return next(infh), next(infh), next(infh), next(infh)

def getBC(fourlines):
    return fourlines[1].rstrip()

def writelines(outfh, lines):
    for line in lines:
        outfh.write(line)

def remove(filepath):
    try:
        os.remove(filepath)
    except:
        pass

def step1(fastq_R1_in, fastq_R2_in, fastq_R1_out, fastq_R2_out, fastq_rowBC_out, fastq_colBC_out, use_gzip, minlen):
    
    open_func0 = gzip.open if use_gzip[0] else open
    open_func1 = open_check_gzip if use_gzip[1] else open_check
    
    # seems to be artificial and should be cut
    R1_queries = ['CCAGGGTTTTCCCAGTCACGAC'] # hard-coded later on to only work for one sequence
    R1_query_dist = 14
    R1_barcode_size = 12
    
    # seem to be biological and be kept in read
    R2_queries = ['GTCACTGGATTTAGAGTCTCTCAG', 'GAGATCTCTGCTTCTGATGGCTC'] # hard-coded later on to only work for two sequences
    R2_query_dist = 7
    R2_barcode_size = 5
    
    # load lists of barcodes
    row_barcodes = parse_pasted_step1(row_barcodes_pasted)
    col_barcodes = parse_pasted_step1(col_barcodes_pasted)
    
    # find reverse complements
    R1_rev_queries = [reverseDNA(seq) for seq in R1_queries]
    R2_rev_queries = [reverseDNA(seq) for seq in R2_queries]
    
    
    no_R1_query_found = 0
    no_R2_query_found = 0
    too_short_R1 = 0
    too_short_R2 = 0
    no_matched_row_barcode = 0
    no_matched_col_barcode = 0
    read_pair_passed = 0
    
    tmp_files_opened = [False, False, False, False]
    try:
        with open_func0(fastq_R1_in, 'rt') as infh1:
            with open_func0(fastq_R2_in, 'rt') as infh2:
                with open_func1(fastq_R1_out, 'wt') as outfhR1:
                    tmp_files_opened[0] = True
                    with open_func1(fastq_R2_out, 'wt') as outfhR2:
                        tmp_files_opened[1] = True
                        with open_func1(fastq_rowBC_out, 'wt') as outfhBC1:
                            tmp_files_opened[2] = True
                            with open_func1(fastq_colBC_out, 'wt') as outfhBC2:
                                tmp_files_opened[3] = True
                                for li, (lineR1, lineR2) in enumerate(zip(infh1, infh2)):
                                    linetype = li % 4
                                    if linetype == 0: # the name line of a fastq read
                                        seqname = lineR1[1:].split()[0]
                                        seqnameR2 = lineR2[1:].split()[0]
                                        if seqname != seqnameR2:
                                            print(seqname)
                                            print(seqnameR2)
                                            raise ValueError
                                    elif linetype == 1: # the sequence line of a fastq read
                                        should_write = False
                                        readR1 = lineR1.rstrip()
                                        
                                        if readR1[R1_query_dist:].startswith(R1_queries[0]): # check for exact match
                                            R1_query_from_5prim = R1_query_dist
                                        else: # check for inexact match, including being one off in start position
                                            R1_query_from_5prim = find_startpoint(readR1, R1_queries[0], R1_query_dist)[1]
                                        
                                        
                                        if R1_query_from_5prim is None:
                                            no_R1_query_found += 1
                                            continue
                                        # found the fixed part of the plate/row sequence
                                        
                                        R1_query_from_3prim = R1_query_from_5prim + len(R1_queries[0])
                                        
                                        readR2 = lineR2.rstrip()
                                        if readR2[R2_query_dist:].startswith(R2_queries[0]): # check for exact match
                                            R2_query_from_5prim = R2_query_dist
                                            R2_query_from_3prim = R2_query_dist + len(R2_queries[0])
                                        elif readR2[R2_query_dist:].startswith(R2_queries[1]): # check for exact match
                                            R2_query_from_5prim = R2_query_dist
                                            R2_query_from_3prim = R2_query_dist + len(R2_queries[1])
                                        else: # check for inexact match, including being one off in start position
                                            Ldist0, R2_query_from_5prim_0, R2_query_from_3prim_0 = find_startpoint(readR2, R2_queries[0], R2_query_dist)
                                            Ldist1, R2_query_from_5prim_1, R2_query_from_3prim_1 = find_startpoint(readR2, R2_queries[1], R2_query_dist)
                                            R2_query_from_5prim, R2_query_from_3prim = (R2_query_from_5prim_0, R2_query_from_3prim_0) if Ldist0 < Ldist1 else (R2_query_from_5prim_1, R2_query_from_3prim_1)
                                        
                                        if R2_query_from_5prim is None:
                                            no_R2_query_found += 1
                                            continue
                                        # found the fixed part of the column sequence
                                        
                                        # cut out the sequence of the barcodes
                                        barcodeR1_raw = readR1[R1_query_from_5prim-R1_barcode_size:R1_query_from_5prim]
                                        barcodeR2_raw = readR2[R2_query_from_5prim-R2_barcode_size:R2_query_from_5prim]
                                        
                                        # get the sequence of the matching barcode of the list
                                        try:
                                            barcodeR1_matched = match_barcode(barcodeR1_raw, row_barcodes)
                                        except ValueError:
                                            no_matched_row_barcode += 1
                                            continue
                                        try:
                                            barcodeR2_matched = match_barcode(barcodeR2_raw, col_barcodes)
                                        except ValueError:
                                            no_matched_col_barcode += 1
                                            continue
                                        
                                        
                                        # remove artifical sequences
                                        readR1_cut = cut_read(readR1[R1_query_from_3prim:], R2_rev_queries)
                                        if len(readR1_cut) < minlen:
                                            too_short_R1 += 1
                                            continue
                                        readR2_cut = cut_read(readR2[R2_query_from_5prim:], R1_rev_queries)
                                        if len(readR2_cut) < minlen:
                                            too_short_R2 += 1
                                            continue
                                        
                                        # tell to output barcodeR1_raw, barcodeR2_raw, readR1_cut, readR2_cut to four different files
                                        should_write = True
                                        read_pair_passed += 1
                                    elif linetype == 3 and should_write: # the quality line of a fastq read, only run if told to output read pair
                                        qualR1 = lineR1.rstrip()
                                        qualR2 = lineR2.rstrip()
                                        
                                        write_fastq_read(outfhR1, seqname, readR1_cut, qualR1[R1_query_from_3prim:R1_query_from_3prim+len(readR1_cut)])
                                        write_fastq_read(outfhR2, seqname, readR2_cut, qualR2[R2_query_from_5prim:R2_query_from_5prim+len(readR2_cut)])
                                        write_fastq_read(outfhBC1, seqname, barcodeR1_matched, qualR1[R1_query_from_5prim-R1_barcode_size:R1_query_from_5prim])
                                        write_fastq_read(outfhBC2, seqname, barcodeR2_matched, qualR2[R2_query_from_5prim-R2_barcode_size:R2_query_from_5prim])
        
        print('No R1 match to query sequence:', no_R1_query_found)
        print('No R2 match to query sequence:', no_R2_query_found)
        print('No row barcode match:', no_matched_row_barcode)
        print('No column barcode match:', no_matched_col_barcode)
        print('Too short R1 insert:', too_short_R1)
        print('Too short R2 insert:', too_short_R2)
        print('Written read pairs to intermediary files:', read_pair_passed)
    except:
        for was_opened, path in zip(tmp_files_opened, tmp_fastqs):
            if was_opened:
                remove(path)
        raise


def open_check(path, mode):
    if os.path.exists(path) and 'w' in mode: raise Exception('Delete %s and try again'%path)
    return open(path, mode)

def open_check_gzip(path, mode):
    if os.path.exists(path) and 'w' in mode: raise Exception('Delete %s and try again'%path)
    return gzip.open(path, mode)

def step2(fastq_R1_in, fastq_R2_in, fastq_rowBC_in, fastq_colBC_in, fastq_prefix_out, concurrent, use_gzip):
    open_func1 = gzip.open if use_gzip[1] else open
    open_func2 = gzip.open if use_gzip[2] else open
    
    # load the barcode tables
    barcode_combinations_left = [(rowname +'_'+ colname, rowBC, colBC) for colname, colBC in parse_pasted_step2(col_barcodes_pasted) for rowname, rowBC in parse_pasted_step2(row_barcodes_pasted)]
    
    R1filesuffix = '_R1.fastq' + ('.gz' if use_gzip[2] else '')
    R2filesuffix = '_R2.fastq' + ('.gz' if use_gzip[2] else '')
    
    read_pairs_written = 0
    
    while barcode_combinations_left:
        # pick a limited number of barcode combinations to write fastq files for, to avoid having more files open than the OS ulimit allows
        barcode_combinations = barcode_combinations_left[:concurrent//2]
        barcode_combinations_left = barcode_combinations_left[concurrent//2:]
        
        # create empty fastq files for them
        open_files = {(rowBC, colBC):(open_func2(fastq_prefix_out + name + R1filesuffix, 'wt'), open_func2(fastq_prefix_out + name + R2filesuffix, 'wt')) for name, rowBC, colBC in barcode_combinations}
        
        with open_func1(fastq_R1_in, 'rt') as infhR1:
            with open_func1(fastq_R2_in, 'rt') as infhR2:
                with open_func1(fastq_rowBC_in, 'rt') as infhBC1:
                    with open_func1(fastq_colBC_in, 'rt') as infhBC2:
                        while True:
                            try:
                                R1lines, R2lines, rowBClines, colBClines = load4lines(infhR1), load4lines(infhR2), load4lines(infhBC1), load4lines(infhBC2)
                            except StopIteration:
                                break
                            
                            barcode_pair = getBC(rowBClines), getBC(colBClines)
                            
                            if barcode_pair in open_files:
                                outfhR1, outfhR2 = open_files[barcode_pair] # find the current file pair to write to
                                writelines(outfhR1, R1lines)
                                writelines(outfhR2, R2lines)
                                read_pairs_written += 1
        
        for outfhR1, outfhR2 in open_files.values():
            # flush to files and free the file handles
            outfhR1.close()
            outfhR2.close()
            if os.path.getsize(outfhR1.name) <= 80 and os.path.getsize(outfhR2.name) <= 80:
                # delete empty files, including empte gzipped files which are 60 bytes
                remove(outfhR1.name)
                remove(outfhR2.name)
    
    print('Written read pairs to output files:', read_pairs_written)

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_R1_in')
    parser.add_argument('fastq_R2_in')
    parser.add_argument('fastq_prefix_out')
    parser.add_argument('--concurrent_files_open', type=int, default=100, help="this option sets number of open files at the same time, higher should mean faster, but if it's too close or above 'ulimit -n' it will crash")
    parser.add_argument('--min_read_length_out', type=int, default=5, help="how short output reads are allowed to be")
    parser.add_argument('--gzip_output', choices=['yes', 'no'], help="whether to gzip the output fastq files")
    parser.add_argument('--row_barcodes_file')
    parser.add_argument('--col_barcodes_file')
    o = parser.parse_args()
    
    if o.row_barcodes_file is not None:
        with open(o.row_barcodes_file, 'rt') as infh:
            row_barcodes_pasted = '\n'+infh.read()
    if o.col_barcodes_file is not None:
        with open(o.col_barcodes_file, 'rt') as infh:
            col_barcodes_pasted = '\n'+infh.read()
    
    # decision whether to load/write gzipped fastq files
    use_gzip = [False, False, False] # [input, intermediary files, output]
    use_gzip[0] = o.fastq_R1_in.endswith('.gz') and o.fastq_R2_in.endswith('.gz')
    use_gzip[2] = True if o.gzip_output == 'yes' else False if o.gzip_output == 'no' else use_gzip[0]
    
    # creating 4 separate files for the corrected barcodes and the two split-out inserts
    tmp_formula = o.fastq_prefix_out + '_tmp_%s.fastq'+ ('.gz' if use_gzip[1] else '')
    tmp_fastqs = [tmp_formula%namepart for namepart in ('R1', 'R2', 'rowBC', 'colBC')]
    
    # split into four files, for the corrected barcodes and the two split-out inserts, remove non-matching reads
    step1(o.fastq_R1_in, o.fastq_R2_in, *tmp_fastqs, use_gzip, o.min_read_length_out)
    
    try:
        # go through each read and put it into a fastq file for its row-column combination
        step2(*tmp_fastqs, o.fastq_prefix_out, max(2, o.concurrent_files_open-4), use_gzip)
    finally:
        # delete temporary files
        for path in tmp_fastqs:
            remove(path)