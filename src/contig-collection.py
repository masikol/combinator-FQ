# -*- coding: utf-8 -*-


# Readable indices
FULL_NAME, NAME, LEN, COV, GC, START, RC_START, END, RC_END, START_MATCH, END_MATCH, MULT = range(12)


def get_contig_collection(infpath):

    # Iterate over contigs and form contig_collection dictionary
    for i, contig_name in enumerate(contigs_names):

        contig_name = contig_name.strip("_")

        contig_len = len(contigs_seqs[i])
        total_length += contig_len

        # Retrieve coverage information from SPAdes fasta header
        if not re.search(spades_patt, contig_name) is None:

            # Parse fasta header:
            cov = round(float(contig_name.split('_')[5]), 2) # get coverage

            # Collecting coverage statistics
            avg_coverage += cov
            min_coverage = min(min_coverage, cov)
            max_coverage = max(max_coverage, cov)

            name = 'NODE_' + contig_name.split('_')[1] # get name in 'NODE_<NUMBER>' format

        # 'combinator-FQ' works for SPAdes and A5 assemblies for now,
        #   and only SPAdes specifies coverage in fasta header.
        # Therefore we'll just write minus character fro coverage,
        #   if contigs were assembled not by SPAdes:
        else:
            cov = '-'
            name = contig_name # use full header as name
        # end if

        # Calculate GC-content
        gc_content = 0
        for up_base, low_base in zip(('G', 'C', 'S'),('g', 'c', 's')):
            gc_content += contigs_seqs[i].count(up_base) + contigs_seqs[i].count(low_base)
        # end for
        gcContent = round((gc_content / len(contigs_seqs[i]) * 100), 2)

        contig_collection.append([
            contig_name,
            name,
            contig_len,
            cov,
            gcContent,
            contigs_seqs[i][:maxk],
            rc(contigs_seqs[i][:maxk]),
            contigs_seqs[i][-maxk:],
            rc(contigs_seqs[i][-maxk:]),
            list(),
            list(),
            None
        ])
    # end for

# end def get_contig_collection


def _fasta_generator(infpath):

    curr_seq_name = ''
    curr_seq = ''

    with open(infpath, 'r') as infile:

        for line in infile: # line by line
            line = fmt_func(line)
            if line[0] == '>':
                contigs_names.append(line[1:])
                contigs_seqs.append(seq.upper())
                seq = ""
            else:
                seq += line
            # end if
        # end for
    # end with

# end def _fasta_generator
