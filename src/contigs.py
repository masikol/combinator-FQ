# -*- coding: utf-8 -*-

import re
import gzip
from functools import partial
from typing import List, Tuple, Dict, Generator, Union, NewType
from typing import Callable, ContextManager, TextIO, BinaryIO


# Dictionary maps complementary bases according to IUPAC:
_COMPL_DICT: Dict[str, str] = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
    'H': 'D', 'V': 'B', 'U': 'A', 'N': 'N'
}


class Contig:
    # Container class representing contig.
    # Fields:
    #  1. `name` -- name.
    #  2. `length` -- length in bp.
    #  3. `cov` -- coverage.
    #  4. `gc_content` -- GC-content.
    #  5. `start` -- prefix of length k of this contig.
    #  6. `rcstart` -- reverse-complement of `start`.
    #  7. `end` -- suffix of length k of this contig.
    #  8. `rcend` -- reverse-complement of `end`.
    #  9. `multplty` -- multiplicity (copies of this contig in the genome).

    def __init__(self, name: str, length: int,
                 cov: float, gc_content: float,
                 start: str, rcstart: str, end: str, rcend: str) -> None:
        self.name = name
        self.length = length
        self.cov = cov
        self.gc_content = gc_content
        self.start = start
        self.rcstart = rcstart
        self.end = end
        self.rcend = rcend
        self.multplty = None
    # end end __init__

#     def __repr__(self):
#         return '<Contig: {}; {} bp; coverage {}; GC content {}%. multiplicity {};\n\
# {}\n\
# {}\n\
# {}\n\
# {}>\n'.format(self.name, self.length, self.cov, self.gc_content, self.multplty,
#     self.start, self.rcstart, self.end, self.rcend)
#     # end def __repr__

# end class Contig

# Custom types declaration
ContigCollection = NewType('ContigCollection', List[Contig])
ContigIndex = NewType('ContigIndex', int)


def get_contig_collection(infpath: str, maxk: int) -> ContigCollection:
    # Function parses a collection of `Contig`s (see above) from a given fasta file.
    # Returns instance of `ContigCollection` (see above).

    # Initialize result colletion
    contig_collection: ContigCollection = list()

    # Pattern that matches SPAdes's header of a seqeunce in FASTA file
    spades_patt: str = r'^NODE_[0-9]+_length_[0-9]+_cov_[0-9,\.]+'

    # Iterate over contigs and form contig_collection
    contig_name: str
    contig_seq:  str
    for contig_name, contig_seq in _fasta_generator(infpath):

        # Remove extra underscores from the name
        contig_name = contig_name.strip("_")

        # Save contig length
        contig_len: int = len(contig_seq)

        # Retrieve coverage information from SPAdes fasta header
        # `combinator-FQ` can parse SPAdes and A5 assemblies correctly,
        #   and only SPAdes specifies coverage in fasta header.
        # Therefore we'll just assign None to coverage ,
        #   if contigs were not assembled with SPAdes.
        cov: float
        name: str
        if not re.search(spades_patt, contig_name) is None:
            # Parse fasta header:
            cov = round(float(contig_name.split('_')[5]), 2) # get coverage
            name = 'NODE_' + contig_name.split('_')[1]       # get name in 'NODE_<NUMBER>' format
        else:
            cov = None
            name = contig_name # use full header as name
        # end if

        # Calculate GC-content
        gc_content: float = 0.0

        up_base: str
        low_base: str
        for up_base, low_base in zip(('G', 'C', 'S'),('g', 'c', 's')):
            gc_content += contig_seq.count(up_base) + contig_seq.count(low_base)
        # end for

        gc_content = round((gc_content / contig_len * 100), 2)

        # Append recently created contig to the `contig_collection`
        contig_collection.append(
            Contig(
                name=name,
                length=contig_len,
                cov=cov,
                gc_content=gc_content,
                start=contig_seq[:maxk].upper(),
                rcstart=_rc(contig_seq[:maxk].upper()),
                end=contig_seq[-maxk:].upper(),
                rcend=_rc(contig_seq[-maxk:].upper())
            )
        )
    # end for

    return contig_collection
# end def get_contig_collection


def assign_multiplty(contig_collection: ContigCollection) -> None:
    # Function assigns multiplicity (copies of this contig in the genome) to contigs.

    # Coverage of 1-st contig can be zero.
    # In this case we cannot calculate multiplicity of contigs.
    # 'calc_multplty' (calculate multiplicity) will indicate whether we can calculate multiplicity.
    calc_multplty: bool = False

    if not contig_collection[0].cov is None:
        if contig_collection[0].cov < 1e-6:
            # Coverage of 1-st contig is zero
            print('\n`{}` has zero coverage (less than 1e-6 actually).'\
                .format(contig_collection[0].name))
            print('Multiplicity of contigs cannot be calculated.\n')
            # We have SPAdes assembly, and coverage of NODE_1 is zero
            calc_multplty = False
        else:
            # Coverage of 1-st contig is not zero
            calc_multplty = True
        # end if
    else:
        # We have non-SPAdes assembly, no coverage provided
        calc_multplty = False
    # end if

    # Calculate multiplicity of contig:
    i: ContigIndex
    if calc_multplty:
        for i in range(len(contig_collection)):
            try:
                multiplicity: float = contig_collection[i].cov / contig_collection[0].cov
            except TypeError:
                contig_collection[i].multplty = 1 # there can be no coverage info for some contigs
            else:
                contig_collection[i].multplty = round(max(multiplicity, 1), 1)
            # end try
        # end if
    else:
        for i in range(len(contig_collection)):
            contig_collection[i].multplty = 1
        # end for
    # end if
# end def assign_multiplty


def _get_compl_base(base: str) -> str:
    # Function returns complement "comrade" of a given base.
    return _COMPL_DICT[base]
# end def _get_compl_base


def _rc(seq: str) -> str:
    # Function returns reverse-complement "comrade" of passed DNA sequence.
    return ''.join(
        map(
         _get_compl_base,
         reversed(seq)
        )
    )
# end def _rc


def _fasta_generator(infpath: str) -> Generator[Tuple[str, str], None, None]:
    # Generator yields "fasta-tuples: 0-th element of such a tuple is sequence name,
    #   and 1-st element if sequence itself.

    curr_seq_name: str = '' # current sequence name
    curr_seq: str = '' # current sequence

    open_func: Callable[[str], ContextManager] # function for opening input file

    # Choose `open_func`
    if infpath.endswith('.gz'):
        open_func = partial(gzip.open, mode='rt', encoding='utf-8')
    else:
        open_func = partial(open, mode='rt', encoding='utf-8')
    # end if

    infile: Union[TextIO, BinaryIO]
    with open_func(infpath) as infile:

        eof: bool = False # indicates of End Of File is reached

        # Get the first sequence name
        curr_seq_name = infile.readline().strip()[1:]

        while not eof:
            # Get next line whatever it is
            line: str = infile.readline().strip()

            if line.startswith('>') or line == '':
                # We reached end of current sequence
                yield curr_seq_name, curr_seq # yield current sequence
                curr_seq_name = line[1:] # read new header
                curr_seq = '' # empty sequence
                if line == '': # no more sequences -- end of file
                    eof = True
                # end if
            else:
                curr_seq += line # new line is a sequence -- append it to `curr_seq`
            # end if
        # end while
    # end with
# end def _fasta_generator
