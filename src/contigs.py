# -*- encoding: utf-8 -*-

import re
import gzip
from functools import partial
from typing import List, Tuple, Dict, Generator, NewType
from typing import Callable, ContextManager, TextIO

from src.platform import platf_depend_exit


# Dictionary maps complementary bases according to IUPAC:
_COMPL_DICT: Dict[str, str] = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
    'H': 'D', 'V': 'B', 'U': 'A', 'N': 'N'
}

# All possible bases for fasta validation
_INVALID_SEQ_PATTERN = r'[^AGTCRYSWKMBDHVUN]+'

# Pattern that matches SPAdes's header of a seqeunce in FASTA file
_SPADES_PATTERN: str = r'^NODE_[0-9]+_length_[0-9]+_cov_[0-9,\.]+'


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
    # :param infpath: path to input fasta file;
    # :param maxk: maximum k-mer length to consider;

    # Initialize result colletion
    contig_collection: ContigCollection = list()

    # Iterate over contigs and form contig_collection
    contig_header: str
    contig_seq: str
    for contig_header, contig_seq in _fasta_generator(infpath):

        # Simplify name
        contig_name: str = _format_contig_name(contig_header)

        # Parse coverage
        cov: float = _parse_coverage(contig_header)

        # Calculate GC-content
        gc_content: float = _calc_gc_content(contig_seq)

        # Append recently created contig to the `contig_collection`
        contig_collection.append(
            Contig(
                name=contig_name,
                length=len(contig_seq),
                cov=cov,
                gc_content=gc_content,
                start=contig_seq[:maxk],
                rcstart=_rc(contig_seq[:maxk].upper()),
                end=contig_seq[-maxk:].upper(),
                rcend=_rc(contig_seq[-maxk:].upper())
            )
        )
    # end for

    return contig_collection
# end def get_contig_collection


def _calc_gc_content(sequence: str) -> float:
    # Function calculates GC-content.
    # :param sequence: sequence to calculate GC content;

    gc_count: int = 0

    base: str
    for base in ('G', 'C', 'S'):
        gc_count += sequence.count(base)
    # end for

    gc_content: float = (gc_count / len(sequence) * 100)

    return gc_content
# end def _calc_gc_content


def _is_spades_header(header: str) -> bool:
    # function determines whether header is created by SPAdes.
    # :param header: fasta header of interest;
    return not re.search(_SPADES_PATTERN, header) is None
# end def _is_spades_header


def _format_contig_name(fasta_header: str) -> str:
    # Function simplifies fasta header.
    # :param fasta_header: fasta header of interest;

    # Remove extra underscores from the name
    fasta_header = fasta_header.strip('_')

    contig_name: str
    if _is_spades_header(fasta_header):
        # Get name in 'NODE_<NUMBER>' format
        contig_name = 'NODE_' + fasta_header.split('_')[1]
    else:
        contig_name = fasta_header # use full header as name
    # end if

    return contig_name
# end def _format_contig_name

def _parse_coverage(fasta_header: str) -> str:
    # Function parses coverage information from fasta header.
    # `combinator-FQ` can parse SPAdes and A5 assemblies correctly,
    #   and only SPAdes specifies coverage in fasta header.
    # Therefore we'll just assign None to coverage ,
    #   if contigs were not assembled with SPAdes.
    # :param fasta_header: fasta header of interest;

    # Remove extra underscores from the name
    fasta_header = fasta_header.strip('_')

    cov: float
    if _is_spades_header(fasta_header):
        # Parse coverage from SPAdes's header
        try:
            cov = float(fasta_header.split('_')[5])
        except ValueError:
            cov = None
        # end try
    else:
        cov = None
    # end if

    return cov
# end def _parse_coverage


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

    infile: TextIO
    with open_func(infpath) as infile:

        eof: bool = False # indicates of End Of File is reached

        # Get the first sequence name
        curr_seq_name = infile.readline().strip()

        while not eof:
            # Get next line whatever it is
            line: str = infile.readline().strip()

            if line.startswith('>') or line == '':
                # We reached end of current sequence

                # Validate parsed sequence
                try:
                    _validate_fasta(curr_seq_name, curr_seq)
                except ValueError:
                    platf_depend_exit(1)
                # end try

                yield curr_seq_name[1:], curr_seq # yield current sequence

                curr_seq_name = line # read next header
                curr_seq = ''        # empty sequence

                if line == '': # no more sequences -- end of file
                    eof = True
                # end if
            else:
                curr_seq += line.upper() # new line is a sequence -- append it to `curr_seq`
            # end if
        # end while
    # end with
# end def _fasta_generator


def _validate_fasta(curr_seq_name: str, curr_seq: str):
    # Function validates fasta given record.
    # :param curr_seq_name: fasta header of current sequence (WITH preceding `>`);
    # :param curr_seq: sequence casted to uppercase;

    # Validate header
    if not curr_seq_name.startswith('>'):
        print('Error: current file is not a fasta file:')
        print('  (putative) header does not start with `>` character.')
        print('Here is this (putative) header: `{}`.'.format(curr_seq_name))
        raise ValueError
    # end if

    # Check if `curr_seq` is not empty string
    if curr_seq == '':
        print('Error: sequence `{}` is empty.'.format(curr_seq_name))
        raise ValueError
    # end if

    # Validate sequence
    invalid_bases: List[str] = re.findall(_INVALID_SEQ_PATTERN, curr_seq)
    if len(invalid_bases) != 0:
        print('Error: invalid bases found in fasta sequence `{}`.'.format(curr_seq_name))
        # Surround invalid data with backticks for convenience
        print('Here they are: {}'.format(
            ', '.join(map(lambda x: '`{}`'.format(x), invalid_bases))
        ))
        raise ValueError
    # end if
# end def _validate_fasta
