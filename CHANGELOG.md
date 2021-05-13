# combinator-FQ changelog

## 2021-05-12 edition

- combinator-FQ now can calculate expected genome length and multiplicity of contigs based on number of overlaps, if no coverage information is present in sequence headers.
- Minor output formatting fixes.

### Version changes:

`1.4.a -> 1.5.a`

## 2021-03-26 edition

- performed total refactoring. Memory usage is reduced. Made some cosmetic fixes.
- Moved combinator-FQ to the separate (this) repository from [https://github.com/masikol/cager-misc](https://github.com/masikol/cager-misc).

### Version changes:

`1.3.f -> 1.4.a`

## 2021-01-21 edition

- bug fix (combinator used to ternimate in the end if it couldn't extract ordinal number of contig from contig's name);

### Version changes:

`1.3.e -> 1.3.f`

## 2020-12-03 edition

- `-o` option addded.

### Version changes:

`1.3.d -> 1.3.e`

## 2020-10-08 edition

- combinator-FQ now rounds coverage with 2 trailing digits.

### Version changes:

`1.3.c --> 1.3.d`

## 2020-06-15 edition.

- errorneous handling of SPADEs's `NODE_1` with zero coverage fixed;

### Version changes:

`1.3.b --> 1.3.c`

## 2020-04-24 edition.

- compatibility with a5 improved;

### Version changes:

`1.3.a --> 1.3.b`

## 2020-04-15 edition.

- calculating of expected genome length is improved;

### Version changes:

`1.2.b --> 1.3.a`;

## 2020-04-10 evening edition.

- a5-compatibility bug fixed;

### Version changes:

`1.2.a --> 1.2.b`

## 2020-04-10 edition.

- calculating of expected genome length is embedded once again -- not considering contigs multiplicity;
- fasta-GC-content is added;

### Version changes:

`1.1.c --> 1.2.a`

## 2020-03-18 edition.

- calculating of expected genome length is disabled because it's impossible to handle high-copy replicons properly;

### Version changes:

`1.1.b --> 1.1.c`

## 2020-03-09 edition.

- calculation of expected genome length fixed;

### Version changes:

`1.1.a --> 1.1.b`

## 2020-03-07 edition.

- output file separated into three: `adjacent_contigs`. `full_matching_log` and `summary`;
- A5-coverage summary bug fixed;

### Version changes:

`1.0.a --> 1.1.a`

## 2020-03-06 edition.

- combinator-FQ version `1.0.a` added;
