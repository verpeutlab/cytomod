# Cytomod #

by Coby Viner <cviner at cs dot toronto dot edu>

 All scripts contain a detailed description, explaining their purpose and usage. [`cytomod.py`](src/cytomod.py) is the main script.

 Cytomod itself can be used by providing it with an unmodified genome assembly and track datasets indicating which nucleobases to modify.
 This is described in detail from its Usage documentation, which we also provide below:

```
 usage: cytomod.py [-h]
                  (-G GENOMEDATAARCHIVEFULLNAME | -d GENOME_DIR TRACKS_DIR)
                  [--archiveOutDir ARCHIVEOUTDIR]
                  [--archiveOutName ARCHIVEOUTNAME] [-r REGION]
                  [-c [CENTEREDREGION]] [-R [RANDOMREGION]] [-A {m,M,u,l}]
                  [-p PRIORITY] [-b | -B] [--BEDOutDir BEDOUTDIR]
                  [-f [FASTAFILE]] [-I] [--mh] [--fC] [--fc]
                  [-M [MASKREGIONS]] [--maskAllUnsetRegions] [-v] [-V]

optional arguments:
  -h, --help            show this help message and exit
  -G GENOMEDATAARCHIVEFULLNAME, --genomeDataArchiveFullname GENOMEDATAARCHIVEFULLNAME
                        The genome data archive. It must contain all needed
                        sequence and track files. If one is not yet created,
                        use '-d' instead to create it.
  -d GENOME_DIR TRACKS_DIR, --archiveCompDirs GENOME_DIR TRACKS_DIR
                        Two arguments first specifying the directory
                        containing the genome and then the directory
                        containing all modified base tracks. The genome
                        directory must contain (optionally Gzipped) FASTA
                        files of chromosomes and/or scaffolds. The sequence
                        files must end in either a ".fa" or ".mask" extension
                        (with an optional ".gz" suffix). The track directory
                        must contain (optionally Gzipped) genome tracks. They
                        must have an extension describing their format. We
                        currently support: ".wig", ".bed", and ".bedGraph".
                        Provided BED files must have exactly four columns. The
                        fourth column must be numeric. Any rows with a data
                        value of zero will be ignored (this does not apply to
                        masking; see '-M'), except that such positions will be
                        considered as having evidence against that particular
                        modification (e.g. will be considered as 'set' and
                        will therefore not be masked by '--
                        maskAllUnsetRegions'). The filename of each track must
                        specify what modified nucleobase it pertains to; one
                        of: ['5mC', '5hmC', '5fC', '5caC']. Track names can
                        also be of ambiguity codes (e.g. "5xC") on the
                        positive strand only. Such tracks directly specify
                        ambiguous loci. If multiple tracks of the same type
                        are provided, all such tracks will be added to the
                        archive. The output sequence will default to the union
                        of all of the same modification type (but see '-I').
                        Alternatively, the track name can contain MASK, in
                        which case masking can be used via '-M' (refer to that
                        option for details). Instead of a track directory, a
                        single filename that meets the aforementioned
                        requirements may be provided if the archive is to
                        contain only one track. Ensure that all tracks are
                        mapped to the same assembly and that this assembly
                        matches the genome provided. This will create a genome
                        data archive in an "archive" sub-directory of the
                        provided track directory. Use '-G' instead to use an
                        existing archive.
  --archiveOutDir ARCHIVEOUTDIR
                        Only applicable if '-d' is used. The directory in
                        which to save the created genome data archive. If not
                        specified, this defaults to the directory containing
                        the tracks (i.e. the second argument provided to
                        '-d'). This defaults to archive.
  --archiveOutName ARCHIVEOUTNAME
                        Only applicable if '-d' is used. The name of the
                        archive (i.e. the name of the directory which
                        comprises the genome data archive).
  -r REGION, --region REGION
                        Only output the modified genome for the given region.
                        This can either be via a file or a region
                        specification string. In the latter case, the region
                        must be specified in the format:
                        chrm<ID>:<start>-<end> (ex. chr1:500-510). If a file
                        is being provided, it can be in any BEDTools-supported
                        file format (BED, VCF, GFF, and Gzipped versions
                        thereof). The full path to the file should be provided
                        (or just the file name for the current directory).
  -c [CENTEREDREGION], --centeredRegion [CENTEREDREGION]
                        If used in conjunction with '-r', only output the
                        modified genome for the given base pair interval
                        (defaults to 500 bp), centered around the given
                        region. NB: This region does not necessarily
                        correspond to the centre of the peak, since the
                        region's start and end coordinates alone are used to
                        find the centre, as opposed to any peak information
                        (from a narrowPeaks file, for example).
  -R [RANDOMREGION], --randomRegion [RANDOMREGION]
                        Output the modified genome for a random region. The
                        chrmomsome will be randomly selected and its
                        coordinate space will be randomly and uniformly
                        sampled. A length for the random region can either be
                        specified or it will otherwise be set to a reasonably
                        small default. The length chosen may constrain the
                        selection of a chromosome.
  -A {m,M,u,l}, --alterIncludedChromosomes {m,M,u,l}
                        Include or exclude chromosome types. 'u': Use only
                        autosomal chromosomes (excludes chrM). 'l': Use only
                        allosomal chromosomes (excludes chrM). 'm': Use only
                        the mitochondrial chromosome. 'M': Include the
                        mitochondrial chromosome. NB: default chromosomal
                        exclusions include: unmapped data, haplotypes, and
                        chrM. This parameter will be ignored if a specific
                        genomic region is queried via '-r', but will be
                        considered if a file of genomic regions is provided
                        (also via '-r').
  -p PRIORITY, --priority PRIORITY
                        Specify the priority of modified bases. The default
                        is:fhmcwxyz, which is based upon the resolution of the
                        biological protocol (i.e. single-base > any chemical >
                        any DIP).
  -b, --suppressBED     Do not generate any BED tracks.
  -B, --onlyBED         Only generate any BED tracks (i.e. do not output any
                        sequence information). Note that generated BED files
                        are always appended to and created in the CWD
                        irrespective of the use of this option. This parameter
                        is ignored if '-f' is used.
  --BEDOutDir BEDOUTDIR
                        Only applicable if '-b' is not used. The directory in
                        which to save the created BED tracks. If not
                        specified, this defaults to the current working
                        directory.
  -f [FASTAFILE], --fastaFile [FASTAFILE]
                        Output to a file instead of STDOUT. Provide a full
                        path to a file to append the modified genome in FASTA
                        format. If this parameter is invoked without any
                        arguments, a default filename will be used within the
                        current directory. This will override the '-B'
                        parameter (i.e. a FASTA file with always be produced).
                        The output file will be Gzipped iff the path provided
                        ends in ".gz".
  -I, --intersection    If multiple files of the same modification type are
                        given, take their intersection. This option is used to
                        override the default, which is to take their union.
  -v, --verbose         increase output verbosity
  -V, --version         show program's version number and exit

Ambiguous Modification:
  Specify that some of the data provided for a given modified base is unable
  to differentiate between some number of modifications. This ensures that
  Cytomod outputs the correct ambiguity code such that modified genomes do
  not purport to convey greater information than they truly contain. The
  most general applicable ambiguities should be specified. Therefore, each
  modified nucleobase may reside in at most one specified set of
  ambiguities.

  --mh                  Specify that input data is not able to differentiate
                        between 5mC and 5hmC. This would be the case if the
                        data originated from a protocol which only included
                        conventional bisulfite sequencing.
  --fC                  Specify that input data is not able to differentiate
                        between 5fC and C. This would be the case if the data
                        originated from a protocol which only included
                        oxidative bisulfite sequencing.
  --fc                  Specify that input data is not able to differentiate
                        between 5fC and 5caC. This would be the case if the
                        data originated from a protocol which only included
                        M.SssI methylase-assisted bisulfite sequencing.
  -M [MASKREGIONS], --maskRegions [MASKREGIONS]
                        Hard mask C/G nucleobases to unknown state. Assumes
                        that the archive contains or is being built with a
                        trackname containing "MASK". The containing loci will
                        be interpreted as nucleobases of unknown modification
                        state. They will be accordingly set to the appropriate
                        (maximally) ambiguous base. This will override any
                        other modifications at those loci. This parameter can
                        accept an optional argument, indicating a value at and
                        below which the locus is considered ambiguous. If not
                        provided, this defaults to 0. An example use case for
                        this option would be to use a mask file containing
                        coverage information and to mask all bases of
                        insufficient coverage.
  --maskAllUnsetRegions
                        Hard mask all C/G nucleobases without any modification
                        information to unknown state. Masked nucleobases are
                        those lacking data, that is, bases not present in the
                        archive nor in any files used to generate an archive.
                        Set, but unmodified, bases can be provided, within any
                        modified genomic interval file with a value of 0. See
                        '-d' for further details. Unset modifiable bases will
                        be accordingly set to the appropriate (maximally)
                        ambiguous base. This will override any other
                        modifications at those loci. An example use case for
                        this option would be when using array data, for which
                        only a subset of bases are queried.
```
