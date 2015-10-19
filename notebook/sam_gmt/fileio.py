"""
    file_io
    =======

    This submodule implements the functionality for interacting
    with file input and output.


    Links:
    ------
    - http://samtools.github.io/hts-specs/SAMv1.pdf
"""


def read_rawcounts(filename):
    """ Read in raw count files.
    """
    result = {}

    with open(filename, "r") as fh_in:
        for row in fh_in.readlines()[1:]:
            # Read columns from row.
            iden, t1, t11, t10, t101, leng = row.split()

            # Cast row values to propper types.
            value = list(map(int, (t1, t11, t10, t101))) + [float(leng)]

            # Assign dictionary value.
            result[iden] = value

    return result


def read_sam(filename):
    """ Reads a Sam file and returns a list of tuples representing rows.
        Mandatory fields are QNAME, FLAG, RNAME, POS, MAPQ, CIGAR,
            RNEXT, PNEXT, TLEN, SEQ and QUAL.
    """
    with open(filename, "r") as fh_in:
        # Create a list of rows, excluding ones starting with '@'.
        result = [tuple(row.strip().split('\t')) for row in fh_in
                  if not row.startswith("@")]

    return result


def write_rpkm(rpkm_data, filename):
    """ Write out rpkm (reads per kilobase per million)
        data to a file.
    """
    with open(filename, "w") as fh_out:
        for key, value in rpkm_data.iteritems():
            line = "{0}\t{1}\n".format(key, "\t".join(value))
            fh_out.write(line)
