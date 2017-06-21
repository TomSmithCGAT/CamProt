####################### ------------------------- ##############################
# code below 'borrowed' from CGAT.FastaIterator (trying to limit dependencies)
class FastaRecord:
    """a :term:`fasta` record.

    Attributes
    ----------
    title: string
       the title of the sequence

    sequence: string
       the sequence

    fold : int
       the number of bases per line when writing out
    """

    def __init__(self, title, sequence, fold=False):

        self.title = title
        self.sequence = sequence
        self.fold = fold

    def __str__(self):
        ''' str method for writing out'''

        if self.fold:
            seq = [self.sequence[i:i + self.fold]
                   for i in range(0, len(self.sequence), self.fold)]
        else:
            seq = (self.sequence,)

        return ">%s\n%s" % (self.title, "\n".join(seq))


class FastaIterator:
    '''a iterator of :term:`fasta` formatted files.

    Yields
    ------
    FastaRecord

    '''

    def __init__(self, f, *args, **kwargs):
        self.iterator = iterate(f)

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.iterator)

    def next(self):
        return next(self.iterator)


def iterate(infile, comment="#", fold=False):
    '''iterate over fasta data in infile

    Lines before the first fasta record are
    ignored (starting with ``>``) as well as
    lines starting with the comment character.

    Parameters
    ----------
    infile : File
        the input file
    comment : char
        comment character
    fold : int
        the number of bases before line split when writing out

    Yields
    ------
    FastaRecord
    '''

    h = infile.readline()[:-1]

    if not h:
        raise StopIteration

    # skip everything until first fasta entry starts
    while h[0] != ">":
        h = infile.readline()[:-1]
        if not h:
            raise StopIteration
        continue

    h = h[1:]
    seq = []

    for line in infile:

        if line.startswith(comment):
            continue

        if line.startswith('>'):
            yield FastaRecord(h, ''.join(seq), fold)

            h = line[1:-1]
            seq = []
            continue

        seq.append(line[:-1])

    yield FastaRecord(h, ''.join(seq), fold)
####################### ------------------------- ##############################
