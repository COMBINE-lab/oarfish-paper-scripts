#!/usr/bin/env python3

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


if __name__ == "__main__":

    import sys
    import gzip
    targets = sys.argv[1]
    txome = sys.argv[2]


    tgts = set([ l.rstrip() for l in open(targets, 'r').readlines()])

    def open_file(f):
        if f.endswith('.gz'):
            return gzip.open(f, 'rt')
        else:
            return open(f, 'r')

    with open_file(txome) as tfile:
        for name,seq,qual in readfq(tfile):
            if name in tgts:
                print(f">{name}")
                print(f"{seq}")
