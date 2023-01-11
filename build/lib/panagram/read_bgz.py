import sys
import numpy as np

class KmerRef:
    def __init__(self, prefix):
        self.pac_fname = f"{prefix}.pac"
        ann_in = open(f"{prefix}.ann")

        self.ngenomes, self.nbits = map(int, ann_in.readline().split())
        self.nbytes = self.nbits // 8

        self.offsets = dict()
        for line in ann_in:
            nt_offs, name = line.split()[:2]
            self.offsets[name] = int(nt_offs)

        ann_in.close()

    def get_counts(self, name, start, end):
        byte_start = self.nbytes * (self.offsets[name] + start)
        length  = end - start

        pac = np.fromfile(self.pac_fname, offset=byte_start, count=length, dtype=f"int{self.nbits}")

        ret = np.zeros((len(pac), self.ngenomes), dtype=bool)
        for i in range(self.ngenomes):
            ret[:,i] = (pac >> i) & 1

        return ret

if __name__ == "__main__":
    ref = KmerRef(sys.argv[1])
    name, coords = sys.argv[2].split(":")
    start, end = map(int, coords.split("-"))

    print(ref.get_counts(name, start, end))
