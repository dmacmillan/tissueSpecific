import sys

class Matrix:

    def __init__(self,m,n,rownames=None,colnames=None):
        self.mtx = [[0]*n for x in xrange(m)]
        self.m = m
        self.n = n
        self.rownames = rownames
        self.colnames = colnames

    def __str__(self):
        if self.rownames is not None and self.colnames is not None:
            s = ['\t' + ('\t').join(self.colnames)]
        elif self.colnames is not None:
            s = [('\t').join(self.colnames)]
        for i,row in enumerate(self.mtx):
            if self.rownames:
                rowname = self.rownames[i]
                s.append(rowname + '\t' + ('\t').join([str(x) for x in row]))
            else:
                s.append(('\t').join([str(x) for x in row]))
        return ('\n').join(s)

    def col(self,c):
        return [x[c] for x in self.mtx]

    def row(self,r):
        return self.mtx[r]

    @staticmethod
    def readMatrix(_file, has_rownames=True, has_colnames=True, delim='\t'):
        mtx = []
        m = n = 0
        with open(_file, 'r') as f:
            if has_colnames:
                colnames = []
            if has_rownames:
                rownames = f.readline().strip().split(delim)
            for line in f:
                m += 1
                line = line.strip().split(delim)
                colnames.append(line[0])
                try:
                    mtx.append([float(x) for x in line[1:]])
                except ValueError:
                    sys.exit(line)
                if not n:
                    n = len(mtx[0])
        m = Matrix(m,n,rownames,colnames)
        m.mtx = mtx
        return m
