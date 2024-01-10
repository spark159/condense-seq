class Matrix (object):
    def __init__ (self, m, n, value=0.0):
        self.m=m; self.n=n;
        self.matrix=[]
        for i in range(m):
            self.matrix.append([value]*n)
    
    def __getitem__ (self, idx):
        return self.matrix[idx]
    
    def __setitem__ (self, idx, value):
        self.matrix[idx]=value
    
    def getrank (self):
        return (self.m,self.n)
        
    def __add__ (self, other):
        if type(other) == Matrix and self.getrank() == other.getrank():
            new_one = Matrix(self.m,self.n)
            for i in range(self.m):
                for j in range(self.n):
                    new_one.matrix[i][j] = self.matrix[i][j] + other.matrix[i][j]
        elif type(other) == int or type(other) == float:
            new_one = Matrix(self.m,self.n)
            for i in range(self.m):
                for j in range(self.n):
                    new_one[i][j] = self.matrix[i][j] + other
        else:
            print "Not possible"
            return None
        return new_one
    
    def __sub__ (self, other):
        if type(other) == Matrix and self.getrank() == other.getrank():
            new_one = Matrix(self.m,self.n)
            for i in range(self.m):
                for j in range(self.n):
                    new_one.matrix[i][j] = self.matrix[i][j] - other.matrix[i][j]
        elif type(other) == int or type(other) == float:
            new_one = Matrix(self.m,self.n)
            for i in range(self.m):
                for j in range(self.n):
                    new_one[i][j] = self.matrix[i][j] - other
        else:
            print "Not possible"
            return None
        return new_one

    def __mul__ (self, other):
        if type(other) == Matrix and self.n == other.m:
            new_one = Matrix(self.m, other.n)
            for i in range(self.m):
                for j in range(other.n):
                    value = 0
                    for k in range(self.n):
                        value += self.matrix[i][k]*other.matrix[k][j]
                    new_one.matrix[i][j] = value
        elif type(other) == int or type(other) == float:
            new_one = Matrix(self.m,self.n)
            for i in range(self.m):
                for j in range(self.n):
                    new_one[i][j] = self.matrix[i][j]*other
        else:
            print "Not possible"
            return None
        return new_one
    
    def __str__ (self):
        s=''
        for i in range(self.m):
            s += str(self.matrix[i]) + '\n'
        return s[0:len(s)-1]
