import sys

def binary_search (sortlist, target):
    st, ed = 0, len(sortlist)-1
    while st <= ed:
        mid = int((st+ed) / 2)
        if sortlist[mid] == target:
            return mid
        elif sortlist[mid] > target:
            ed = mid - 1 
        elif sortlist[mid] < target:
            st = mid + 1
    return st


# simple hash function when binning size/step is constant
class bin_hash:

    def __init__(self,
                 ID_interval,
                 bin_size,
                 bin_step,
                 max_pos):

        self.ID_value = {}    
        self.bin_size = bin_size
        self.bin_step = bin_step
        self.max_pos = max_pos

        # map bin idx to bin ID
        self.idx_ID = {}
        self.ID_idx = {}
        for ID in ID_interval:
            st, ed = ID_interval[ID]
            assert st % self.bin_step == 0
            assert ed == st + self.bin_size
            idx = int(st / self.bin_step)
            assert idx not in self.idx_ID
            self.idx_ID[idx] = ID
            self.ID_idx[ID] = idx
            
        print ("hash function is built", file=sys.stderr)
        
    def find(self, pos):
        find_IDs = []
        idx = int(pos / self.bin_step)
        st = self.bin_step*idx
        ed = st + self.bin_size
        while pos >= st and pos < ed:
            try:
                ID = self.idx_ID[idx]
                find_IDs.append(ID)
            except:
                None
            idx -= 1
            if idx < 0:
                break
            st = self.bin_step*idx
            ed = st + self.bin_size
        return find_IDs

    def insert (self, pos, value):
        find_IDs = self.find(pos)
        for ID in find_IDs:
            if ID not in self.ID_value:
                self.ID_value[ID] = 0.0
            self.ID_value[ID] += value
        return find_IDs
        
    def find_range(self, rst, red):
        find_IDs = []

        idx = int(rst / self.bin_step)
        min_idx = idx
        st = self.bin_step*idx
        ed = st + self.bin_size
        while rst >= st and rst < ed:
            if idx < min_idx:
                min_idx = idx
            idx -= 1
            if idx < 0:
                break
            st = self.bin_step*idx
            ed = st + self.bin_size

        red = min(red, self.max_pos + 1)
        max_idx = int((red - 1) / self.bin_step)

        for idx in range(min_idx, max_idx+1):
            try:
                ID = self.idx_ID[idx]
                find_IDs.append(ID)
            except:
                None        
        return find_IDs

    def insert_range (self, rst, red, value):
        find_IDs = self.find_range(rst, red)
        for ID in find_IDs:
            idx = self.ID_idx[ID]
            st = self.bin_step*idx
            ed = st + self.bin_size
            a, b = max(st, rst), min(ed, red)
            length = b - a
            if ID not in self.ID_value:
                self.ID_value[ID] = 0.0
            self.ID_value[ID] += value*length
        return find_IDs

    def keys (self):
        return self.ID_value.keys()

    def values (self):
        return self.ID_value.values()

    def ID (self, id):
        return self.ID_value[id]
        
    def get (self):
        return self.ID_value

# build interval dictionary by using double hashing
class double_hash:
    def __init__(self,
                 ID_interval,
                 domain_size,
                 max_pos):
        
        self.ID_value = {}
        self.ID_interval = ID_interval

        edID = []
        for ID, interval in ID_interval.items():
            st, ed = interval
            edID.append([ed, ID])
        edID = sorted(edID)

        edlist, IDlist = [], []
        for ed, ID in edID:
            edlist.append(ed)
            IDlist.append(ID)
        
        self.domain_size = domain_size
        self.max_pos = max_pos
        self.domain_IDs = {}
        self.domain_num = max_pos // domain_size + 1

        for i in range(self.domain_num):
            self.domain_IDs[i] = []
            dst = i * self.domain_size
            ded = min(dst + self.domain_size, max_pos+1)
            idx1 = binary_search(edlist, dst)
            if idx1 == len(edlist):
                continue
            for j in range(idx1, len(edlist)):
                ID = IDlist[j]
                st, ed = self.ID_interval[ID]
                if st < ded:
                    self.domain_IDs[i].append(ID)
                
        print ("hash fucntion is built", file=sys.stderr)

    def __str__ (self):
        print ("%s\t%s\t%s\t%s" % ("ID", "st", "ed", "value"))
        for ID, value in self.ID_value.items():
            st, ed = ID_interval[ID]
            print ("%d\t%d\t%d\t%f" % (ID, st, ed, value))
        return

    def find (self, pos):
        find_IDs = []
        
        domain = pos // self.domain_size
        IDs = self.domain_IDs[domain]

        edID = []
        for ID in IDs:
            st, ed  = self.ID_interval[ID]
            edID.append([ed,ID])
        edID = sorted(edID)

        edlist, IDlist = [], []
        for ed, ID in edID:
            edlist.append(ed)
            IDlist.append(ID)
            
        idx = binary_search(edlist, pos)
        if idx == len(edlist):
            return find_IDs

        for i in range(idx, len(edlist)):
            ID = IDlist[i]
            st, ed = self.ID_interval[ID]
            if pos >= st and pos < ed:
                find_IDs.append(ID)
            
        return find_IDs

    def insert (self, pos, value):
        find_IDs = self.find(pos)
        for ID in find_IDs:
            if ID not in self.ID_value:
                self.ID_value[ID] = 0.0
            self.ID_value[ID] += value
        return find_IDs

    def find_range (self, rst, red):
        find_IDs = []
        domain1 = rst // self.domain_size
        domain2 = red // self.domain_size

        IDs = set([])
        for i in range(domain1, domain2 + 1):
            IDs |= set(self.domain_IDs[i])
        IDs = list(IDs)

        edID = []
        for ID in IDs:
            st, ed = self.ID_interval[ID]
            edID.append([ed, ID])
        edID = sorted(edID)

        edlist, IDlist = [], []
        for ed, ID in edID:
            edlist.append(ed)
            IDlist.append(ID)

        idx1 = binary_search(edlist, rst)
        if idx1 == len(edlist):
            return find_IDs
        
        for j in range(idx1, len(edlist)):
            ID = IDlist[j]
            st, ed = self.ID_interval[ID]
            if st < red:
                find_IDs.append(ID)

        return find_IDs
    
    def insert_range (self, rst, red, value):
        find_IDs = self.find_range(rst, red)
        for ID in find_IDs:
            st, ed = self.ID_interval[ID]
            a, b = max(st, rst), min(ed, red)
            length = b - a
            if ID not in self.ID_value:
                self.ID_value[ID] = 0.0
            self.ID_value[ID] += value*length
        return find_IDs

    def keys (self):
        return self.ID_value.keys()

    def values (self):
        return self.ID_value.values()

    def ID (self, id):
        return self.ID_value[id]
        
    def get (self):
        return self.ID_value
