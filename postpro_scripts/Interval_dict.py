import sys
import math

def binary_search (sortlist, target):
    st, ed = 0, len(sortlist)-1
    while st <= ed:
        mid = (st+ed) / 2
        if sortlist[mid] == target:
            return mid
        elif sortlist[mid] > target:
            ed = mid - 1 
        elif sortlist[mid] < target:
            st = mid + 1
    return st

def tuple_cmp (a,b):
    if a[0] < b[0]:
        return -1
    elif a[0] > b[0]:
        return 1
    else:
        if a[1] < b[1]:
            return -1
        elif a[1] > b[1]:
            return 1
        else:
            return 0

# simple hash function when binning size/step is constant
class bin_hash:

    def __init__(self,
                 bin_size,
                 bin_step,
                 max_pos=None,
                 ID_interval=None,
                 silent=False):

        self.ID_value = {}    
        self.bin_size = bin_size
        self.bin_step = bin_step

        # if ID:interval not defined
        if ID_interval == None:
            assert max_pos != None # max_pos must be provided
            ID_interval = {}
            for i in range(max_pos/bin_step + 1):
                ID = i
                st = bin_step * i
                ed = st + bin_size
                ID_interval[ID] = (st, ed)

        # map bin idx to bin ID
        self.idx_ID = {}
        self.ID_idx = {}
        for ID in ID_interval:
            st, ed = ID_interval[ID]
            assert st % self.bin_step == 0
            assert ed == st + self.bin_size
            idx = st / self.bin_step
            assert idx not in self.idx_ID
            self.idx_ID[idx] = ID
            self.ID_idx[ID] = idx

        # if max_pos not provided
        if max_pos == None:
            max_idx = max(self.idx_ID.keys())
            max_pos = max_idx * self.bin_step + self.bin_size
        self.max_pos = max_pos
        
        if not silent:
            print >> sys.stderr, "hash function is built"
        
    def find(self, pos):
        find_IDs = []
        idx = pos / self.bin_step
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

        idx = rst / self.bin_step
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
        max_idx = (red - 1) / self.bin_step

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

    def clear (self):
        self.ID_value = {}

# build interval dictionary by using double hashing
class double_hash:
    def __init__(self,
                 ID_interval,
                 domain_size=None,
                 max_pos=None,
                 silent=False):
        
        self.ID_value = {}
        self.ID_interval = ID_interval

        # record the boundaries of intervals
        edID = []
        for ID, interval in ID_interval.items():
            st, ed = interval
            edID.append([ed, ID])
        edID = sorted(edID, cmp=tuple_cmp)

        edlist, IDlist = [], []
        for ed, ID in edID:
            edlist.append(ed)
            IDlist.append(ID)

        # if max_pos not provided
        if max_pos == None:
            max_pos = edlist[-1]
        self.max_pos = max_pos

        # if domain_size not provided:
        if domain_size == None:
            domain_size = 10**(int(math.log10(max_pos))/2)
        self.domain_size = domain_size

        # categorize the IDs into each domains
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

        if not silent:
            print >> sys.stdout, "hash fucntion is built"

    def __str__ (self):
        print "%s\t%s\t%s\t%s" % ("ID", "st", "ed", "value")
        for ID, value in self.ID_value.items():
            st, ed = ID_interval[ID]
            print "%d\t%d\t%d\t%f" % (ID, st, ed, value)
        return

    def find (self, pos):
        find_IDs = []
        
        domain = pos // self.domain_size
        try:
            IDs = self.domain_IDs[domain]
        except:
            IDs = []

        edID = []
        for ID in IDs:
            st, ed  = self.ID_interval[ID]
            edID.append([ed,ID])
        edID = sorted(edID, cmp=tuple_cmp)

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
            try:
                IDs |= set(self.domain_IDs[i])
            except:
                pass
        IDs = list(IDs)

        edID = []
        for ID in IDs:
            st, ed = self.ID_interval[ID]
            edID.append([ed, ID])
        edID = sorted(edID, cmp=tuple_cmp)

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

    def clear (self):
        self.ID_value = {}
