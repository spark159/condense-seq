import numpy as np
import analysis
#import maxcal_Dixit as maxcal
from sklearn.neighbors import KernelDensity

"""
def read_DNAshape(fname, id_seq):
    names = ['MGW', 'HelT', 'ProT', 'Roll']
    dic_list = [{} for i in range(len(names))]
    for i in range(len(names)):
        data = []
        for line in open(fname+"."+names[i]):
            line = line.strip()
            if line.startswith('>'):
                if data:
                    assert seq not in dic_list[i]
                    dic_list[i][seq] = data
                id = int(line[1:].strip())
                seq = id_seq[id]
                data =[]
                continue
            if line:
                temp = line.split(',')
                for k in range(len(temp)):
                    try:
                        temp[k] = float(temp[k])
                    except Exception:
                        pass
                data += temp
        assert seq not in dic_list[i]
        dic_list[i][seq] = data
    return dic_list
#seq_MGW, seq_HelT, seq_ProT, seq_Roll = read_DNAshape('phpoXp6Za', analysis.all_path(16))    
"""
class Slider:
    
    def __init__(self,
                id,
                ref_length,
                dyad_axis,
                left_offset,
                right_offset,
                seq,
                dyadmap,
                left_cutmap,
                right_cutmap,
                MGW,
                HelT,
                ProT,
                Roll):
        
        self.id = id
        self.ref_length = ref_length
        self.dyad_axis = dyad_axis
        self.left_offset = left_offset
        self.right_offset = right_offset
        self.seq = seq
        self.dyadmap = dyadmap
        self.left_cutmap = left_cutmap
        self.right_cutmap = right_cutmap
        self.MGW = MGW
        self.HelT = HelT
        self.ProT = ProT
        self.Roll = Roll

        
    def __add__ (self, other, norm_choice=True):
        assert self.ref_length == other.ref_length
        assert self.dyad_axis == other.dyad_axis
        assert self.left_offset == other.left_offset
        assert self.right_offset == other.right_offset
        if self.id == other.id:
            new_id = self.id
        else:
            new_id = (self.id, other.id)
        if self.seq == other.seq:
            new_seq = self.seq
        else:
            new_seq = (self.seq, other.seq)

        if norm_choice:
            dyadmap1 = analysis.norm(self.dyadmap)
            dyadmap2 = analysis.norm(other.dyadmap)
            left_cutmap1 = analysis.norm(self.left_cutmap)
            left_cutmap2 = analysis.norm(other.left_cutmap)
            right_cutmap1 = analysis.norm(self.right_cutmap)
            right_cutmap2 = analysis.norm(other.right_cutmap)
            
        new_dyadmap = [dyadmap1[i] + dyadmap2[i] for i in range(self.ref_length)]
        new_left_cutmap = [left_cutmap1[i] + left_cutmap2[i] for i in range(self.ref_length)]
        new_right_cutmap = [right_cutmap1[i] + right_cutmap2[i] for i in range(self.ref_length)]
        return Slider(new_id,
                      self.ref_length,
                      self.dyad_axis,
                      self.left_offset,
                      self.right_offset,
                      new_seq,
                      new_dyadmap,
                      new_left_cutmap,
                      new_right_cutmap)

    def __sub__ (self, other, norm_choice=True):
        assert self.ref_length == other.ref_length
        assert self.dyad_axis == other.dyad_axis
        assert self.left_offset == other.left_offset
        assert self.right_offset == other.right_offset
        if self.id == other.id:
            new_id = self.id
        else:
            new_id = (self.id, other.id)
        if self.seq == other.seq:
            new_seq = self.seq
        else:
            new_seq = (self.seq, other.seq)

        if norm_choice:
            #dyadmap1 = analysis.norm(self.dyadmap)
            #dyadmap2 = analysis.norm(other.dyadmap)
            dyadmap1 = self.KDE(band_width=0.5)
            dyadmap2 = other.KDE(band_width=0.5)
            left_cutmap1 = analysis.norm(self.left_cutmap)
            left_cutmap2 = analysis.norm(other.left_cutmap)
            right_cutmap1 = analysis.norm(self.right_cutmap)
            right_cutmap2 = analysis.norm(other.right_cutmap)
            
        new_dyadmap = [dyadmap1[i] / dyadmap2[i] for i in range(self.ref_length)]
        new_left_cutmap = [left_cutmap1[i] - left_cutmap2[i] for i in range(self.ref_length)]
        new_right_cutmap = [right_cutmap1[i] - right_cutmap2[i] for i in range(self.ref_length)]
        return Slider(new_id,
                      self.ref_length,
                      self.dyad_axis,
                      self.left_offset,
                      self.right_offset,
                      new_seq,
                      new_dyadmap,
                      new_left_cutmap,
                      new_right_cutmap)
        
    def get_dyadmap(self):
        return self.dyadmap

    def GC_content(self):
        num=0.0
        for nt in self.seq:
            if nt in 'GC':
                num+=1
        return (num/float(len(self.seq)))*100

    def Amer_len_detail(self, pos=True):
        seq = self.seq
        num = []
        num_pos = {}
        i = 0
        while i < len(seq):
            if seq[i] in 'AT':
                nt = seq[i]
                count = 1
                j = i + 1
                while j < len(seq):
                    if seq[j] != nt:
                        break
                    count +=1
                    j +=1
                num.append(count)
                if count not in num_pos:
                    num_pos[count] = []
                num_pos[count].append(i)
                i = j
            else:
                i +=1
        if pos:
            return num_pos
        if len(num) == 0:
            return 0
        return max(num)
    
    def Amer_len(self):
        #win,st = self.id.split('-')
        loc, mtype, nts = self.id.split('-')
        return len(nts)
    
    def weighted_GC (self):
        result = 0.0
        freq = analysis.norm(self.dyadmap)
        seq = self.seq
        for i in range(147/2, len(freq)-147/2):
            NCPseq = seq[i-147/2:i+147/2 + 1]
            GC = analysis.GC_content(NCPseq)
            result += GC*freq[i]
        return result

    def weighted_Amer_len (self):
        result = 0.0
        freq = analysis.norm(self.dyadmap)
        seq = self.seq
        for i in range(147/2, len(freq)-147/2):
            NCPseq = seq[i-147/2:i+147/2 + 1]
            num_pos = analysis.Amer_len_detail(NCPseq)
            total = 0.0
            mean_num = 0.0
            for num, pos in num_pos.items():
                total += len(pos)
                mean_num += len(pos)*num
            mean_num = mean_num / total
            result += mean_num*freq[i]
        return result

    def MGW_profile(self):
        return seq_MGW[self.seq]

    def HelT_profile(self):
        return seq_HelT[self.seq]

    def ProT_profile(self):
        return seq_ProT[self.seq]

    def Roll_profile(self):
        return seq_Roll[self.seq]

    def MGW(self):
        return np.mean(MGW_profile[self.seq])

    def HelT(self):
        return np.mean(HelT_profile[self.seq])

    def ProT(self):
        return np.mean(ProT_profile[self.seq])

    def Roll(self):
        return np.mean(Roll_profile[self.seq])

    def read_counts (self, choice):
        if choice == 'dyad':
            map = self.dyadmap
        elif choice == 'left':
            map = self.left_cutmap
        elif choice == 'right':
            map = self.right_cutmap
        total = 0.0;
        for i in range(len(map)):
            total += map[i]
        return total

    def rel_bin_counts (self, choice, bin_num):
        if choice == 'dyad':
            map = self.dyadmap
        elif choice == 'left':
            map = self.left_cutmap
        elif choice == 'right':
            map = self.right_cutmap
        assert bin_num >= 0 and bin_num < len(map)
        return map[bin_num] / self.read_counts(choice)

    def find_peaks(self, choice, num, back= False):
        def key_cmp(a, b):
            if a[0] <= b[0]:
                return -1
            else:
                return 1

        def value_cmp(a, b):
            if a[1] <= b[1]:
                return -1
            else:
                return 1
    
        def sub_background (map, frac=0.1):
            thres= min(map) + frac*(max(map)-min(map))
            new = [0 for i in range(len(map))]
            for i in range(len(map)):
                if map[i] > thres:
                    new[i] = map[i]
            return new
        if choice == 'dyad':
            map = self.dyadmap
        elif choice == 'left':
            map = self.left_cutmap
        elif choice == 'right':
            map = self.right_cutmap
        nmap = map
        if back:
            nmap = sub_background(map)
        peak_sig={}
        for i in range(1, len(nmap)-1):
            if nmap[i] > nmap[i-1] and nmap[i] > nmap[i+1]:
                peak_sig[i] = nmap[i]

        peak_sig=[[peak, sig]  for peak, sig in peak_sig.items()]
        peak_sig=sorted(peak_sig, cmp=value_cmp, reverse=True)
        peaks=[]
        for i in range(min(len(peak_sig), num)):
            peaks.append(peak_sig[i])
        return peaks

        
    def median_pos (self, scale=1.0, expos=None, selected=None):
        data = []
        if selected == None:
            selected = range(len(self.dyadmap))
        for i in selected:
            if i == expos:
                continue
            for k in range(int(scale*self.dyadmap[i])):
                data.append(i-self.dyad_axis)
        return np.median(data)

    
    def mean_pos (self, expos=None, selected=None):
        mean_pos, total = 0.0, 0.0
        if selected == None:
            selected = range(len(self.dyadmap))
        for i in selected:
            if i == expos:
                continue
            mean_pos += (i-self.dyad_axis)*self.dyadmap[i]
            total += self.dyadmap[i]
        return mean_pos/total

    def max_dis (self, peak_num = 5):
        def find_max (L):
            index = 0
            for i in range(len(L)):
                if L[index] < L[i]:
                    index = i
            return index, L[index]
        dyadpeak = analysis.find_peaks(self.dyadmap, num = peak_num)
        abspeak = [abs(peak - self.dyad_axis) for peak, sig in dyadpeak]
        idx, _ = find_max(abspeak)
        mdis = dyadpeak[idx][0] - self.dyad_axis
        return mdis

    def minor_peak_sig (self, minor_loc):
        if choice == 'dyad':
            map = self.dyadmap
        elif choice == 'left':
            map = self.left_cutmap
        elif choice == 'right':
            map = self.right_cutmap
        map_peaks = analysis.find_peaks(map, num=2)
        if len(map_peaks) < 2:
            return 0.0
        loc, sig = map_peaks[1][0] - map_peaks[0][0], map_peaks[1][1]/map_peaks[0][1]
        if loc != minor_loc:
            return 0.0
        return sig

    def peak_signal (self, choice='dyad', num=15, broden=1):
        if choice == 'dyad':
            map = self.dyadmap
        elif choice == 'left':
            map = self.left_cutmap
        elif choice == 'right':
            map = self.right_cutmap
        loc_sig_list = self.find_peaks(choice=choice, num=num, back= False)
        peak_signal = [0.0]*self.ref_length
        for loc_sig in loc_sig_list:
            loc, sig = loc_sig
            peak_signal[loc] += sig
            for offset in range(1, broden+1):
                peak_signal[loc-offset] += sig
                peak_signal[loc+offset] += sig
        return peak_signal

    def KDE (self, scale=100.0, band_width=1.0):
        X, X_plot = [], []
        for k in range(len(self.dyadmap)):
            for num in range(int(scale*self.dyadmap[k])):
                X.append([k])
            X_plot.append([k])
        #X_plot = np.linspace(0,len(self.dyadmap)-1, num = len(self.dyadmap)*5)[:,np.newaxis]
        kde = KernelDensity(kernel="gaussian", bandwidth=band_width).fit(X)
        log_density = kde.score_samples(X_plot)
        return np.exp(log_density)

    def entropy (self, scale=1.0, band_width=0.5):
        X, X_plot = [], []
        for k in range(len(self.dyadmap)):
            for num in range(int(scale*self.dyadmap[k])):
                X.append([k])
            X_plot.append([k])
        #X_plot = np.linspace(0,len(self.dyadmap)-1, num = len(self.dyadmap)*5)[:,np.newaxis]
        kde = KernelDensity(kernel="gaussian", bandwidth=band_width).fit(X)
        log_density = kde.score_samples(X_plot)
        kde = np.exp(log_density)
        entropy = 0.0
        for i in range(len(log_density)):
            entropy += -kde[i]*log_density[i]
        return entropy

    def energy_profile (self, kT=1, scale=1.0, band_width=0.5):
        X, X_plot = [], []
        for k in range(len(self.dyadmap)):
            for num in range(int(scale*self.dyadmap[k])):
                X.append([k])
            X_plot.append([k])
        #X_plot = np.linspace(0,len(self.dyadmap)-1, num = len(self.dyadmap)*5)[:,np.newaxis]
        kde = KernelDensity(kernel="gaussian", bandwidth=band_width).fit(X)
        log_density = kde.score_samples(X_plot)
        return -kT * log_density

    def force_profile (self, kT=1, band_width=0.5):
        X = []
        for k in range(len(self.dyadmap)):
            for num in range(int(self.dyadmap[k])):
                X.append([k])
        #X_plot = np.linspace(0,len(self.dyadmap), num = (len(self.dyadmap)+1)*5)[:,np.newaxis]
        X_plot = np.arange(0,len(self.dyadmap), 0.1)[:,np.newaxis]
        kde = KernelDensity(kernel="gaussian", bandwidth=band_width).fit(X)
        energy_profile = -kT* kde.score_samples(X_plot)    
        X_axis = list(X_plot[:,0])
        force_profile = []
        for i in range(len(self.dyadmap)):
            idx = X_axis.index(i)
            dX = X_axis[idx+1] - X_axis[idx]
            denergy = energy_profile[idx+1] - energy_profile[idx]
            force = - denergy / dX
            force_profile.append(force)
        return np.asarray(force_profile)

    def eqm_maxcal(self):
        kde = self.KDE(band_width=0.5)
        kde = kde[147/2:len(kde)-147/2]
        total = sum(kde)
        eqm_dist = [float(value)/total for value in kde]
        N = len(eqm_dist)
        edge_matrix = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                if abs(i-j) > 1:
                    continue
                edge_matrix[i][j] = 1.0
        eqm_matrix = maxcal.eqm_maxcal(eqm_dist, edge_matrix=edge_matrix)
        return eqm_matrix

    def eqm_flux (self):
        eqm_matrix = self.eqm_maxcal()
        N = len(eqm_matrix)
        net_flow = []
        for i in range(N):
            if i == 0:
                net_flow.append(eqm_matrix[i][i+1])
            elif i == N-1:
                net_flow.append(-eqm_matrix[i][i-1])
            else:
                net_flow.append(eqm_matrix[i][i+1] - eqm_matrix[i][i-1])
        return net_flow
