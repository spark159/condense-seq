import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import random
import copy
import heapq
import bisect
import pickle
import cProfile, pstats, StringIO
import numpy as np
from numpy import linalg as LA
#import scipy.special as special
#from scipy.integrate import quad

def tuple_cmp (a, b):
    if a[0] < b[0]:
        return -1
    elif a[0] > b[0]:
        return 1
    else:
        if a[1] < b[1]:
            return -1
        return 1

# bisection search
def bisect (sL, target):
    st, ed = 0, len(sL)
    while st != ed:
        mid = (st + ed) / 2
        if target < sL[mid]:
            ed = mid
        elif target > sL[mid]:
            st = mid + 1
        elif target == sL[mid]:
            return mid
    return None

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


def half_right_normal (mu, sigma):
    value = random.gauss(mu, sigma)
    if value < mu:
        return 2*mu - value
    return value
        
class Graph (object):
    def __init__ (self,
                  nodes, # node set
                  edges, # edge set
                  node_value={}
                  edge_weight={}): 

        # combine all information into single data structure
        self.node_info = {}
        
        for node in nodes:
            try:
                value = node_value[node]
            except:
                value = None
            self.node_info[node] = {}
            self.node_info[node]['value'] = value # node value
            self.node_info[node]['next'] = {} # next node
            
        for edge in edges:
            node1, node2 = edge
            try:
                weight = edge_weight[node1][node2]
            except:
                weight = None
            self.node_info[node1]['next'][node2] = weight # edge weight
            
    def add_node (self, node, value=None):
        self.node_info[node] = {'value':value, 'next':{}}
        return

    def add_edge (self, node1, node2, weight=None):
        # two nodes already in the graph
        assert node1 in self.node_info
        assert node2 in self.node_info
        # update edge
        self.node_info[node1]['next'][node2] = weight
        return

    def remove_node (self, node):
        # delete the node and all edges from the node
        del self.node_info[node]
        # delete all edges toward the node
        for start in self.node_info:
            try:
                del self.node_info[start]['next'][node]
            except:
                pass
        return 

    def remove_edge (self, node1, node2):
        del self.node_info[node1]['next'][node2]
        return 

    # depth-first search
    def DFS (self, source, target=None):
        visited = []
        stack = [source]
        while len(stack) > 0 :
            v = stack.pop()
            visited.append(v)
            if target and v == target:
                return True
            neighbors = self.node_info[v]['next'].keys()
            for w in neighbors:
                if w not in visited:
                    stack.append(w)
        if target:
            return False
        return visited

    # breadth-first search
    def BFS (self, source, target=None):
        visited = []
        queue = [source]
        while len(queue) > 0 :
            v = queue.pop(0)
            visited.append(v)
            if target and v == target:
                return True
            neighbors = self.node_info[v]['next'].keys()
            for w in neighbors:
                if w not in visited:
                    queue.append(w)
        if target:
            return False
        return visited

    # find the shortest path
    def Dijkstra (self, source, target=None):
        dist, prev = {source:0}, {}
        Q = []
        V = self.node_info.keys()
        for v in V:
            if v != source:
                dist[v] = float("inf")
            prev[v] = None
            heapq.heappush(Q, (dist[v], v))
        while len(Q) > 0:
            _, u = heapq.heappop(Q)
            if target and u == target:
                return dist, prev
            for v in self.node_info[u]['next']:
                alt = dist[u] + self.node_info[u]['next'][v]
                if alt < dist[v]:
                    Q.remove((dist[v],v))
                    heapq.heapify(Q)
                    heapq.heappush(Q, (alt, v))
                    dist[v] = alt
                    prev[v] = u
        return dist, prev


class DirectedGraph (Graph):
    # check the graph has cycle (directed graph)
    def has_cycle (self):
        # depth-first search to detect cycle
        def DFS (v):
            V[v] = 0
            for w in self.node_info[v]['next']:
                if V[w] == -1:
                    DFS(w)
                elif V[w] == 0:
                    return True
            V[v] = 1
            return False
        # node state (-1: not, 0: being, 1: fully explored)
        V = {v:-1 for v in self.node_info}
        # start explore nodes
        for v in V:
            if V[v] == -1:
                if DFS(v):
                    return True
        return False
    
    # find all strongly connected components
    def Tarjan (self):
        index = 0
        stack = []
        # node state
        V = {v:{'index':None,
                'lowlink':None,
                'onStack':None} for v in self.node_info}
        SCC_list = []

        # pick unidexed node
        for v in V:
            if V[v]['index'] == None:
                strongconnect(v)

        # recursive depth-first search starting from v
        def strongconnect(v):
            V[v]['index'] = index
            V[v]['lowlink'] = index
            index +=1
            stack.append(v)
            V[v]['onStack'] = True

            # check all neighbors of v
            for w in self.node_info[v]['next']:
                if V[w]['index'] == None:
                    strongconnect(w)
                    V[v]['lowlink'] = min(V[v]['lowlink'], V[w]['lowlink'])
                elif V[w]['onStack'] == True:
                    V[v]['lowlink'] = min(V[v]['lowlink'], V[w]['index'])

            # if v is a root node, pop all nodes in the SCC
            if V[v]['lowlink'] == V[v]['index']:
                SCC = set([])
                while len(stack) > 0:
                    w = stack.pop()
                    V[w]['onStack'] = False
                    assert V[w]['lowlink'] == V[v]['lowlink']
                    SCC.add(w)
                    if w == v:
                        break
                SCC_list.append(SCC)
                    
        return SCC_list
    
class UndirectedGraph (Graph):
    # override
    def __init__ (self,
                  nodes, # node set
                  edges, # edge set
                  node_value={}
                  edge_weight={}): 

        # combine all information into single data structure
        self.node_info = {}
        
        for node in nodes:
            try:
                value = node_value[node]
            except:
                value = None
            self.node_info[node] = {}
            self.node_info[node]['value'] = value # node value
            self.node_info[node]['next'] = {} # next node
            
        for edge in edges:
            node1, node2 = edge
            try:
                weight = edge_weight[node1][node2]
            except:
                try:
                    weight = edge_weight[node2][node1]
                except:
                    weight = None
            self.node_info[node1]['next'][node2] = weight # edge weight
            self.node_info[node2]['next'][node1] = weight # bi-directionality

    # override
    def add_edge (self, node1, node2, weight=None):
        super(UndirectedGraph, self).add_edge(node1, node2, weight=weight)
        super(UndirectedGraph, self).add_edge(node2, node1, weight=weight)

    # override
    def remove_edge (self, node1, node2):
        super(UndirectedGraph, self).remove_edge(node1, node2)
        super(UndirectedGraph, self).remove_edge(node2, node1)

    # check the graph has cycle (undirected graph)
    def has_cycle (self):
        # depth-first search to detect cycle
        def DFS (v):
            V[v] = 0
            for w in self.node_info[v]['next']:
                if V[w] == -1:
                    DFS(w)
                elif w!= v and V[w] == 0:
                    return True
            V[v] = 1
            return False
        # node state (-1: not, 0: being, 1: fully explored)
        V = {v:-1 for v in self.node_info}
        # start to explore nodes
        for v in V:
            if V[v] == -1:
                if DFS(v):
                    return True
        return False

    # find the connected components of graph
    def clustering (self):
        cluster_list = []
        nodes = set(self.node_info.keys())
        while len(nodes) > 0:
            node = nodes.pop()
            cluster = set(self.BFS(node))
            cluster_list.append(cluster)
            nodes = nodes - cluster
        return cluster_list


# load parameter file
def read_para_file(fname, target_names=None, num_max=sys.maxint):
    ID_chr, ID_pos = {}, {}
    name_ID_value = {}
    First = True
    count = 0
    for line in open(fname):
        if count > num_max:
            break
        cols = line.strip().split()
        if First:
            names = cols[3:]
            First = False
            continue
        ID, chr, pos = int(cols[0]), cols[1], int(cols[2])
        ID_chr[ID] = chr
        ID_pos[ID] = pos
        cols = cols[3:]
        for i in range(len(cols)):
            name = names[i]
            if target_names and name not in target_names:
                continue
            if name not in name_ID_value:
                name_ID_value[name] = {}
            assert ID not in name_ID_value[name]
            try:
                value = float(cols[i])
            except:
                value = cols[i]
            name_ID_value[name][ID] = value
        count += 1
    return ID_chr, ID_pos, name_ID_value

def Metropolis_move (move_num, sigma, state, volume, gamma, record=None):
    if record != None:
        note, counter, pre_cycle, edge_prob, ctnum_list, size_cnum, node_size_cnum = record
        f = open('temp_' + str(note) + '_rgs_trace.txt', 'w')
        print >> f, "\n@" + str(counter)
        print >> f, "edges:" + ",".join([str(edge[0]) + '-' + str(edge[1]) for edge in state.edges])

    i = 0
    while i < move_num:
        num = min(int(round(half_right_normal(0, sigma))) + 1, len(state.nodes)*(len(state.nodes)-1)/2)
        #print num
        
        target_edges = {'-':set([]), '+':set([])}
        target_nodes = set([])

        while len(target_edges['+'] | target_edges['-']) < num:
            edge = tuple(sorted(random.sample(state.nodes, 2)))
            if edge in state.edges:
                target_edges['-'].add(edge)
                target_nodes |= set(edge)
            else:
                target_edges['+'].add(edge)
                target_nodes |= set(edge)

        cluster_num = len(state.cluster_nodes)
        node_deg = {node:len(state.node_neighbors[node]) for node in target_nodes}
        potential = state.potential

        # trial change
        #print target_edges
        state.remove_edges(target_edges['-'])
        temp = set([])
        error = False
        for edge in target_edges['+']:
            try:
                state.add_edge(*edge)
                temp.add(edge)
            except:
                state.remove_edges(temp)
                state.add_edges(target_edges['-'])
                error = True
                break
        if error:
            continue

        new_cluster_num = len(state.cluster_nodes)
        new_node_deg = {node:len(state.node_neighbors[node]) for node in target_nodes}
        new_potential = state.potential

        logAcc = 0.0
        for node in target_nodes:
            if new_node_deg[node] >= node_deg[node]:
                start = node_deg[node]
                end = new_node_deg[node]
                sign = 1.0
            else:
                start = new_node_deg[node]
                end = node_deg[node]
                sign = -1.0
            for k in range(start, end):
                logAcc += sign*math.log(State.node_deglimit[node]-k)
        logAcc += (new_cluster_num - cluster_num)*math.log(volume)
        logAcc += -(new_potential - potential)
        logAcc = gamma*logAcc
        
        # go back to original state
        state.remove_edges(target_edges['+'])
        state.add_edges(target_edges['-'])
        
        if record != None:
            counter +=1
            #print >> sys.stderr, "cycle num", counter

        #if record and int(math.log10(counter)) > order:
        #    order = int(math.log10(counter))
        #    print >> sys.stderr, "cycle num", counter

        if random.random() < min(1, math.exp(logAcc)):
            if record != None:
                print >> f, "\n@" + str(counter)
                if len(target_edges['-']) > 0:
                    print >> f, "-:" + ",".join([str(edge[0]) + '-' + str(edge[1]) for edge in target_edges['-']])
                if len(target_edges['+']) > 0:
                    print >> f, "+:" + ",".join([str(edge[0]) + '-' + str(edge[1]) for edge in target_edges['+']])
                
                cur_cycle = counter
                repeat = cur_cycle - max(pre_cycle, burn_in +1)
                if repeat > 0:
                    for edge in state.edges:
                        if edge not in edge_prob:
                            edge_prob[edge] = 0.0
                        edge_prob[edge] += repeat
                    ctnum_list += [len(state.edges)]*repeat
                    for cluster in state.cluster_nodes.values():
                        size = len(cluster)
                        size_cnum[size-1] += repeat
                    for ID, cID in state.node_cluster.items():
                        size = len(state.cluster_nodes[cID])
                        node_size_cnum[ID][size-1] += repeat

                pre_cycle = cur_cycle

                # permanent change
                state.remove_edges(target_edges['-'])
                state.add_edges(target_edges['+'])

        i +=1

    if record != None:
        f.close()
        return state, (counter, pre_cycle, edge_prob, ctnum_list, size_cnum, node_size_cnum)
    return state

def Metropolis_wrapper (arg):
    args, kwargs = arg
    return Metropolis_move(*args, **kwargs)

def concatenate_files (filenames):
    assert len(filenames) > 1
    with open(filenames[0], 'a') as outfile:
        for fname in filenames[1:]:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            os.remove(fname)

def random_graph_simulation (cycle_num,
                             burn_in,
                             sigma,
                             replica_num,
                             root_max,
                             swap_num,
                             move_num,
                             extra_move_num,
                             volume,
                             indis,
                             out_fname):

    # random graph simulation
    print >> sys.stderr, "Start"

    # start writing the trace file    
    f = open(out_fname + '_rgs_trace.txt', 'w')
    print >> f, "IDs" + ":" + ",".join(str(ID) for ID in State.nodes)
    print >> f, "CycleNum:" + str(cycle_num)
    print >> f, "SwapNum:" + str(swap_num)
    print >> f, "MoveNum:" + str(move_num)
    f.close()

    # cycle counter
    counter = 0
    pre_cycle = 1
    order = -1

    # contact probability
    edge_prob = {}
    
    # total contact number
    ctnum_list = []

    # cluster size distribution
    size_cnum = [0]*len(State.nodes) # cluster size distribution
    node_size_cnum = {ID:[0]*len(State.nodes) for ID in State.nodes} # cluster size distribution by particle

    gammas = [1.0]
    for i in range(1, replica_num):
        gamma = (1.0 / float(root_max))**(float(i) / (replica_num-1))
        gammas.append(gamma)
    
    replicas = []
    for i in range(replica_num):
        replicas.append(State(edges = set([]), edge_weight = {}))

    if swap_num > 0:
        assert replica_num > 1
        import multiprocessing as mp
        pool = mp.Pool()

    for i in range(swap_num):
        arg_list = [((move_num, sigma, replicas[0], volume, gammas[0]),
                     {"record":(out_fname + "_" + str(i), counter, pre_cycle, edge_prob, ctnum_list, size_cnum, node_size_cnum)})]
        arg_list += [((move_num, sigma, replicas[k], volume, gammas[k]), {"record":None}) for k in range(1, replica_num)]
        result = pool.map(Metropolis_wrapper, arg_list)
        counter, pre_cycle, edge_prob, ctnum_list, size_cnum, node_size_cnum = result[0][1]
        replicas = [result[0][0]] + result[1:]
        
        concatenate_files([out_fname + '_rgs_trace.txt', 'temp_' + out_fname + "_" + str(i) + '_rgs_trace.txt'])

        for j in range(replica_num-1):
            logAcc = 0.0
            for node in State.nodes:
                if len(replicas[j+1].node_neighbors[node]) >= len(replicas[j].node_neighbors[node]):
                    start = len(replicas[j].node_neighbors[node])
                    end = len(replicas[j+1].node_neighbors[node])
                    sign = 1.0
                else:
                    start = len(replicas[j+1].node_neighbors[node])
                    end = len(replicas[j].node_neighbors[node])
                    sign = -1.0
                for k in range(start, end):
                    logAcc += sign*math.log(State.node_deglimit[node]-k)
            logAcc += (len(replicas[j+1].cluster_nodes) - len(replicas[j].cluster_nodes))*math.log(volume)
            logAcc += - (replicas[j+1].potential - replicas[j].potential)
            logAcc = (gammas[j] - gammas[j+1])*logAcc

            if j == 0:
                #print counter
                counter +=1
                if int(math.log10(counter)) > order:
                    order = int(math.log10(counter))
                    print >> sys.stderr, "cycle num", counter
            
            if random.random() < min(1, math.exp(logAcc)):
                #print "Swap!"
                if j == 0:
                    cur_cycle = counter
                    repeat = cur_cycle - max(pre_cycle, burn_in+1)
                    if repeat > 0:                    
                        for edge in replicas[0].edges:
                            if edge not in edge_prob:
                                edge_prob[edge] = 0.0
                            edge_prob[edge] += repeat
                        ctnum_list += [len(replicas[0].edges)]*repeat
                        for cluster in replicas[0].cluster_nodes.values():
                            size = len(cluster)
                            size_cnum[size-1] += repeat
                        for ID, cID in replicas[0].node_cluster.items():
                            size = len(replicas[0].cluster_nodes[cID])
                            node_size_cnum[ID][size-1] += repeat

                    pre_cycle = cur_cycle
                    
                replicas[j], replicas[j+1] = replicas[j+1], replicas[j]
    
    result = Metropolis_move(extra_move_num, sigma, replicas[0], volume, gammas[0], record=(out_fname + '_extra', counter, pre_cycle, edge_prob, ctnum_list, size_cnum, node_size_cnum))
    counter, pre_cycle, edge_prob, ctnum_list, size_cnum, node_size_cnum = result[1]

    concatenate_files([out_fname + '_rgs_trace.txt', 'temp_' + out_fname + '_extra_rgs_trace.txt'])
    
    assert counter == cycle_num

    cur_cycle = counter
    repeat = cur_cycle - max(pre_cycle, burn_in+1) + 1
    if repeat > 0:
        for edge in result[0].edges:
            if edge not in edge_prob:
                edge_prob[edge] = 0.0
            edge_prob[edge] += repeat
        ctnum_list += [len(result[0].edges)]*repeat
        for cluster in result[0].cluster_nodes.values():
            size = len(cluster)
            size_cnum[size-1] += repeat
        for ID, cID in result[0].node_cluster.items():
            size = len(result[0].cluster_nodes[cID])
            node_size_cnum[ID][size-1] += repeat

    assert len(ctnum_list) == max(0, cycle_num - burn_in)

    for edge in edge_prob:
        edge_prob[edge] = edge_prob[edge] / float(cycle_num - burn_in)
    for k in range(len(size_cnum)):
        size_cnum[k] = size_cnum[k] / float(cycle_num - burn_in)
    for ID in node_size_cnum:
        for k in range(len(node_size_cnum[ID])):
            node_size_cnum[ID][k] = node_size_cnum[ID][k] / float(cycle_num - burn_in)
    

    # writing the averaged data
    f = open(out_fname + '_rgs_ctprob.txt', 'w')
    print >> f, "IDs" + ":" + ",".join(str(ID) for ID in State.nodes)

    for edge in edge_prob:
        print >> f, str(edge[0]) + '-' + str(edge[1]) + ':' + str(edge_prob[edge])

    f.close()

    f = open(out_fname + '_rgs_ctnum.txt', 'w')
    print >> f, "CycleNum:" + str(cycle_num)
    print >> f, ",".join(str(value) for value in ctnum_list)

    f.close()

    f = open(out_fname + '_rgs_cluster.txt', 'w')
    print >> f, "IDs" + ":" + ",".join(str(ID) for ID in State.nodes)
    print >> f, "@total"
    print >> f, ",".join(str(cnum) for cnum in size_cnum)
    for ID in State.nodes:
        print >> f, "@" + str(ID)
        print >> f, ",".join(str(cnum) for cnum in node_size_cnum[ID])

    f.close()

    print >> sys.stderr, "Done"
    print >> sys.stderr


if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = ArgumentParser(description='Random graph simulation')
    parser.add_argument(metavar='--para',
                        dest="para_fname",
                        type=str,
                        help='parameter file')
    parser.add_argument('--cycle',
                        dest="cycle_num",
                        type=int,
                        default=1000,
                        help='iteration cycle number')
    parser.add_argument('--burn',
                        dest="burn_in",
                        type=int,
                        help='burn-in cycle number')
    parser.add_argument('--sigma',
                        dest="sigma",
                        type=float,
                        default=1,
                        help='Metropolis step size STD')
    parser.add_argument('--replica',
                        dest="replica_num",
                        type=int,
                        default=5,
                        help='number of replicas')
    parser.add_argument('--max-root',
                        dest="root_max",
                        type=float,
                        default=10,
                        help='maximum power of rooting')
    parser.add_argument('--swap',
                        dest="swap_freq",
                        type=float,
                        default=0.1,
                        help='replica exchange frequency')
    parser.add_argument('--ptls',
                        dest="particles",
                        type=str,
                        default="",
                        help='particle selection (, or -)')
    parser.add_argument('--volume',
                        dest="volume",
                        type=float,
                        default=10**7,
                        help='Solution volume in unit of particle volume')
    parser.add_argument('--bethe',
                        dest="bethe",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='Bethe-lattice constraint')
    parser.add_argument('--indis-site',
                        dest="indis",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='indistinguishable binding sites')    
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    if args.burn_in == None:
        burn_in = int(args.cycle_num*0.1) # trash first 10% of total cycle
    else:
        burn_in = args.burn_in

    # if replica num is 1, then swap num is 0
    if args.replica_num < 2:
        swap_num = 0
    else:
        swap_num = int(args.swap_freq*args.cycle_num)

    # if swap num is 0, then replica num is 1
    if swap_num <= 0:
        replica_num = 1
        move_num = 0
        extra_move_num = args.cycle_num
    else:
        replica_num = args.replica_num
        move_num = args.cycle_num / swap_num - 1
        extra_move_num = args.cycle_num % swap_num

    if len(args.particles) > 0:
        particles = []
        for word in args.particles.split(','):
            if '-' in word:
                st, ed = word.split('-')
                st, ed = int(st), int(ed)
                particles += range(st, ed+1)
            particles.append(int(word))
    else:
        particles = None

    ID_chr, ID_pos, name_ID_value = read_para_file(args.para_fname)
    if particles == None:
        particles = ID_chr.keys()


    # Bethe-lattice Cluster object
    class BetheCluster (UndirectedGraph):
        # override
        def __init__ (self,
                      nodes, # node set
                      edges, # edge set
                      node_value={}
                      edge_weight={}):

            # record node and edge information
            super(BetheCluster, self).__init__(nodes,
                                               edges,
                                               node_value=node_value,
                                               edge_weight=edge_weight)
            # get the potential of graph
            self.potential = self.get_potential()

        # to do
        def __add__ (self, other):
            assert type(other) == BetheCluster
            return

        def __len__ (self):
            return len(self.nodes)

        def get_potential (self):
            potential = 0.0
            for node in self.node_info:
                potential += 0.5*self.node_info[node]['value']*len(self.node_info[node]['next'])
            return potential


    # Gaussian Network object
    class GaussianNet (UndirectedGraph):
        # override
        def __init__ (self,
                      nodes, # node set
                      edges, # edge set
                      node_value={}
                      edge_weight={}):

            # save as node_info dictionary
            super(GaussianNet, self).__init__(nodes,
                                              edges,
                                              node_value=node_value,
                                              edge_weight=edge_weight)

            # get the partition funcition
            self.partition = self.get_partition()

        # get the Kirchhoff matrix of Gaussian network
        def get_Kirchhoff_matrix (self):
            nodes = self.node_info.keys()
            N = len(nodes)
            matrix = np.zeors((N, N))
            for i in range(N):
                for j in range(N):
                    node1, node2 = nodes[i], nodes[j]
                    if i == j:
                        matrix[i][j] = sum(self.node_info[node1]['next'].values())
                    else:
                        matrix[i][j] = -self.node_info[node1]['next'][node2]
            return matrix

        # get the partition funcition of Gaussian network
        def get_partition (self):
            matrix = self.get_Kirchhoff_matrix()[1:,1:] # remove first row and column

            # eigendecomposition of the matrix
            eigenvalues, eigen_matrix = LA.eigh(matrix)

            # compute the determinant of inverse matrix
            det = 1.0
            for eigenvalue in eigenvalues:
                assert eigenvalue != 0
                assert np.isreal(eigenvalue)
                det *= 1.0/eigenvalue

            # compute the parition function
            N = len(self.node_info)
            GN_partition = ((2*math.pi)**(1.5*(N-1)))*(det**1.5)        

            return GN_partition

        
    # Worm-like Network object (Todo)
    #class WormLikeNet (object)


    # Cluster State object
    class State (object):
        # make static variables shared with all instances
        # the detail information of particles
        node_stinfo = {}
        for ID in particles:
            node_stinfo[ID] = {}
            node_stinfo[ID]['score'] = name_ID_value['Score'][ID]
            node_stinfo[ID]['valency'] = name_ID_value['Valency'][ID]
            node_stinfo[ID]['pos'] = name_ID_value['Position'][ID]

        def __init__ (self):

            # cluster information
            self.cID_cluster, self.node_cID = {}, {}
            for node in node_stinfo:
                nodes = set([node])
                node_value = {node:copy.copy(node_stinfo[ID]['score'])} # shallow copy
                cluster = BetheCluster(nodes=nodes, node_value=node_value)
                cID = cluster.name
                self.cID_cluster[cID] = cluster
                self.node_cID[node] = cID
                

            # polymer network information
            self.polymer_network = GaussianNet(nodes=set([]),
                                               edges=set([]))

            # sorted list of contacting nodes
            self.sorted_contact_nodes = []
            self.idx_cnode, self.cnode_idx = [], {}

            # total potential of the state
            self.potential = 0.0

        # merge two clusters
        def merge (self, cID1, cID2):
            # update cluster information
            cluster1, cluster2 = self.cID_cluster[cID1], self.cID_cluster[cID2]
            new_cluster = cluster1 + cluster2
            new_cID = cluster.name
            del self.cID_cluster[cID1]
            del self.cID_cluster[cID2]
            self.cID_cluster[new_cID] = new_cluster
            for node in cluster.nodes:
                self.node_cID[node] = new_cID

            # update contact node information
            for cID in [cID1, cID2]:
                if len(self.cID_cluster[cID]) < 2:
                    new_cnode = self.cID_cluster[cID]
                    idx = bisect.bisect_left (self.idx_cnode, new_cnode)
                    self.idx_cnode.insert(idx, new_cnode)
                    for cnode in self.idx_cnode[idx+1:]:
                        self.cnode_idx[cnode] +=1
            
            # update polymer network 
            self.polymer_network.add_node(new_cID)
            for cID in [cID1, cID2]:
                if len(self.cID_cluster[cID]) < 2:
                    cnode = self.cID_cluster[cID]
                    idx = cnode_idx[cnode]
                    nn_cnodes = self.idx_cnode[max(0, idx-1):idx] + self.idx_cnode[idx+1:idx+2]
                    for nn_cnode in nn_cnodes:
                        next_cID = self.node_cID[nn_code]
                        weight = abs(self.node_stinfo[nn_cnode]['pos'] - self.node_stinfo[cnode]['pos'])
                        

                try:
                    for next_cID, weight in self.polymer_network.node_info[cID]['next'].items():
                        if next_cID not in self.polymer_network.node_info[new_cID]['next']:
                            self.polymer_network.add_edge(new_cID, next_cID, weight=0)
                        self.polymer_network.node_info[new_cID]['next'][next_cID] += weight
                except:
                    for node in self.cID_cluster[cID].nodes:
                        idx = binary_search (self.sorted_contact_nodes, node)
                        self.sorted_contanct_nodes.insert(idx, node)
                        nn_nodes = self.sorted_contact_nodes[max(0, idx-1):idx]
                        nn_nodes += self.sorted_contact_nodes[idx+1:idx+2]
                        for nn_node in nn_nodes:
                            next_cID = self.node_cID[nn_node]
                            weight = abs(self.node_stinfo[nn_node]['pos'] - self.node_stinfo[node]['pos'])
                            if next_cID not in self.polymer_network.node_info[new_cID]['next']:
                                self.polymer_network.add_edge(new_cID, next_cID, weight=0)
                            self.polymer_network.node_info[new_cID]['next'][next_cID] += weight

            self.polymer_network.remove_node(cID1)
            self.polymer_network.remove_node(cID2)
            
            return

        # split the cluster into two
        def split (self, cID):
            # update cluster information
            cluster1, cluster2 = self.cID_cluster[cID]/2
            cID1, cID2 = cluster1.name, cluster2.name
            del self.cID_cluster[cID]
            self.cID_cluster[cID1] = cluster1
            self.cID_cluster[cID2] = cluster2

            # update polymer network
            self.polymer_network.add_node(cID1)
            for node in self.cID_cluster[cID1].nodes:
                idx = binary_search (self.sorted_contact_nodes, node)
                self.sorted_contanct_nodes.insert(idx, node)
                nn_nodes = self.sorted_contact_nodes[max(0, idx-1):idx]
                nn_nodes += self.sorted_contact_nodes[idx+1:idx+2]
                for nn_node in nn_nodes:
                    next_cID = self.node_cID[nn_node]
                    weight = abs(self.node_stinfo[nn_node]['pos'] - self.node_stinfo[node]['pos'])
                    if next_cID not in self.polymer_network.node_info[new_cID]['next']:
                        self.polymer_network.add_edge(new_cID, next_cID, weight=0)
                    self.polymer_network.node_info[new_cID]['next'][next_cID] += weight

            
            
            return 

        # get polymer network
        def get_polymer_network (self):
            return

            
            

    # Contact State object
    class State (UndirectedGraph):
        # make static variables shared with all instances
        node_stinfo = {}
        for ID in particles:
            node_stinfo[ID] = {}
            node_stinfo[ID]['score'] = name_ID_value['Score'][ID]
            node_stinfo[ID]['valency'] = name_ID_value['Valency'][ID]
            node_stinfo[ID]['pos'] = name_ID_value['Position'][ID]
            
        no_cycle = args.bethe
        
        def __init__ (self, edges):

            # record the information of nodes and edges
            self.node_info = {}
            for node in node_stinfo:
                self.node_info[node] = {}
                # shallow copy of static variables (score)
                self.node_info[node]['value'] = copy.copy(node_stinfo[node]['score'])
                # shallow copy of static variables (valency)
                self.node_info[node]['deglimit'] = copy.copy(node_stinfo[node]['valency'])
                self.node_info[node]['next'] = {}

            for edge in edges:
                node1, node2 = edge
                self.node_info[node1]['next'][node2] = None
                self.node_info[node2]['next'][node1] = None

            # check the degree limit violation
            for node in self.node_info:
                assert len(self.node_info[node]['next']) <= self.node_info[node]['deglimit']

            # check the no cycle violation
            if self.no_cycle:
                assert not self.has_cycle()

            # find the connected components
            self.cID_nodes = self.clustering()

            # sorted contact node list
            self.cnodes = []

            # build the polymer network
            self.polymer_net = self.get_polymer_network()

            # compute the potential of the state
            self.potential = self.get_potential()

        def get_polymer_network (self):
            PN_nodes, PN_edges = set([]), set([])
            PN_edge_weight = {}
            cnode_pt = None
            for node in sorted(self.node_info.keys()):
                cluster_size = self.node_info[node]['cID'][0]
                if cluster_size < 2:
                    continue
                if cnode_pt == None:
                    cnode_pt = node
                    continue
                cID1, cID2 = self.node_info[cnode_pt]['cID'], self.node_info[node]['cID']
                if cID1 == cID2:
                    cnode_pt = node
                    continue
                PN_nodes.add(cID1)
                PN_nodes.add(cID2)
                edge = (cID1, cID2)
                PN_edges.add(edge)
                weight = self.node_sinfo[node] - self.node_sinfo[cnode_pt]
                if edge not in PN_edge_weight:
                    PN_edge_weight[edge] = 0
                PN_edge_weight[edge] += weight
                cnode_pt = node
            return GaussianNet(PN_nodes, PN_edges, edge_weight=PN_edge_weight)

        def find_nn_clusters (self, cID, min_size=2):
            def linear_search (start_node, step):
                find, weight = None, None
                pt = start_node + step
                while pt >= 0 and pt < len(self.node_info):
                    next_cID = self.node_info[pt]['cID']
                    if len(self.cID_nodes['cID']) >= min_size:
                        find = next_cID
                        weight = self.node_sinfo[pt]['pos'] - self.node_sinfo[start_node]['pos']
                        weight = abs(weight)
                        break
                    pt += step
                return find, weight

            next_weight = {}
            start_nodes = sorted(self.cID_nodes[cID])
            pt = 0
            while pt < len(start_nodes)*2:
                start_node = start_nodes[pt/2]
                if pt % 2 == 0:
                    step = -1 # left search
                else:
                    step = 1 # right search
                next_cID, weight = linear_search (start_node, step)
                if next_cID == None:
                    pt +=1
                    continue
                if next_cID == cID:
                    assert pt != 0 and len(start_nodes)*2 - 1 
                    assert step == 1
                    pt +=2
                    continue
                if next_cID not in next_weight:
                    next_weight[next_cID] = 0
                next_weight[next_cID] += weight
                pt +=1

            return next_weight
                    
                
            


        
        def add_edge (self, node1, node2, weight=None):
            # check the degree limit violation
            assert len(self.node_info[node1]['next']) < self.node_info[node1]['deglimit']
            assert len(self.node_info[node2]['next']) < self.node_info[node2]['deglimit']

            # check the no-cycle violation
            if self.no_cycle:
                #assert self.node_cluster[node1] != self.node_cluster[node2]
                assert self.node_info[node1]['cID'] != self.node_info[node2]['cID']

            # add edge
            super(State, self).add_edge(node1, node2, weight=weight)

            # check clustering change
            merge = False
            cID1, cID2 = self.node_info[node1]['cID'], self.node_info[node2]['cID']
            if self.no_cycle or cID1 != cID2:
                merge = True

            if merge:

                # update clustering
                cluster1, cluster2 = self.cID_nodes[cID1], self.cID_nodes[cID2] 
                new_cluster = cluster1 | cluster2
                new_cID = (len(new_cluster), min(new_cluster))
                del self.cID_nodes[cID1]; del self.cID_nodes[cID2]
                self.cID_nodes[new_cID] = new_cluster
                for node in new_cluster:
                    self.node_info[node]['cID'] = new_cID

                # update polymer-network
                if len(cluster1) > 1:
                    next_weight1 = self.polymer_net.node_info[cID1]['next']
                else:
                    next_weight1 = self.find_nn_cnodes(cluster1)
                    if len(next_weight1) > 1:
                        assert len(next_weight1) == 2
                        n1, n2 = next_weight1.keys()
                        self.polymer_net.remove_edge(n1, n2)
                    
                if len(cluster2) > 1:
                    next_weight2 = self.polymer_net.node_info[cID2]['next']
                else:
                    next_weight2 = self.find_nn_cnodes(cluster2)
                    if len(next_weight2) > 1:
                        assert len(next_weight2) == 2
                        n1, n2 = next_weight2.keys()
                        self.polymer_net.remove_edge(n1, n2)
                
                next_weight = add_dicts(next_weight1, next_weight2)
                
                self.polymer_net.add_node(new_cID)
                for next, weight in next_weight.items():
                    self.polymer_net.add_edge(new_cID, next, weight=weight)
                
                self.polymer_net.remove_node(cID1)
                self.polymer_net.remove_node(cID2)
            

        def remove_edge (self, node1, node2):
            super(State, self).remove_edge(node1, node2)

            # check clustering change
            split = False
            cluster1 = set(self.BFS(node1))
            if self.no_cycle or node2 not in cluster1:
                split = True

            if split:

                # update clustering
                cID = self.node_info[node1]['cID']
                cluster = self.cID_nodes[cID]
                cluster2 = cluster - cluster1
                cID1 = (len(cluster1), min(cluster1))
                cID2 = (len(cluster2), min(cluster2))
                del self.cID_nodes[cID]
                self.cID_nodes[cID1] = cluster1
                self.cID_nodes[cID2] = cluster2
                for node in cluster1:
                    self.node_info[node]['cID'] = cID1
                for node in cluster2:
                    self.node_info[node]['cID'] = cID2

                # update polymer-network                
                next_weight = self.polymer_net.node_info[cID]['next']
                    
                if len(cluster1) <=  len(cluster2):
                    next_weight1 = find_nns (cluster1)
                    next_weight2 = subtract_dicts (next_weight, next_weight1)
                else:
                    next_weight2 = find_nns (cluster2)
                    next_weight1 = subtract_dicts (next_weight, next_weight2)

                if len(cluster1) > 1:
                    self.polymer_net.add_node(cID1)
                    for next, weight in next_weight1.items():
                        self.polymer_net.add_edge(cID1, next, weight=weight)
                else:
                    if len(next_weight1) > 1:
                        assert len(next_weight1) == 2
                        n1, n2 = next_weight1.keys()
                        weight = sum(next_weight1.values())
                        self.polymer_net.add_edge(n1, n2, weight=weight)

                if len(cluster2) > 1:
                    self.polymer_net_add_node(cID2)
                    for next, weight in next_weight2.items():
                        self.polymer_net.add_edge(cID2, next, weight=weight)
                else:
                    if len(next_weight2) > 1:
                        assert len(next_weight2) == 2
                        n1, n2 = next_weight2.keys()
                        weight = sum(next_weight2.values())
                        self.polymer_net.add_edge(n1, n2, weight=weight)
                
                self.polymer_net.remove_node(cID)
                    

            return 
        
        def add_node (self, node):
            raise AttributeError("'State' object has no attribute 'add_node'")
        def remove_node (self, node):
            raise AttributeError("'State' object has no attribute 'remove_node'")

        def get_potential (self):
            potential = 0.0
            for edge in self.edges:
                node1, node2 = edge
                try:
                    value1, value2 = self.node_value[node1], self.node_value[node2]
                except:
                    pass
                potential += 0.5*(value1 + value2)
            return potential


    print >> sys.stderr
    print >> sys.stderr, "Random graph simulation"
    print >> sys.stderr
    print >> sys.stderr, "input file: " + str(args.para_fname)
    print >> sys.stderr, "cycle number: " + str(args.cycle_num)
    if swap_num > 0:
        print >> sys.stderr, "replica number: " + str(replica_num)
        print >> sys.stderr, "%s repeat: %s move and 1 swap" % (str(swap_num), str(move_num))
    if len(args.particles) > 0:
        print >> sys.stderr, "particles: " + str(particles)
    else:
        print >> sys.stderr, "particle number: " + str(len(particles))
    print >> sys.stderr, "simulation volume: 10^%s" % str(math.log10(args.volume))
    print >> sys.stderr, "Bethe lattice: " + str(args.bethe)
    print >> sys.stderr, "identical binding sites: " + str(args.indis)
    print >> sys.stderr

    pr = cProfile.Profile()
    pr.enable()  # start profiling
    
    random_graph_simulation (args.cycle_num,
                             burn_in,
                             args.sigma,
                             replica_num,
                             args.root_max,
                             swap_num,
                             move_num,
                             extra_move_num,
                             args.volume,
                             args.indis,
                             args.out_fname)

    
    pr.disable()  # end profiling
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print s.getvalue()

