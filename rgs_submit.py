import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import random
import copy
import heapq
import pickle
import cProfile, pstats, StringIO
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

def half_right_normal (mu, sigma):
    value = random.gauss(mu, sigma)
    if value < mu:
        return 2*mu - value
    return value
        
class Graph (object):
    def __init__ (self,
                  nodes,
                  edges,
                  node_value = {},
                  edge_weight = {},
                  directed = False):
        
        self.nodes = nodes # node set
        self.edges = edges # edge set
        self.node_value = node_value # node value
        self.edge_weight = edge_weight # edge weight
        self.directed = directed # directed graph option

        # find neighbors for each node
        self.node_neighbors = self.get_neighbors()

        # find the connected components
        self.cluster_nodes, self.node_cluster = self.clustering()

        # get the potential of graph
        self.potential = self.get_potential()

    def add_node (self, node):
        self.nodes.add(node)
        return self.nodes

    def add_nodes (self, nodes):
        self.nodes |= set(nodes)
        return self.nodes

    def add_edge (self, node1, node2, weight=None):
        # two nodes should be in node set
        assert set([node1, node2]) <= self.nodes
        # if already edge between them, do nothing
        if self.directed:
            if (node1, node2) in self.edges:
                return self.edges
        else:
            if tuple(sorted((node1, node2))) in self.edges:
                return self.edges
        # update edge set
        if self.directed:
            edge = (node1, node2)
        else:
            edge = tuple(sorted((node1, node2)))
        self.edges.add(edge)
        # update edge weight
        if weight != None:
            if node1 not in self.edge_weight:
                self.edge_weight[node1] = {}
            self.edge_weight[node1][node2] = weight
            if not self.directed:
                if node2 not in self.edge_weight:
                    self.edge_weight[node2] = {}
                self.edge_weight[node2][node1] = weight
        # update neighbors
        self.node_neighbors[node1].add(node2)
        self.node_neighbors[node2].add(node1)
        # update clusters
        cID1, cID2 = self.node_cluster[node1], self.node_cluster[node2]
        if cID1 != cID2:
            new_cluster = self.cluster_nodes[cID1] | self.cluster_nodes[cID2]
            new_cID = str(len(new_cluster)) + '-' + str(min(new_cluster))
            del self.cluster_nodes[cID1]; del self.cluster_nodes[cID2]
            self.cluster_nodes[new_cID] = new_cluster
            for node in new_cluster:
                self.node_cluster[node] = new_cID
        # update the potential
        try:
            value1, value2 = self.node_value[node1], self.node_value[node2]
        except:
            pass
        self.potential += 0.5*(value1 + value2)
        return self.edges

    def add_edges (self, edges, weight=None):
        for edge in edges:
            node1, node2 = edge
            self.add_edge(node1, node2, weight=weight)
        return self.edges

    def add_edges_with_weight (self, edge_weight):
        for node1 in edge_weight:
            for node2 in edge_weight[node1]:
                weight = edge_weight[node1][node2]
                self.add_edge(node1, node2, weight=weight)
        return self.edge_weight

    #To do
    def remove_node (self, node):
        self.nodes -= set([node])
        for node1 in self.edge_weight:
            if node1 == node:
                del self.edge_weight[node1]
                continue
            for node2 in self.edge_weight[node1]:
                if node2 == node:
                    del self.edge_weight[node1][node2]
        return self.nodes

    #To do
    def remove_nodes (self, nodes):
        self.nodes -=  set(nodes)
        for node1 in self.edge_weight:
            if node1 in nodes:
                del self.edge_weight[node1]
                continue
            for node2 in self.edge_weight[node1]:
                if node2 in nodes:
                    del self.edge_weight[node1][node2]
        return self.nodes

    def remove_edge (self, node1, node2):
        # two nodes should be in node set
        assert set([node1, node2]) <= self.nodes
        # if no edge between them, do nothing
        if self.directed:
            if (node1, node2) not in self.edges:
                return self.edges
        else:
            if tuple(sorted((node1, node2))) not in self.edges:
                return self.edges
        # update edge set
        if self.directed:
            self.edges.remove((node1, node2))
        else:
            self.edges.remove(tuple(sorted([node1, node2])))
        # update edge_weight
        try:
            del self.edge_weight[node1][node2]
            if len(self.edge_weight[node1]) <= 0:
                del self.edge_weight[node1]
        except:
            pass
        if not self.directed:
            try:
                del self.edge_weight[node2][node1]
                if len(self.edge_weight[node2]) <= 0:
                    del self.edge_weight[node2]
            except:
                pass
        # update neighbors
        self.node_neighbors[node1].remove(node2)
        self.node_neighbors[node2].remove(node1)
        # update clusters
        cluster1 = self.BFS(node1)
        if node2 not in cluster1:
            cluster2 = self.BFS(node2)
            new_cID1 = str(len(cluster1)) + '-' + str(min(cluster1))
            new_cID2 = str(len(cluster2)) + '-' + str(min(cluster2))
            cID1, cID2 = self.node_cluster[node1], self.node_cluster[node2]
            assert cID1 == cID2
            del self.cluster_nodes[cID1]
            self.cluster_nodes[new_cID1] = cluster1
            self.cluster_nodes[new_cID2] = cluster2
            for node in cluster1:
                self.node_cluster[node] = new_cID1
            for node in cluster2:
                self.node_cluster[node] = new_cID2
        # update the potential
        try:
            value1, value2 = self.node_value[node1], self.node_value[node2]
        except:
            pass
        self.potential -= 0.5*(value1 + value2)
        return self.edges

    def remove_edges (self, edges):
        for edge in edges:
            node1, node2 = edge
            self.remove_edge(node1, node2)
        return self.edges

    def DFS (self, source, target=None):
        visited = set([])
        stack = [source]
        while len(stack) > 0 :
            v = stack.pop()
            visited.add(v)
            if target and v == target:
                return True
            neighbors = self.node_neighbors[v]
            for w in neighbors:
                if w not in visited:
                    stack.append(w)
        if target:
            return False
        return visited

    def BFS (self, source, target=None):
        visited = set([])
        queue = [source]
        while len(queue) > 0 :
            v = queue.pop(0)
            visited.add(v)
            if target and v == target:
                return True
            neighbors = self.node_neighbors[v]
            for w in neighbors:
                if w not in visited:
                    queue.append(w)
        if target:
            return False
        return visited

    def has_cycle (self):
        nodes = copy.deepcopy(self.nodes)
        while len(nodes) > 0:
            source = nodes.pop()
            visited = set([])
            stack = [source]
            previous = None
            while len(stack) > 0 :
                v = stack.pop()
                visited.add(v)
                neighbors = self.node_neighbors[v]
                for w in neighbors:
                    if w != previous and w in visited:
                        return True
                    if w not in visited:
                        stack.append(w)
                previous = copy.deepcopy(v)
            nodes -= visited
        return False

    def Dijkstra (self, source, target=None):
        dist, prev = {source:0}, {}
        Q = []
        V = self.edge_weight.keys()
        for v in V:
            if v != source:
                dist[v] = float("inf")
            prev[v] = None
            heapq.heappush(Q, (dist[v], v))
        while len(Q) > 0:
            _, u = heapq.heappop(Q)
            if target and u == target:
                return dist, prev
            for v in self.edge_weight[u]:
                alt = dist[u] + min(self.edge_weight[u][v])
                #alt = dist[u] + min(v1_v2_weights[u][v].values())
                if alt < dist[v]:
                    Q.remove((dist[v],v))
                    heapq.heapify(Q)
                    heapq.heappush(Q, (alt, v))
                    dist[v] = alt
                    prev[v] = u
        return dist, prev

    def get_neighbors (self):
        node_neighbors = {node:set([]) for node in self.nodes} 
        for edge in self.edges:
            node1, node2 = edge
            node_neighbors[node1] |= set([node2])
            node_neighbors[node2] |= set([node1])
        return node_neighbors

    def clustering (self):
        cluster_nodes = {} # cluster to nodes belong to it
        node_cluster = {} # node to cluster it belong to
        nodes = copy.deepcopy(self.nodes)
        while len(nodes) > 0:
            node = nodes.pop()
            cluster = self.BFS(node)
            size = len(cluster)
            cID = str(size) + '-' + str(min(cluster))
            assert cID not in cluster_nodes
            cluster_nodes[cID] = cluster
            for node in list(cluster):
                node_cluster[node] = cID
            nodes = nodes - cluster
        return cluster_nodes, node_cluster

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

    class State(Graph):
        # make static variables shared with all instances
        nodes = set(particles)
        node_value = {ID:name_ID_value["Score"][ID] for ID in particles}
        directed = False
        node_deglimit = {ID:name_ID_value["Valency"][ID] for ID in particles}
        no_cycle = args.bethe
        #combining_rule = combining_rule
        
        def __init__ (self, edges, edge_weight = {}):
            self.edges = edges
            self.edge_weight = edge_weight
            self.node_neighbors = self.get_neighbors()
            self.saturated_nodes = set([])
            for node in self.node_deglimit:
                assert len(self.node_neighbors[node]) <= self.node_deglimit[node]
                if len(self.node_neighbors[node]) == self.node_deglimit[node]:
                    self.saturated_nodes.add(node)
            if self.no_cycle:
                assert not self.has_cycle()
            self.cluster_nodes, self.node_cluster = self.clustering()
            self.potential = self.get_potential()

        def add_edge (self, node1, node2, weight=None):
            # check the degree limit violation
            if self.node_deglimit.get(node1) != None:
                assert len(self.node_neighbors[node1]) < self.node_deglimit[node1]
            if self.node_deglimit.get(node2) != None:
                assert len(self.node_neighbors[node2]) < self.node_deglimit[node2]
            # check the no cycle violation
            if self.no_cycle:
                assert self.node_cluster[node1] != self.node_cluster[node2]
            super(State, self).add_edge(node1, node2, weight=weight)
            # update saturated nodes
            if self.node_deglimit[node1] == len(self.node_neighbors[node1]):
                self.saturated_nodes.add(node1)
            if self.node_deglimit[node2] == len(self.node_neighbors[node2]):
                self.saturated_nodes.add(node2)
            return self.edges

        def remove_edge (self, node1, node2):
            super(State, self).remove_edge(node1, node2)
            # update saturated nodes
            if self.node_deglimit[node1] == len(self.node_neighbors[node1]) + 1:
                self.saturated_nodes.remove(node1)
            if self.node_deglimit[node2] == len(self.node_neighbors[node2]) + 1:
                self.saturated_nodes.remove(node2)
            return self.edges
        
        def add_node (self, node):
            raise AttributeError("'State' object has no attribute 'add_node'")
        def add_nodes (self, nodes):
            raise AttributeError("'State' object has no attribute 'add_nodes'")
        def remove_node (self, node):
            raise AttributeError("'State' object has no attribute 'remove_node'")
        def remove_nodes (self, nodes):
            raise AttributeError("'State' object has no attribute 'remove_nodes'")

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

    #pr = cProfile.Profile()
    #pr.enable()  # start profiling
    
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

    
    #pr.disable()  # end profiling
    #s = StringIO.StringIO()
    #sortby = 'cumulative'
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    #ps.print_stats()
    #print s.getvalue()

