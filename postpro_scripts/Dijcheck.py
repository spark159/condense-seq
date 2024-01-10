import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import math
import random
import copy
import heapq
import numpy as np

def Dijkstra(v1_v2_weight, source, target=None):
    dist, prev = {source:0}, {}
    Q = []
    V = v1_v2_weight.keys()
    for v in V:
        if v != source:
            dist[v] = np.inf
        prev[v] = None
        heapq.heappush(Q, (dist[v], v))
    while len(Q) > 0:
        _, u = heapq.heappop(Q)
        if target and u == target:
            return dist, prev
        for v in v1_v2_weight[u]:
            alt = dist[u] + v1_v2_weight[u][v]
            if alt < dist[v]:
                Q.remove((dist[v],v))
                heapq.heapify(Q)
                heapq.heappush(Q, (alt, v))
                dist[v] = alt
                prev[v] = u
    return dist, prev

v1_v2_weight = {}

v1_v2_weight[1] = {}
v1_v2_weight[1][2] = 7
v1_v2_weight[1][3] = 9
v1_v2_weight[1][6] = 14

v1_v2_weight[2] = {}
v1_v2_weight[2][1] = 7
v1_v2_weight[2][3] = 10
v1_v2_weight[2][4] = 15

v1_v2_weight[3] = {}
v1_v2_weight[3][1] = 9
v1_v2_weight[3][2] = 10
v1_v2_weight[3][4] = 11
v1_v2_weight[3][6] = 2

v1_v2_weight[4] = {}
v1_v2_weight[4][2] = 15
v1_v2_weight[4][3] = 11
v1_v2_weight[4][5] = 6

v1_v2_weight[5] = {}
v1_v2_weight[5][4] = 6
v1_v2_weight[5][6] = 9

v1_v2_weight[6] = {}
v1_v2_weight[6][1] = 14
v1_v2_weight[6][3] = 2
v1_v2_weight[6][5] = 9


dist, prev = Dijkstra(v1_v2_weight, source=1, target=None)
