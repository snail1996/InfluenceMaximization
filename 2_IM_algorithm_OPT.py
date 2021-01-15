import numpy as np
from tqdm import tqdm
from igraph import *
 
def triggering(G, S, p, mc):

    spread = []
    
    for _ in tqdm(range(mc)):

        new_nodes, RRS0 = S[:], S[:]
        
        while new_nodes:                            #求所选节点的激活节点集合
            
            B = G.neighborhood(vertices=new_nodes, order=1, mode="out", mindist=0)   #新节点的邻居节点集合

            A = [item for subset in B for item in subset] #新节点的邻居节点集合，转化为一个列表

            success = np.random.uniform(0, 1, len(A)) < p

            temp = list(np.extract(success, A))           #新节点以概率p随机激活邻居节点，每个节点只能激活一次

            RRS = list(set(RRS0 + temp))      #list(set())可以去除相同元素

            new_nodes = list(set(RRS) - set(RRS0))
            
            RRS0 = RRS[:]
            
        spread.append(len(RRS))

    return np.mean(spread)


def greedy(G, S, k):
    R = []
    r = []
    A = []
    s = []
    influence = []
    l = len(S)
    for m in range(1, k+1):
        for i in range(l):
            if len(S[i]) > m:
                A.append(S[i][0:m])
            else:
                A.append(S[i])
        s.append(list(set(item for subset in A for item in subset)))
    print(len(s))
    for j in range(k):
        IM = 0
        s[j] = [b for b in s[j] if b not in R]
        print(len(s[j]))
        for e in s[j]:
            r = [_ for _ in R]
            r.append(e)
            im = triggering(G, r, 0.05, 100)
            if im > IM:
                IM = im
                E = e
        print(str(j)+"个节点优化最大影响力"+str(IM))
        influence.append(IM)
        R.append(E)
        print(R)
    print(influence)

    return R
    
    
if __name__ == '__main__':

    open('higgs-social_network_134574.txt').read()
    G = Graph.Read_Edgelist('higgs-social_network_134574.txt', directed=True)
    S = [[4400, 1458, 17915, 718, 6338, 15594, 638, 9796, 9514, 6218, 1450, 3433, 731, 9646, 1353, 40849, 2054, 18278, 6217, 5578, 3430, 6212, 2025, 18523, 21752, 4306, 15592, 2465, 29383, 1105, 1946, 12839, 16911, 2470, 747, 18988, 2458, 9310, 28446, 14906, 13562, 5931, 1851, 2090, 9639, 10549, 1472, 5098, 63338, 16668], [8312, 31011, 8700, 1458, 8764, 17915, 718, 17554, 20912, 12824, 53162, 22450, 13566, 5065, 22359, 9514, 12314, 6218, 1450, 5813, 20128, 2436, 39752, 6982, 20923, 19193, 39108, 23741, 7542, 3433, 2243, 23045, 5965, 29620, 438, 731, 39688, 9646, 56666, 22199, 4317, 938, 21229, 26122, 37642, 13612, 132817, 1353, 19678, 31016], [4400, 731, 18278, 3269, 18523, 28446, 14906, 36827, 41830, 116856, 37216, 57029, 73445, 22559, 34728, 83323, 121844, 110792, 67940, 47310, 59703, 83477, 64126, 126654, 113044, 44378, 12229, 66978, 131810, 84488, 43074, 71087, 19119, 23505, 7393, 64014, 62587, 121180, 94873, 24745, 117742, 115751, 109029, 71486, 9726, 77736, 65619, 13716, 74585, 98522], [8312, 31011, 17915, 17554, 8802, 20128, 23741, 39631, 2243, 18577, 731, 132817, 17495, 28314, 11542, 13732, 108923, 4298, 8643, 47184, 2251, 11927, 94542, 67153, 47735, 121844, 46282, 110792, 42619, 14906, 13716, 115751, 32811, 39028, 74585, 76982, 21909, 48197, 92221, 67003, 66978, 50884, 81676, 57898, 9390, 30101, 67940, 22363, 108352, 17128]]
    R = greedy(G, S, 20)
    print(R)