import numpy as np
from igraph import *
import time
from tqdm import *#进度条显示
import matplotlib.pyplot as plt

# 对字典按值排序，返回元祖列表
def ST(d):
    return sorted(d.items(), key = lambda x:x[1], reverse=True)

#CCA核覆盖算法
def CCA(g, k_seeds, d):
    s_time = time.time()
    #种子节点
    s = []
    #节点索引值
    V = range(1, g.vcount()+1)
    #核数字典，迭代删减
    t_index = dict(zip(V, t[1:]))
    #出度字典，迭代删减
    o_index = dict(zip(V, o[1:]))
    #节点覆盖度
    CO = np.zeros(g.vcount()+1)
    for _ in range(k_seeds):
        tmp = []#候选列表
        Max = ST(t_index)[0][1]#当前最大核数
        for k in t_index.keys():
            if t_index[k] == Max and CO[k] == 0:
                tmp.append(k)
        #print(tmp)
        D = {}#最大核数对应出度字典
        for y in tmp:
            if CO[y] == 0 and y in o_index.keys():
                D[y] = o_index[y]
        #print(D)
        seed = ST(D)[0][0]#当前轮候选种子节点
        #print(seed)
        s.append(seed)
        #CCA1
        # for v in g.neighbors(seed, mode=OUT):
        #         CO[v] = 1
        #CCA2
        # for v1 in g.neighbors(seed, mode=OUT):
        #     for v2 in g.neighbors(v1, mode=OUT):
        #         CO[v2] = 1
        for v in g.neighborhood(seed, order=d, mode=OUT):#当前被覆盖d跳数重置覆盖指示1，且被覆盖的节点退出相应字典
                CO[v] = 1
                if v in t_index.keys():
                    t_index.pop(v)
                if v in o_index.keys():
                    o_index.pop(v)
    #print(sum(CO))#打印当前覆盖范围
    return [s, time.time()-s_time]

def DCA(g, k_seeds, d):
    s_time = time.time()
    #种子节点
    s = []
    #节点索引值
    V = range(1, g.vcount()+1)
    #核数字典，迭代删减
    t_index = dict(zip(V, t[1:]))
    #出度字典，迭代删减
    o_index = dict(zip(V, o[1:]))
    #节点覆盖度
    CO = np.zeros(g.vcount()+1)
    #for _ in tqdm(range(k_seeds)):
    for _ in range(k_seeds):
        tmp = []#候选列表
        Max = ST(o_index)[0][1]#当前最大度数
        for k in o_index.keys():
            if o_index[k] == Max and CO[k] == 0:
                tmp.append(k)
        #print(tmp)
        D = {}#最大度数对应出度字典
        for y in tmp:
            if CO[y] == 0 and y in t_index.keys():
                D[y] = t_index[y]
        #print(D)
        seed = ST(D)[0][0]#当前轮候选种子节点
        #print(seed)
        s.append(seed)
        #CCA1
        # for v in g.neighbors(seed, mode=OUT):
        #         CO[v] = 1
        #CCA2
        # for v1 in g.neighbors(seed, mode=OUT):
        #     for v2 in g.neighbors(v1, mode=OUT):
        #         CO[v2] = 1
        for v in g.neighborhood(seed, order=d, mode=OUT):#当前被覆盖d跳数重置覆盖指示1，且被覆盖的节点退出相应字典
                CO[v] = 1
                if v in o_index.keys():
                    o_index.pop(v)
                if v in t_index.keys():
                    t_index.pop(v)
    #print(sum(CO))#打印当前覆盖范围
    return [s, time.time()-s_time]

def averRuntime_CCA(G, S, d, rep):
    timLis = []
    for _ in range(rep):
        timLis.append(CCA(G, S, d)[1])
    return np.mean(timLis)

def averRuntime_DCA(G, S, d, rep):
    timLis = []
    for _ in range(rep):
        timLis.append(DCA(G, S, d)[1])
    return np.mean(timLis)

if __name__ == '__main__':
    g = Graph.Read_Edgelist("C:\\Users\\Snail\\Desktop\\higgs\\higgs-social_network_134574.txt", directed=True)
    # 计算核数存储在t
    t = g.coreness()
    # 计算入度存储在i
    i = g.indegree()
    # 计算出度存储在o
    o = g.outdegree()
    k_l = range(1, 51)
    t_CCA1 = []
    t_CCA2 = []
    t_DCA1 = []
    t_DCA2 = []
    rep = 100
    for k in tqdm(k_l):
        # print(k)
        t_CCA1.append(averRuntime_CCA(g, k, 1, rep))
        t_CCA2.append(averRuntime_CCA(g, k, 2, rep))
        t_DCA1.append(averRuntime_DCA(g, k, 1, rep))
        t_DCA2.append(averRuntime_DCA(g, k, 2, rep))
    plt.figure(1)
    plt.plot(k_l, t_CCA1, label='CCA1', marker='o', color='r', linestyle=':')
    plt.plot(k_l, t_CCA2, label='CCA2', marker='+', color='y', linestyle='-')
    plt.plot(k_l, t_DCA1, label='DCA1', marker='*', color='b', linestyle='--')
    plt.plot(k_l, t_DCA2, label='DCA2', marker='D', color='g', linestyle='-.')
    plt.legend()
    plt.xlabel('seed_set size')
    plt.ylabel('time (s)')
    plt.title('time (s) vs seed_set size')
    plt.show()