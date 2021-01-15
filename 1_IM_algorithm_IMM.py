import numpy as np
from igraph import *
from collections import Counter
from scipy.special import comb
from tqdm import tqdm


def get_RRS(G, p, l):                           #获取一个反向可达集
    source = np.random.randint(0, l)            #随机选点ID
    #print('source:' + str(source))
    new_nodes, RRS0 = [source], [source]
    #print("进入while循环")
    while new_nodes:                            #求所选节点的反向可达集
        #print(len(new_nodes))
        B = G.neighborhood(vertices=new_nodes, order=1, mode="in", mindist=0)   #新节点的邻居节点集合
        #print("B:"+str(B))
        #print("B len:"+str(len(B)))

        A = [item for subset in B for item in subset] #新节点的邻居节点集合，转化为一个列表
        #print("A len:"+str(len(A)))

        success = np.random.uniform(0, 1, len(A)) < p
        #print(success)
        #print("success len:" + str(len(success)))

        temp = list(np.extract(success, A))           #新节点以概率p随机激活邻居节点，每个节点只能激活一次
        #print(temp)
        #print("temp len:" + str(len(temp)))

        RRS = list(set(RRS0 + temp))      #list(set())可以去除相同元素
        #print("RRS len:" + str(len(RRS)))

        new_nodes = list(set(RRS) - set(RRS0))
        RRS0 = RRS[:]

    #print('RRS:' + str(RRS))
    #print("反向可达集长度" + str(len(RRS)))
    #print("结束while循环")
    return (RRS)

#def RIS(G, k, p, mc, l, S, times):
def RIS(G, k, p, mc, l, R):                                 #选取一次种子集合
    #R = [get_RRS(G, p, l) for _ in tqdm(range(mc))]        #获得指定个数的反向可达集
    if mc>len(R):
        for _ in tqdm(range(mc - len(R))):
            R.append(get_RRS(G, p, l))
    SEED = []
    T = []
    #print('R:' + str(R))

    for _ in range(k):                                     #根据贪心策略选取k个节点，使覆盖反向可达集最多
        flat_map = [item for subset in R for item in subset]
        #print('flat_map:' + str(flat_map))
        seed = Counter(flat_map).most_common()[0][0]
        t = Counter(flat_map).most_common()[0]
        SEED.append(seed)
        T.append(t)
        #j = 1
        #while seed in S:
        #    seed = Counter(flat_map).most_common()[j][0]
        #    t = Counter(flat_map).most_common()[j]
        #    j = j + 1
        #S.append(seed)
        #times.append(t)

        R = [rrs for rrs in R if seed not in rrs]
    #print('seed:' + str(seed))
    #print(S)
    #print(times)
    return SEED


def get_mc(G, k, p, l, Rtemp):                #求对于网络G，IC模型下，在传播概率为p时，选取k个影响力最大节点，需要生成的反向可达集个数
    LB = 1
    q = 0.1                            #利用贪心算法得到的结果，与最优结果的偏差指数q
    n = int(math.log(l-1, 2))
    epsilon = math.sqrt(2) * q
    e = math.e
    c = comb(l, k)
    lamb = (2 + 2/3*epsilon) * (math.log(c, e) + math.log(l, e) + math.log(math.log(l, 2), e)) * l / (epsilon*epsilon)
    setn = 0
    for i in range(1, n+1):            #经估算，传播概率为0.05时，i = 4时求得的反向可达集个数可达到要求，所以直接从4开始循环节省时间
        x = l/(2**i)
        setN = int(lamb/x) + 1
        print("测试反向可达集总数：" + str(setN))
        #Rtemp = [get_RRS(G, p, l) for _ in tqdm(range(setN))]   # 耗时，可以继承前一轮循环的列表节省时间，由于时间成本尚可这里没有继承
        for _ in tqdm(range(setN-setn)):
            Rtemp.append(get_RRS(G, p, l))
        print("Rtemp长：" + str(len(Rtemp)))
        print("测试RR sets 创建完成")
        FR = 0
        Rt = Rtemp
        for _ in range(k):             # 根据贪心策略选取k个节点，使覆盖反向可达集最多，并计算k个节点的影响力
            flat_map = [item for subset in Rt for item in subset]
            seed = Counter(flat_map).most_common()[0][0]
            f = Counter(flat_map).most_common()[0][1]
            FR = FR + f
            Rt = [rrs for rrs in Rt if seed not in rrs]
        FR = FR/setN
        print("i等于" + str(i) + "贪心子集影响力:" + str(FR))
        if (l * FR) >= ((1 + epsilon) * x):           # 判断影响力是否能够满足阈值
            LB = l * FR/(1 + epsilon)
            break
        setn = setN

    #求解需要反向可达集个数num
    a = math.sqrt(math.log(l, e) + math.log(2, e))
    b = math.sqrt((1 - 1/e) * (math.log(l, e) + math.log(2, e) + math.log(c, e)))
    m = ((1 - 1/e)*a + b)**2
    Lambda = 2*l*m / (q*q)
    num = int(Lambda/LB) + 1
    return num


if __name__ == '__main__':
    open('higgs-social_network_134574.txt').read()
    G = Graph.Read_Edgelist('higgs-social_network_134574.txt', directed=True)
    L = len(G.vs)

    MC = []
    S = []
    for i in range(1, 21):                             #按种子集合体量为1到20，选取20次种子集合
        R = []
        mc = get_mc(G, i, 0.01, L, R)                  #返回需要生成反向可达集数目
        print("R长："+str(len(R)))
        print(str(i)+"个节点需要反向可达集子集数:"+str(mc))
        m = RIS(G, i, 0.01, mc, L, R)                  #返回所选种子集合
        S.append(m)
        MC.append(mc)
    print("种子集合"+str(S))
    print("需要反向可达集子集数：" + str(MC))

