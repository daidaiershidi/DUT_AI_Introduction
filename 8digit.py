import numpy as np

toset = {0:[1, 3], 1:[0, 2, 4], 2:[1, 5],
         3:[0,4,6], 4:[1,3,5,7], 5:[2,4,8],
         6:[3,7],  7:[4,6,8], 8:[5,7]}
'''
根据空格位置的不同，区别出下一步所能移动的位置
例如，九个数的序号分别为0-8
0 1 2
3 4 5
6 7 8
若空格在0，下一步所能移动的位置为1,3
'''

def swap_num(state, i, j, g, end):
    '''
    state:当前状态
    i:空格索引序号
    j:要交换的索引序号
    g:新状态的G
    end:目的状态
    '''
    if(i > j):
        i, j = j, i
    next_state = state[:i] + state[j] + state[i+1:j] + state[i] + state[j+1:]#交换，得到新状态
    F = get_H(next_state, end) + g#得到F
    return next_state, F

def get_H(state, end):
    '''
    H函数
    state是当前状态
    end是目标状态
    '''
    H = 0

    # 值与序号的差的总和
    # index0 = state.index("0")
    # for i in range(0, 9):
    #     if(i!=index0):
    #         H = H + abs(i - end.index(state[i]))
    # return H

    # 曼哈顿距离
    def get_site(value):
        for i in range(3):
            for j in range(3):
                if(end_arr[i][j] == value):
                    return i, j
    state_arr = np.array(list(state)).reshape((3,3))
    end_arr = np.array(list(end)).reshape((3,3))
    for i in range(3):
        for j in range(3):
            if(state_arr[i][j]!=0):
                x, y = get_site(state_arr[i][j])
                H = H + abs(x-i) + abs(y-j)
    return H

def A_star(start, end):
    '''
    A*算法
    start表示初始状态，end表示目的状态
    '''
    # 先进行判断srcLayout和end_str逆序值是否同是奇数或偶数
    # 这是判断起始状态是否能够到达目标状态，同奇同偶时才是可达 
    start_i = 0
    end_i = 0
    for i in range(0, 9):
        t = 0
        for j in range(0, i):
            if(start[i]!='0' and start[j]>start[i]):
                t = t + 1
        start_i = start_i + t
    for i in range(0, 9):
        t = 0
        for j in range(0, i):
            if(end[i]!='0' and end[j]>end[i]):
                t = t + 1
        end_i = end_i + t
    if(start_i%2 != end_i%2):
        return -1, -1
    
    #数据结构
    Father = {}
    Open = []
    Open_G = {}
    Open_F = {}
    Close = []
    Close_G = {}
    Close_F = {}
    '''
    Father相当于储存了一个所以已经搜索到的状态的树状图，这个树状图open和close通用
    Father = {状态1:-1, 状态2:状态1, 状态3:状态1}
    Open分为三个部分,分别储存每一个状态及其对应的G和F
    Open = [状态1，状态2， 状态3]
    Open_G = {状态1:0，状态2:1， 状态3:1}
    Open_F = {状态1:12，状态2:12， 状态3:13}
    Close分为三个部分,分别储存每一个状态及其对应的G和F
    Close = [状态1，状态2， 状态3]
    Close_G = {状态1:0，状态2:1， 状态3:1}
    Close_F = {状态1:12，状态2:12， 状态3:13}    
    '''

    # start存入全局树状图
    Father[start] = -1
    # 状态1:start进入open
    Open_G[start] = 0
    Open_F[start] = Open_G[start] + get_H(start, end)
    Open.append(start)

    # open≠[],即open长度不等于0,进入循环
    while(len(Open) > 0):
        # 从open表取出F最小的状态，称之为cur_state当前状态
        cur_state = min(Open_F, key=Open_F.get)# 返回最小值的key
        G = Open_G[cur_state] + 1# 取得当前的G
        
        # 存入close
        Close.append(cur_state)
        Close_G[cur_state] = Open_G[cur_state]
        Close_F[cur_state] = Open_F[cur_state]

        # 从open表中删除第一个状态(其实是每次删除最小,避免了排序)
        del Open_F[cur_state]
        del Open_G[cur_state]
        Open.remove(cur_state)

        # 判断最小状态是否是目的状态
        if(cur_state == end):
            break

        # 返回空格（0）的索引序号,并得到下一次空格所能移动的位置
        index0 = cur_state.index("0")
        slides = toset[index0]

        # 根据下一次空格所能移动的位置，生成子状态和对应F（每个状态和G不需要生成,已知）
        for slide in slides:
            new_state, F = swap_num(cur_state, index0, slide, G, end)
            
            # 如果新状态不在open和close，加入open
            if(Open_F.get(new_state) == None and Close_F.get(new_state) == None):
                # 加入open
                Open_G[new_state] = G
                Open_F[new_state] = F
                Open.append(new_state)
                # 存入全局树状图
                Father[new_state] = cur_state
            # 如果新状态在open，并且G小于open中的值
            elif(Open_F.get(new_state) != None and G < Open_G[new_state]):
                # 更新open
                Open_F[new_state] = F
                # 存入全局树状图
                Father[new_state] = cur_state
            # 如果新状态在close，并且G小于close中的值
            elif(Close_F.get(new_state) != None and G < Close_G[new_state]):
                # 从close中取出
                del Close_F[new_state]
                del Close_G[new_state]
                Close.remove(new_state)
                # 将新的值放入open
                Open.append(new_state)
                Open_G[new_state] = G
                Open_F[new_state] = F
                # 存入全局树状图
                Father[new_state] = cur_state              

    #通过目的状态寻找上一个状态，整理路径存入一个列表
    how2end = []
    how2end.append(cur_state)
    while(Father[cur_state] != -1):
        cur_state = Father[cur_state]
        how2end.append(cur_state)
    how2end.reverse()#将路径反过来
    return 0, how2end
if __name__ == "__main__":
    '''
    start_str = '123456780'表示初始状态为
    1 2 3
    4 5 6
    7 8 0（空格）
    end_str目标状态也同理
    '''
    start_str  = "013726548"
    end_str = "123647850"
    retCode, how2end = A_star(start_str, end_str)
    if retCode != 0:
        print("目标不可达")
    else:
        for i in range(len(how2end)):
            print("第" + str(i + 1) + '步')
            arr = np.array(list(how2end[i])).reshape((3,3))#字符串变数组
            print(arr)
