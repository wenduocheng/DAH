def StringReconstruction(file):
    with open("week3/dataset_635620_8.txt") as file:
        a = int(file.readline())
        Dna = file.readlines()
        Dna = [line.rstrip() for line in Dna]

    dB = KmerDeBruijin(Dna)

    with open('dB.txt', 'w') as f:
        for k in dB.keys():
            f.write(k + ' -> ')
            f.write(','.join(dB[k]))
            f.write('\n')

    input = {}
    with open ('dB.txt', 'r') as f:
        lines = f.readlines()

        for line in lines:
            l = line.rstrip()
            sliced = l.split(' -> ')
            input[sliced[0]] = []

            nodes = sliced[1].split(',')

            input[sliced[0]] = nodes

    path = EulerianCycle(input)
    text = SimpleReconstruction(path)
    return text




#Debruijin
def KmerDeBruijin(Dna):
    num = 0
    result = {}
    k = len(Dna[0])
    for num in range(0,len(Dna)):
        string = Dna[num]
        if string[0:k-1] not in result:
            result[string[0:k-1]] = []
        result[string[0:k-1]].append(string[1:k])
        num += 1

    return result


##################################################
#EulerianPath
def empty_edge(d):
    for k in d.values():
        if not (k == []):
            return False
    return True

def dfs(graph, startV):
    path = [startV]
    while not (graph[startV] == []):
        path.append(graph[startV][0])
        startV = graph[startV].pop(0)
    return graph, path

def select_startV(graph, path):
    for node in path:
        if not graph[node] == []:
            for edge in graph.values():
                if node in edge:
                    return node
    return False

def findStartEnd(graph):

    start = 0
    end = 0
    for node in graph.keys():
        numIn = 0
        numOut = 0
        numOut = len(graph[node])
        for node2 in graph.keys():
            if node in graph[node2]:
                numIn +=1
        if numOut > numIn:
            start = node
        if numOut < numIn:
            end = node
    for node in graph.keys():
        for x in graph[node]:
            if x not in graph.keys():
                end = x

    return start,end

def EulerianCycle(graph):

    startV,endV = findStartEnd(graph)
    startVcopy = startV

    graph[endV] = [startV]
    path = []

    curr_graph = graph.copy()
    while not empty_edge(curr_graph):
        curr_graph, found_path = dfs(curr_graph, startV)

        if path == []:
            path = found_path
        else:
            path = found_path[0:-1] + path

        new_start = select_startV(curr_graph, path)

        if not new_start:
            break

        x = path.index(new_start)

        # while x >= 0:
        #     if not curr_graph[path[x]] == []:
        #         break
        #     x -= 1

        path = (path[x:] + path[1:x])
        path.append(path[0])
        startV = path[0]

    #pop the end one to avoid repeat
    path.pop()
    num = 0
    for x in range(0,len(path)-1):
        if path[x] == endV and path[x+1] == startVcopy:
            num = x+1
            break

    path = path[num::]+path[0:num]
    return path


def SimpleReconstruction(Dna):
    Dna_list = Dna
    string = ""
    k = len(Dna_list[0])
    for n in range(0,len(Dna_list)):
        if (n == 0):
            string = string + Dna_list[n]
        else:
            string = string + Dna_list[n][k-1]

    return string



file = "week3/dataset_635620_8.txt"
print(StringReconstruction(file))
