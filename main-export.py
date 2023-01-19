# %%
from collections import defaultdict
from collections import Counter
import random
import numpy as np
import networkx as nx
import os
import matplotlib.pyplot as plt
import random as rd
import statistics as st

# %%
dataset = 'data2.txt'

# %%
# def get_edge_list(file_name):
#     with open(file_name) as ifs:
#         lines = ifs.readlines()
#         edge_list = map(lambda line: map(int, line.strip().split()), filter(lambda ele: '#' not in ele, lines))
#     return edge_list


def get_edge_list(file_name):
    # Open the file and read its contents
    with open(file_name) as file:
        lines = file.readlines()

    # Remove lines that are comments (start with "#")
    lines = [line for line in lines if '#' not in line]

    # Convert each line to a pair of integers
    edge_list = []
    for line in lines:
        edge = tuple(map(int, line.strip().split()))
        edge_list.append(edge)

    return edge_list


def get_undirected_graph_info(file_name):
    my_edge_list = get_edge_list(file_name)
    undirected_graph = nx.Graph()
    undirected_graph.add_edges_from(my_edge_list)
    file_info = file_name.split(os.sep)[-1].split('_')[0]
    str_list = [file_info, 'nodes num:' + str(undirected_graph.number_of_nodes()), 'edges num:'
                + str(undirected_graph.number_of_edges())]
    print(' | '.join(str_list))


def get_dir_info(dir_name):
    my_walk = os.walk(dir_name)
    my_root, sub_root_list, file_list = list(my_walk)[0]
    path_list = map(lambda ele: my_root + os.sep + ele, file_list)
    for my_path in path_list:
        get_undirected_graph_info(my_path)


# %%
get_undirected_graph_info(dataset)

# %%
edgeList = get_edge_list(dataset)

# %%
# print(edgeList)

# %%
undirected_graph = nx.Graph()
undirected_graph.add_edges_from(edgeList)

# %%

# nx.draw_networkx(undirected_graph,pos=nx.random_layout(undirected_graph))

# %%
# rd.choice (list(undirected_graph.edges(1)))
# undirected_graph.edges((1,2))


# %%
graphEdge = undirected_graph.edges()
graphNode = undirected_graph.nodes()


# %%
number_of_edges = len(graphEdge)
number_of_nodes = len(graphNode)
print(number_of_nodes, number_of_edges)

# %%
initialGene = np.full(number_of_edges, -1)

# %%
locus = [i for i in range(1, number_of_edges+1)]

# %%
len(initialGene)

# %%
labledEdge = {t: i+1 for i, t in enumerate(graphEdge)}
# labledEdge

# %%


def getEdgeLabel(labledEdge, edge):
    try:
        return labledEdge[edge]
    except:
        try:
            return labledEdge[(edge[1], edge[0])]
        except:
            print('NO KEY FOUND')


# %%
def getEdgeValue(labledEdge, value):
    return next((k for k, v in labledEdge.items() if v == value), None)


getEdgeValue(labledEdge, 5)

# %%
# make initial communities

gene = []
for edge, l in labledEdge.items():
    connectedEdges = list(undirected_graph.edges(edge))
    geneValue = getEdgeLabel(labledEdge, rd.choice(connectedEdges))

    if (geneValue != None):
        gene.append(geneValue)

# print(locus)
# print(gene)


# %%
def drawLayout(labledEdge):
    # layout = nx.random_layout(undirected_graph)
    layout = nx.circular_layout(undirected_graph)
    plt.figure(3, figsize=(12, 12))
    nx.draw_networkx(undirected_graph, pos=layout, node_size=500)
    nx.draw_networkx_edge_labels(
        undirected_graph, pos=layout, edge_labels=labledEdge)
# drawLayout(labledEdge)

# %%
# nx.draw_networkx_edge_labels(undirected_graph,pos=nx.random_layout(undirected_graph),edge_labels=labledEdge)

# %%
# finished but not tested with real gene


def calculateFitness(communities: list, undirected_graph):
    D = 0
    E = undirected_graph.number_of_edges()
    for i in communities:
        community = undirected_graph.subgraph(i)
        mc = community.number_of_edges()
        nc = community.number_of_nodes()
        if nc >= 3:
            D += mc * ((mc - (nc-1)) / (((nc-2)*(nc-1))))
            # print(D)

    return 2*D/E

# calculateFitness([[1, 5], [2, 3, 4], [6, 8, 9], [7]],undirected_graph)


# %%


def most_frequent(lst):
    # Count the number of occurrences of each element in the list
    count = Counter(lst)

    # Get a list of tuples with the elements and their frequencies, sorted by frequency in descending order
    most_common = count.most_common()

    # Select the most frequent item and, if there are multiple items with the same frequency, choose one at random
    selected = None
    frequency = 0
    for element, freq in most_common:
        if freq > frequency:
            selected = random.sample([element], 1)[0]
            frequency = freq
        elif freq == frequency:
            selected = random.sample([selected, element], 1)[0]

    return selected


# Test the function
print(most_frequent([1, 2, 3, 3, 3, 4, 4, 5, 5]))  # 3
print(most_frequent([1, 2, 3, 3, 4, 4, 5, 5]))  # 3 or 5
print(most_frequent([1, 2, 3, 4, 5]))  # 1 or 2 or 3 or 4 or 5

# def most_common(List):
#     return(st.mode(List))

# List = [2, 1, 2, 2, 1, 3]
# print(most_common(List))

# %%


def has_repetitions(lst):
    return len(set(lst)) < len(lst)


has_repetitions([1, 1, 2, 2])


# %%


def createCommunities(x: dict):
    groups = defaultdict(list)

    # Group the keys based on the values
    for key, value in x.items():
        groups[value].append(key)

    # Get the list of lists with the keys
    y = list(groups.values())

    return y


x = {1: 2, 2: 2, 3: 2, 4: 8, 5: 1, 6: 8, 7: 8, 8: 8, 9: 8}
createCommunities(x)

# %%
# Label Propagation


def labelProp(undirected_graph):

    graphNode = undirected_graph.nodes()
    asList = sorted(list(graphNode))
    randomizedList = asList[:]
    rd.shuffle(randomizedList)
    # print(randomizedList)
    labledNodes = {i+1: t for i, t in enumerate(asList)}
    prevLabledNodes = labledNodes
    # print(labledNodes)
    stop = False
    while not stop:
        for v in randomizedList:
            # v = rd.choice(asList)
            # v = 5
            neighbors = [i for i in nx.all_neighbors(undirected_graph, v)]
            labledNeighbors = [labledNodes[i]
                               for i in nx.all_neighbors(undirected_graph, v)]
            # print(neighbors)

            # if(has_repetitions(labledNeighbors)):
            #     freq = most_frequent(labledNeighbors)
            #     labledNodes[v] = freq
            # else:
            #     freq = rd.choice(labledNeighbors)
            #     labledNodes[v] = freq

            # print(labledNodes)
            freq = most_frequent(labledNeighbors)
            labledNodes[v] = freq

            if (prevLabledNodes == labledNodes):
                stop = True
    return labledNodes


createCommunities(labelProp(undirected_graph))

# %%
# Local Expansion


def localExpansion(undirected_graph):
    counter = 0
    communities = []
    graphNode = undirected_graph.nodes()
    asList = list(graphNode)
    V = len(asList)
    tempList = asList[:]
    while counter < V:
        community = []
        seed = rd.choice(tempList)
        community.append(seed)
        tempList.remove(seed)
        counter += 1
        naturalCommNodes = list(nx.all_neighbors(undirected_graph, seed))
        for i in naturalCommNodes:
            # naturalCommNodes.remove(seed)
            if i not in tempList:
                continue
            community.append(i)
            try:
                tempList.remove(i)
            except ValueError:
                print('E')
            counter += 1
        communities.append(community)
    return communities

# localExpansion(undirected_graph)


# %%
centrality = nx.eigenvector_centrality(undirected_graph)
sorted_centrality = sorted(
    centrality.items(), key=lambda x: x[1], reverse=True)
sorted_centrality[0:11]

# %%
# Local Expansion


def localExpansionWithEigen(undirected_graph):
    counter = 0
    communities = []
    graphNode = undirected_graph.nodes()
    asList = list(graphNode)
    V = len(asList)
    tempList = asList[:]

    centrality = nx.eigenvector_centrality(undirected_graph)
    sorted_centrality = sorted(
        centrality.items(), key=lambda x: x[1], reverse=True)

    while counter < V:
        community = []
        seed = sorted_centrality[0][0]
        sorted_centrality.pop(0)
        if seed not in tempList:
            continue
        community.append(seed)
        tempList.remove(seed)
        counter += 1
        naturalCommNodes = list(nx.all_neighbors(undirected_graph, seed))
        for i in naturalCommNodes:
            # naturalCommNodes.remove(seed)
            if i not in tempList:
                continue
            community.append(i)
            tempList.remove(i)
            counter += 1
        communities.append(community)
    return communities

# localExpansionWithEigen(undirected_graph)

# %%


for i in range(7):
    le = localExpansion(undirected_graph)
    fit = calculateFitness(le, undirected_graph)
    # print(f'fitness : {fit} | CommunityList : {le}')


# %%
for i in range(7):
    lp = labelProp(undirected_graph)
    cmnti = createCommunities(lp)
    # print(cmnti)
    fit = calculateFitness(cmnti, undirected_graph)
    # print(f'fitness : {fit} | CommunityList : {cmnti}')


# %%
# for i in range(7):
#     le = localExpansionWithEigen(undirected_graph)
#     fit = calculateFitness(le,undirected_graph)
#     print(f'fitness : {fit} | CommunityList : {le}')

# %%
# undirected_graph.edges([8,9])

# %%
g = [-1]*5

g[4] = 'a'
g

# %%


def has_common_element(list1, list2):
    for element in list1:
        if element in list2:
            return True
    return False

# %%
# # this takes a list of community and returns a list of locus and gene
# def communityToGene(communities:list,number_of_nodes:int,number_of_edges:int,labledEdge:dict):
#     locus = [i for i in range(1,number_of_edges+1)]
#     gene = [-1] * number_of_edges
#     # overlappingEdges = []
#     for community in communities:
#         edgesAsTuples = list(undirected_graph.edges(community))
#         edges=[getEdgeLabel(labledEdge,i) for i in edgesAsTuples]
#         # print(edgesAsTuples)
#         for n,i in enumerate(edgesAsTuples):
#             edgeLabel = getEdgeLabel(labledEdge,i)
#             ngEdges = list(undirected_graph.edges(i))
#             ngEdges.remove(i)
#             filtered_ngEdges = [t for t in ngEdges if t[0] in community and t[1] in community]
#             ngEdgesLabeled=[getEdgeLabel(labledEdge,t) for t in filtered_ngEdges]
#             # print(ngEdgesLabeled)
#             for element in ngEdgesLabeled:
#                 if element in edges:
#                     if(gene[edgeLabel-1]!=-1 or gene[element-1]== edgeLabel):
#                         continue
#                     gene[edgeLabel-1] = element
#                     edges.remove(element)
#                     break
#     # for n,i in enumerate(gene):
#     #     if i == -1:
#     #         l= locus[n]
#     #         edge = getEdgeValue(labledEdge,l)
#     #         for k in communities:
#     #             if edge[0] in k or edge[1] in k:


#     return locus,gene


# # locus,gene = [],[]
# # locus,gene = communityToGene([[3, 1, 2, 4], [8, 6, 7, 9], [5]],number_of_nodes,number_of_edges,labledEdge)
# locus,gene = communityToGene([[1,2,3,4,5], [5,6,7,8,9]],number_of_nodes,number_of_edges,labledEdge)
# print((locus))
# print((gene))

# # locus,gene = communityToGene(localExpansionWithEigen(undirected_graph),number_of_nodes,number_of_edges,labledEdge)
# # print((gene))
# # decode_communities(locus,gene)


# %%
def removeEdge(adj_list, i):
    for edge in adj_list:
        try:
            adj_list[edge].remove(i)
        except:
            pass
    return adj_list


# %%

def communityToGene(communities: list, number_of_nodes: int, number_of_edges: int, labledEdge: dict):
    locus = [i for i in range(1, number_of_edges+1)]
    gene = [-1] * number_of_edges
    # overlappingEdges = []
    notVisited = set()
    for community in communities:
        subGraph = undirected_graph.subgraph(community)
        subGraphEdges = list(subGraph.edges())
        # print(subGraphEdges)
        adj_list = defaultdict(set)
        adL = defaultdict(set)
        for i in subGraphEdges:
            nbr = list(subGraph.edges(i))
            nbr.remove(i)
            for j in nbr:
                adj_list[getEdgeLabel(labledEdge, i)].add(
                    getEdgeLabel(labledEdge, j))
                adL[getEdgeLabel(labledEdge, i)].add(
                    getEdgeLabel(labledEdge, j))
        old_i = -1

        for edge in adj_list:
            for i in adj_list[edge]:
                if old_i != edge and gene[edge-1] == -1:
                    # print(i)
                    adj_list = removeEdge(adj_list, i)
                    gene[edge-1] = i
                    old_i = i
                    break
                else:
                    notVisited.add(edge)
        # print(notVisited)
        # print(adL)
        while notVisited:
            value = notVisited.pop()
            for i in adL[value]:
                # print(locus[gene.index(value)])
                if locus[gene.index(value)] != i:
                    gene[value-1] = i
                    break
        for i in range(len(gene)):
            if gene[i] == -1:
                gene[i] = locus[i]

    return locus, gene

# locus,gene = communityToGene([[3, 1, 2, 4], [8, 6, 7, 9], [5]],number_of_nodes,number_of_edges,labledEdge)
# print((locus))
# print((gene))


# %%


# %%
def locusNgeneToTuples(locus: list, gene: list):
    pass

# %%


def crossOver(parent1: list, parent2: list, cp: float):
    offspring = [-1]*len(parent1)
    for g in range(len(offspring)):
        if rd.uniform(0, 1) < cp:
            offspring[g] = parent1[g]
        else:
            offspring[g] = parent2[g]
    return offspring

# testParent1 = [2, 3, 1, 5, 6, 4, 8, 9, 7, 8, 10, 11, 9, 13]
# testParent2 = [5, 1, 1, 6, 2, 7, 4, 9, 10, 8, 12, 14, 11, 11]
# crossOver(testParent1,testParent2, 0.5)


# %%
# drawLayout(labledEdge)

# %%
def mutate(gene: list, mp: float, labledEdge: dict, undirected_graph) -> list:
    individual = gene[:]
    # print(individual)
    for g in range(len(individual)):
        if rd.uniform(0, 1) < mp:
            edge = getEdgeValue(labledEdge, g+1)
            # print(edge)
            neigEdges = list(undirected_graph.edges(edge))
            neigEdgesLabeled = [getEdgeLabel(labledEdge, i) for i in neigEdges]
            # oldGene = getEdgeValue( labledEdge ,individual[g])
            # print("In Mutate",individual[g])
            try:
                neigEdgesLabeled.remove(individual[g])
                neigEdgesLabeled.remove(g+1)
            except ValueError:
                print('Ignoring ValueError In mutate function for ',
                      individual[g])
            if ((neigEdgesLabeled) == 0):
                continue
            mutatedEdge = rd.choice(neigEdgesLabeled)
            # print(randomEdge)
            # mutatedEdge = getEdgeLabel(labledEdge,randomEdge)
            individual[g] = mutatedEdge

    return individual


# print("Before Mutation ⬇")
# print([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
# print([5, 3, 1, 6, 6, 7, 4, 9, 7, 8, 12, 14, 9, 11])
# mutated = mutate([5, 3, 1, 6, 6, 7, 4, 9, 7, 8, 12, 14, 9, 11],0.5,labledEdge,undirected_graph)
# print([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
# print(mutated)
# print("After Mutation ⬆")

# %%
def localSearch(individual: list, labledEdge: dict, undirected_graph) -> list:
    locus, individualEdges = communityToGene(individual)
    for cluster in individual:
        for node in cluster:
            edge = getEdgeValue(labledEdge, g+1)
            # print(edge)
            neigEdges = list(undirected_graph.edges(edge))
            neigEdgesLabeled = [getEdgeLabel(labledEdge, i) for i in neigEdges]
            oldEdge = getEdgeValue(labledEdge, individual[g])
            # neigEdgesLabeled.remove(individual[g])
            # neigEdgesLabeled.remove(g+1)
            # if((neigEdgesLabeled)==0): continue
            # mutatedEdge = rd.choice(neigEdgesLabeled)
            # # print(randomEdge)
            # # mutatedEdge = getEdgeLabel(labledEdge,randomEdge)
            # individual[g]= mutatedEdge
            try:
                neigEdges.remove(oldEdge)
            except KeyError:
                neigEdges.remove((oldEdge[1], oldEdge[0]))
            for i in neigEdges:
                pass


# %%
# def decode_communities(locus, gene):
#     # Initialize an empty list to store the communities of locuss
#     communities = []

#     # Initialize a dictionary to store the mapping from locuss to their communities
#     locus_to_community = {}
#     loc = []
#     gen = []
#     s = -1
#     # Iterate through each locus and its corresponding gene value
#     for l,g  in zip(locus, gene):
#         if len(loc)==0:
#             loc.append(l)
#             gene.append(g)
#             continue

#         # If the locus is not already assigned to a community
#         if l not in gene:
#             # If the gene value for this locus is also not in the locus_to_community dictionary,
#             # this means we have found a new community of locuss
#             if g not in loc:
#                 # Create a new community and add the current locus to it
#                 new_community = [locus,community]
#                 locus_to_community[locus] = new_community
#                 communities.append(new_community)
#                 # print(locus_to_community)
#             # Otherwise, the gene value corresponds to another locus that is already in a community,
#             # so we add the current locus to that community
#             else:
#                 loc.append(l)
#                 gene.append(g)
#         # If the locus is already assigned to a community, we check if it is in the same community
#         # as the locus corresponding to its gene value. If it is not, we merge the two communities.
#         else:
#             # If the gene value corresponds to an locus that is not in the same community,
#             # we merge the two communities by adding the locuss in the community of the gene value
#             # to the community of the current locus
#             if (locus in locus_to_community) and (community in locus_to_community) and (locus_to_community[locus] != locus_to_community[community]):
#                 community_to_merge = locus_to_community[community]
#                 current_community = locus_to_community[locus]
#                 current_community.extend(community_to_merge)
#                 for locus_in_community in community_to_merge:
#                     locus_to_community[locus_in_community] = current_community
#                 # Remove the merged community from the list of communities
#                 communities.remove(community_to_merge)

#     # Return the list of communities of locuss
#     # print(locus_to_community)
#     return communities

# locus = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
# gene =  [2, 3, 6, 1, 4, 7, 5, 11, 10, 12, 12, 9, 14, 11]
# communities = decode_communities(locus, gene)
# print(communities)


# %%
# from google_play_scraper import app

# result = app(
#     'com.whatsapp'
# )
# result['genreId']

# %%
def generateAgents(undirected_graph):
    # Generate a random float between 0 and 1

    agents = []
    fitness = []
    # Call the corresponding function
    indi = []
    fit = -1
    for i in range(7):
        x = rd.randint(1, 3)

        if x == 1:
            indi = localExpansion(undirected_graph)
            fit = calculateFitness(indi, undirected_graph)
            # print(f'fitness : {fit} | CommunityList : {indi}')
        elif x == 2:
            lp = labelProp(undirected_graph)
            indi = createCommunities(lp)
            # print(cmnti)
            fit = calculateFitness(indi, undirected_graph)
            # print(f'fitness : {fit} | CommunityList : {indi}')
        else:
            indi = localExpansionWithEigen(undirected_graph)
            fit = calculateFitness(indi, undirected_graph)
            # print(f'fitness : {fit} | CommunityList : {indi}')
        agents.append(indi)
        fitness.append(fit)

    pairs = zip(fitness, agents)
    # Sort the pairs using the values in the first array as the key
    sorted_pairs = sorted(pairs, key=lambda x: x[0], reverse=True)

    # Unpack the sorted pairs back into two separate arrays
    fitness, agents = zip(*sorted_pairs)
    return list(fitness), list(agents)
    # return dict(sorted_pairs)

# generateAgents(undirected_graph)

# %%


class Agent:
    pocket = None
    pocketFitness = None
    individuals = None
    individualsFitness = None
    std = None

    def print_values(self):
        print(f'pocket = {self.pocket}')
        print(f'pocketFitness = {self.pocketFitness}')
        print(f'individuals = {self.individuals}')
        print(f'individualsFitness = {self.individualsFitness}')
        print(f'std = {self.std}')

# %%


def calculate_std(agent):
    agent.std = st.pstdev(agent.individualsFitness)
    return agent.std

# %%


def lostDiversity(agent):
    current_std = st.pstdev(agent.individualsFitness)
    if (current_std <= agent.std/2):
        return True
    else:
        return False

# %%


def generateTreeBase():
    treeBase = []
    for n in range(3):
        agentPool = []
        for i in range(3):
            agent = Agent()
            agent.individualsFitness, agent.individuals = generateAgents(
                undirected_graph)
            agent.pocket = agent.individuals[0]
            agent.pocketFitness = agent.individualsFitness[0]
            calculate_std(agent)
            agentPool.append(agent)
        treeBase.append(agentPool)
    return treeBase


treeBase = generateTreeBase()


# %%
def sort_and_store(agent):
    paired_list = list(zip(agent.individuals, agent.individualsFitness))
    sorted_paired_list = sorted(paired_list, key=lambda x: x[1], reverse=True)
    agent.individuals, agent.individualsFitness = zip(*sorted_paired_list)

    agent.pocketFitness = agent.individualsFitness[0]
    agent.pocket = agent.individuals[0]

# %%


def generateTreeMiddle(treeBase):
    treeMiddle = []
    for n in treeBase:
        agent = Agent()
        agent.individuals = []
        agent.individualsFitness = []
        for i in n:
            agent.individuals.append(i.pocket)
            agent.individualsFitness.append(i.pocketFitness)
            added = []
        while len(added) < 4:
            supporter = rd.choice(n)
            cell = supporter.individuals.index(
                rd.choice(supporter.individuals))
            if cell not in added:
                added.append(cell)
                agent.individuals.append(supporter.individuals[cell])
                agent.individualsFitness.append(
                    supporter.individualsFitness[cell])
        sort_and_store(agent)
        calculate_std(agent)
        treeMiddle.append(agent)
    return treeMiddle


treeMiddle = generateTreeMiddle(treeBase)


# %%
# x = [[6,7],[1,2,3]]

# x[0][]


# %%
def generateTreeTop(treeMiddle):
    agent = Agent()
    agent.individuals = []
    agent.individualsFitness = []
    for i in treeMiddle:
        agent.individuals.append(i.pocket)
        agent.individualsFitness.append(i.pocketFitness)
        added = []
    while len(added) < 4:
        supporter = rd.choice(treeMiddle)
        cell = supporter.individuals.index(rd.choice(supporter.individuals))
        if cell not in added:
            added.append(cell)
            agent.individuals.append(supporter.individuals[cell])
            agent.individualsFitness.append(supporter.individualsFitness[cell])
    sort_and_store(agent)
    treeMiddle.append(agent)
    return agent


treeTop = generateTreeTop(treeMiddle)

# %%
for i in treeMiddle:
    i.print_values()
    print("---------------")

# %%
treeTop.print_values()

# %%
# locus = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
# gene =  [2, 3, 6, 1, 4, 7, 5, 11, 10, 12, 12, 9, 14, 11]
# # communities = decode_communities(locus, gene)
# # print(communities)
# # [[2,3,6,1],[]]
# geneTuple=[]
# for u,v in zip(locus,gene):
#     geneTuple.append((u,v))
# geneTuple


# %%
def decode_communities(locus, gene):

    edges = []
    for u, v in zip(locus, gene):
        edges.append((u, v))

    clusters = []
    adj_list = defaultdict(set)

    for u, v in edges:

        adj_list[u].add(v)

        adj_list[v].add(u)
    # adj_list = defaultdict(set,labledEdgeR)
    # adj_list = defaultdict(labledEdge)

    visited = set()

    for node in adj_list:
        if node in visited:
            continue
        visited.add(node)
        cluster = {node}
        queue = list(adj_list[node])
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue

            visited.add(current)

            cluster.add(current)

            # if len(cluster) >= k:

            #     break

            for neighbor in adj_list[current]:

                if neighbor not in visited:

                    queue.append(neighbor)

        clusters.append(cluster)
    listCluster = []
    for i in clusters:
        listCluster.append(list(i))
    # print(listCluster)

    return listCluster
# locus,gene = communityToGene([[1,2,3,4,5], [4,6,7,8, 9]],number_of_nodes,number_of_edges,labledEdge)
# print(locus)

# gene = [2, 5, 3, 1, 4, 5, 7, 8, 9, 11, 12, 10, 10, 11]
# print(gene)
# decode_communities(locus,gene)


# %%


def decode_communities2(edges):

    clusters = []
    adj_list = defaultdict(set)

    for u, v in edges:

        adj_list[u].add(v)

        adj_list[v].add(u)
    # adj_list = defaultdict(set,labledEdgeR)
    # adj_list = defaultdict(labledEdge)

    visited = set()

    for node in adj_list:
        if node in visited:
            continue
        visited.add(node)
        cluster = {node}
        queue = list(adj_list[node])
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue

            visited.add(current)

            cluster.add(current)

            # if len(cluster) >= k:

            #     break

            for neighbor in adj_list[current]:

                if neighbor not in visited:

                    queue.append(neighbor)

        clusters.append(cluster)
    listCluster = []
    for i in clusters:
        listCluster.append(list(i))
    # print(listCluster)

    # print (adj_list)
    return listCluster
# locus,gene = communityToGene2([[1,2,3,4,5], [4,6,7,8, 9]],number_of_nodes,number_of_edges,labledEdge)
# print(locus)

# gene = [2,3,7,6,4,5,3,9,13,11,14,12,10,12]
# print(gene)
# decode_communities2([(1, 2), (1, 3), (1, 5), (2, 4), (2, 3), (3, 4), (5, 4)])

# %%
# # print(localExpansionWithEigen(undirected_graph))
# l,g = communityToGene(localExpansionWithEigen(undirected_graph),number_of_nodes,number_of_edges,labledEdge)
# print(l)
# print(g)
# geneTuple=[]
# for u,v in zip(l,g):
#     geneTuple.append((u,v))
# # decode_communities(geneTuple)

# %%
# drawLayout(labledEdge)

# %%


def initializePopulation():
    base = generateTreeBase()
    middle = generateTreeMiddle(base)
    top = generateTreeTop(middle)
    return base, middle, top

# %%


def updatePopulation(mutatedOffspring, base):
    x = rd.randint(0, 2)
    y = rd.randint(0, 2)
    baseAgent = base[x][y]
    # baseAgent.print_values()
    baseAgent.individuals = list(baseAgent.individuals)
    baseAgent.individuals[-1] = list(mutatedOffspring)
    baseAgent.individualsFitness = list(baseAgent.individualsFitness)
    baseAgent.individualsFitness[-1] = calculateFitness(
        mutatedOffspring, undirected_graph)
    sort_and_store(baseAgent)
    middle = generateTreeMiddle(base)
    top = generateTreeTop(middle)
    return base, middle, top


# %%
def M_Link():
    base, middle, top = initializePopulation()
    count = 0
    topScore = -1
    while count < 5:
        parent1 = top.pocket
        calculate_std(top)
        if not lostDiversity(top):
            parent2Agent = rd.choice(middle)
        else:
            parent2Agent = rd.choice(rd.choice(base))
        parent2 = rd.choice(parent2Agent.individuals)
        print('parents')
        print(parent1, parent2)
        locus, parent1Edge = communityToGene(
            parent1, number_of_nodes, number_of_edges, labledEdge)
        _, parent2Edge = communityToGene(
            parent2, number_of_nodes, number_of_edges, labledEdge)
        print("parentsEdge")
        print(parent1Edge, parent2Edge)

        offspring = crossOver(parent1Edge, parent2Edge, .5)
        print("offspring")
        print(offspring)

        #
        mutatedOffspringEdge = mutate(
            offspring, .5, labledEdge, undirected_graph)
        # print(mutatedOffspringEdge)
        mutatedOffspringNode = decode_communities(locus, mutatedOffspringEdge)
        # print(mutatedOffspringNode)
        base, middle, top = updatePopulation(mutatedOffspringNode, base)
        if (topScore < top.pocketFitness):
            topScore = top.pocketFitness
            count = 0
        else:
            count += 1
        print('*************************************')
        print(top.pocket, '|| Fitness:', top.pocketFitness)


M_Link()

# %%
drawLayout(labledEdge)
