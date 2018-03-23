"""
Created on Fri Mar 24 15:36:30 2017

@author: Priyanka Mocherla
@version: 2.7

This code contains functions to model a random network (ER model) and a preferential attachment network (BA model). The code also produces a data collapse of the simulation data and uses statistical methods to extract features of the network such as the cut-off degree value.

Also includes the code required to complete the project tasks.

Functions in this module:
    - initial_graph
    - BA_model_visual
    - BA_model_degree
    - data_prob
    - random_model
    - random_walk
    - join_BA
    - join_random
    - random_theory_func
    - preferential_theory_func
    - average_height
    - stand_dev
    - cut_off_theoreticalpref
    - cut_off_theoreticalrand

"""

#------------------------------- Functions -----------------------------#
import networkx as nx
import numpy as np
import random as r
import matplotlib.pyplot as plt
from log_bin import *
import collections
from scipy.stats import chisquare



def initial_graph(start_nodes):
    """ 
       returns a networkx inital sparse graph
        
        Args: 
        start_nodes: integer, number of nodes starting in initial graph
    """
    G = nx.MultiGraph()
    G.add_node(1)
    for i in np.arange(2,start_nodes+1):
        G.add_node(i)
        trial = [(1,1)]
        while trial[0][0] == trial[0][1]:
            trial = [(i,r.randint(1,i))]
        else:
            G.add_edges_from(trial)

    return G


def BA_model_visual(start_nodes,time,m):
    """ 
       returns a network connected by preferential attachment according to the BA model. For visual use only
        
        Args: 
        start_nodes: integer, number of nodes starting in initial graph
        time: integer, time the network is run for, aka total number of nodes
        m: integer, number of edges per new node
    """
    G = initial_graph(start_nodes)
    m_data = []
    new_edge = list(sum(G.edges(),()))
    for i in np.arange(1,time-start_nodes+1):
        G.add_node(start_nodes+i)
        while len(m_data) < m:
            m_data.append(r.choice(new_edge))
        else:
            for j in m_data:
                G.add_edge(j,start_nodes+i)
                m_data = []
            new_edge = list(sum(G.edges(),()))
    return G
    
def BA_model_degree(start_nodes,time,m):
    """ 
        returns a dictionary of the degree of each node in the network for the BA model
        
        Args: 
        start_nodes: integer, number of nodes starting in initial graph
        time: integer, time the network is run for, aka total number of nodes
        m: integer, number of edges per new node
    """
    G = initial_graph(start_nodes)
    degrees = list(sum(G.edges(),()))
    for i in np.arange(1,time-start_nodes+1):
        for j in range(m):
            degrees.append(r.choice(degrees))
            degrees.append(start_nodes+i)
    node_degrees = collections.Counter(degrees)
    return node_degrees

def data_prob(data):
    """ 
       returns a dictionary of data keys and corresponding frequency probability.
        
        Args: 
        data: array eg. np.array([x1,x2,x3...])
    """
    data = np.array(data)
    dataprob = {}
    for datum in data:
        if datum not in dataprob:
            dataprob[datum] = 0
        dataprob[datum] = dataprob[datum] + 1./data.size
    return dataprob
    
def random_model(start_nodes,time,m):
    """ 
        returns a dictionary of the degree of each node in the network for the random model
        
        Args: 
        start_nodes: integer, number of nodes starting in initial graph
        time: integer, time the network is run for, aka total number of nodes
        m: integer, number of edges per new node
    """
    G = initial_graph(start_nodes)
    degrees = list(sum(G.edges(),()))
    for i in np.arange(1,time-start_nodes+1):
        for j in range(m):
            degrees.append(r.randint(1,start_nodes+i-1))
            degrees.append(start_nodes+i)
    node_degrees = collections.Counter(degrees)
    return node_degrees

def random_walk(start_nodes,time,m,steps):
    """ 
        returns a networkx graph of the network for the random walk attachment model
        
        Args: 
        start_nodes: integer, number of nodes starting in initial graph
        time: integer, time the network is run for, aka total number of nodes
        m: integer, number of edges per new node
        steps: integer, number of steps per random walk
    """
    G = initial_graph(start_nodes)
    for i in np.arange(1,time-start_nodes+1):
        G.add_node(start_nodes+i)
        for j in range(m):
            start = r.randint(1,start_nodes+i-1)
            if steps == 0:
                G.add_edge(start_nodes+i,start)
            else:
                for k in range(steps):
                    poss = G.neighbors(start)
                    start = r.choice(poss)
                G.add_edge(start_nodes+i,start)
    return G
    
def join_BA(time,runs,start_nodes,m):
    """ 
        returns a concatenated list of the degree distribution for multiple runs and the cut-off degree of each run for BA
        
        Args: 
        time: integer, time the network is run for, aka total number of nodes
        runs: integer, number of runs per network
        start_nodes: integer, number of nodes starting in initial graph
        m: integer, number of edges per new node
    """   
    cut_off = np.zeros(runs)
    start = np.array(BA_model_degree(start_nodes,time,m).values())
    for i in np.arange(runs):
        print i
        dist = BA_model_degree(start_nodes,time,m)
        add = np.array(dist.values())
        start = np.concatenate((start,add))
        cut_off[i] = max(dist.values())
        
    return start,cut_off
    
def join_random(time,runs,start_nodes,m):
    """ 
        returns a concatenated list of the degree distribution for multiple runs and the cut-off degree of each run for random
        
        Args: 
        time: integer, time the network is run for, aka total number of nodes
        runs: integer, number of runs per network
        start_nodes: integer, number of nodes starting in initial graph
        m: integer, number of edges per new node
    """   
    cut_off = np.zeros(runs)
    start = np.array(random_model(start_nodes,time,m).values())
    for i in np.arange(runs):
        print i
        dist = random_model(start_nodes,time,m)
        add = np.array(dist.values())
        start = np.concatenate((start,add))
        cut_off[i] = max(dist.values())
        
    return start,cut_off
    
def random_theory_func(m,k):
    """ 
        returns the theoretical function for the random attachment degree distribution
        
        Args: 
        k: array eg. np.array([x1,x2,x3...]), the degree values to evaluate
        m: integer, number of edges per new node
    """   
    k = np.array(k)
    values = np.zeros(k.size)
    for i in range(k.size):
        if k[i] < m:
            continue
        else:
            values[i] = m**(k[i]-m)/(1+m)**(1+k[i]-m)
    return values

def preferential_theory_func(m,k):
    """ 
        returns the theoretical function for the preferential attachment degree distribution
        
        Args: 
        k: array eg. np.array([x1,x2,x3...]), the degree values to evaluate
        m: integer, number of edges per new node
    """   
    k = np.array(k)
    values = np.zeros(k.size)
    for i in range(k.size):
        if k[i] < m:
            continue
        else:
            values[i] = 2*m*(m+1)/((k[i]+2)*(k[i]+1)*k[i])
    return values
    
def average_height(heights):
    """ 
       returns the average value of an input array.
        
        Args: 
        heights: array eg. np.array([x1,x2,x3...])
    """
    return np.sum(heights)/heights.size
    
def stand_dev(heights):
    """ 
       returns the standard deviation of an input array.
        
        Args: 
        heights: array eg. np.array([x1,x2,x3...])
    """
    return np.sqrt(np.sum(heights*heights)/heights.size-(np.sum(heights)/heights.size)*(np.sum(heights)/heights.size))
    
def cut_off_theoreticalpref(N,m):
    """ 
       returns theoretical cut-off value for preferential attachment
        
        Args: 
        N: integer, system size
        m: integer, edges per node
    """
    k1 = 0.5*(-1+np.sqrt(1+4*N*m*(m+1)))
    return k1

def cut_off_theoreticalrand(N,m):
    """ 
       returns theoretical cut-off value for random attachment
        
        Args: 
        N: integer, system size
        m: integer, edges per node
    """
    k1 = m - (np.log(N)/(np.log(m)-np.log(m+1)))
    return k1


#----------------------------------- Project Code -----------------------------------#
"""
#TASK 2B: DATACOLLAPSE Random

time = [100,1000,10000,100000,1000000]
m = 4
runs = [10000,1000,100,50,25]

#DO NOT TOUCH - PLOTTING
plt.figure()

colours = ["red","blue","green", "orange", "purple","black"]# enter same number of colours as m size

for i in range(len(runs)):
    data = join_random(time[i],runs[i],m+1,m)[0]
    distribution = data_prob(data)
    k1 = cut_off_theoreticalrand(time[i],m)
    x,y = distribution.keys(),distribution.values()
    b,c = log_bin(np.array(data), 1.,1.5, 1.5, debug_mode=True, drop_zeros=False)
    plt.loglog(b/k1,c/random_theory_func(m,b),'-', color = colours[i], label = 'N = '+ str(time[i]))
    data = 0
    
    
#plt.title("Random attachment Data collapse: m = " +str(m) )
plt.xlabel("$k/k_1$", fontsize = 22)
plt.ylabel("$p(k)/p_\infty(k)$", fontsize = 22)
plt.legend(loc = 3)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.show()    
        
"""  
###############################################################################             
"""
#Print the number of duplicates
m = 2
time = 50000

counter = 0

edges = []
dups = []

degrees = BA_model_degree(m+1,time,m)[1]
            
for i in np.arange(0,(len(degrees)/2.),1):
    edges.append((degrees[2*int(i)],degrees[2*int(i)+1]))

    
for j in edges:
    if j in dups:
        counter = counter+1
    else:
        dups.append(j)

print counter
"""
"""
#Plotting the BA Model graph
G = BA_model_visual(start,time,3)
r.seed(2)

color_map = []
for node in G:
    if node <=start:
        color_map.append('#FF8200')#starting
    else: color_map.append('#22BDE6')#blue
    
plt.figure()
nx.draw_shell(G,node_color = color_map,with_labels = True)
"""
###############################################################################
"""
#TASK 1C: PREFERENTIAL ATTACHMENT - FIXED N VARYING m

#Initialise parameters
time = 100000
m = [1,2,4,16,64]
runs = [1,1,1,1,1]

#DO NOT TOUCH - PLOTTING
plt.figure()

colours = ["red","blue","green", "orange", "purple","black"]# enter same number of colours as m size

for i in range(len(m)):
    data = join_BA(time, runs[i], m[i]+1,m[i])[0]
    distribution = data_prob(data)
    
    x,y = distribution.keys(),distribution.values()
    plt.loglog(x,y,'o', color = colours[i], label = 'm = '+ str(m[i]))
    
    #b,c = log_bin(np.array(data), 1.,1.5, 1.3, debug_mode=True, drop_zeros=False)
    #plt.loglog(b,c, 'o', color = colours[i])
    k = np.arange(1,x[-1],0.5)
    plt.loglog(k, preferential_theory_func(m[i],k),'--', color = colours[i])
    data = 0

#plt.title("Preferential attachment: N = " +str(time) )
plt.xlabel("$k$", fontsize = 22)
plt.ylabel("$p(k)$", fontsize = 22)
plt.legend(loc = 1)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.show()
"""
###################################################################################
"""
#TASK 1C: FITTING DATA

#Initialise parameters
time = 1000000
m = [1,2,4,16,64]

#DO NOT TOUCH - PLOTTING
start = [3,5,7,9,11]

colours = ["red","blue","green", "orange", "purple","black"]# enter same number of colours as m size

for i in range(len(start)):
    plt.figure()
    data = BA_model_degree(m[i]+1,time,m[i])
    raw = data_prob(np.array(data.values()))
    
    x,y = raw.keys(),raw.values()
    b, c = log_bin(np.array(data.values()), 1.,1.5, 1.5, debug_mode=True, drop_zeros=False)
    p,q = np.polyfit(np.log(b)[start[i]:-start[i]],np.log(c)[start[i]:-start[i]],1)
    
    print '#####################  Gradient: ' + str(p)
    k = np.arange(1,b[-1],0.5)
    plt.loglog(x,y,'o',label = 'Raw data')
    plt.loglog(k,preferential_theory_func(m[i],k),'--', label = 'Theoretical')
    plt.loglog(b,np.exp(p*np.log(b)+q),'-', label = 'Scale fitted')
    plt.title("Fitting data: N = " +str(time) + ' , m = ' +str(m[i]))
    print chisquare(np.exp(p*np.log(b)+q)[start[i]:-start[i]]*time,preferential_theory_func(m[i],b)[start[i]:-start[i]]*time)
    plt.legend(loc = 1)

plt.xlabel("$k$", fontsize = 22)
plt.ylabel("$p(k)$", fontsize = 22)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.show()
"""
################################################################################
"""
#TASK 1D: FIXED m VARYING N

time = [100,1000,10000,100000,1000000]
m = 4
runs = [10000,1000,100,50,25]

#DO NOT TOUCH - PLOTTING
plt.figure()

colours = ["red","blue","green", "orange", "purple","black"]# enter same number of colours as m size

for i in range(len(time)):
    data = join_BA(time[i],runs[i],m+1,m)[0]
    distribution = data_prob(data)
    
    x,y = distribution.keys(),distribution.values()
    b,c = log_bin(np.array(data), 1.,1.5, 1.2, debug_mode=True, drop_zeros=False)
    plt.loglog(b,c,'o', color = colours[i], label = 'N = '+ str(time[i]))
    data = 0
    
    
k = np.arange(1,x[-1],0.5)
plt.loglog(k, preferential_theory_func(m,k),'--', color = "black")
    
#plt.title("Preferential attachment: m = " +str(m) )
plt.xlabel("$k$", fontsize = 22)
plt.ylabel("$p(k)$", fontsize = 22)
plt.legend(loc = 1)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.show()
"""
################################################################################
"""
#TASK 1D: CUT-OFF

#Initialise parameters
time = [100,1000,10000,100000,1000000]
m = 4 
runs = [10000,1000,100,50,25]

#DO NOT TOUCH
k1 = np.zeros(len(time))
error = np.zeros(len(time))

for i in range(len(time)):
    cutoff = join_BA(time[i], runs[i], m+1,m)[1]
    k1[i] = average_height(cutoff)
    error[i] = stand_dev(cutoff)

a,b = np.polyfit(np.log(time),np.log(k1),1)
print a    
test = np.polyfit(np.log(time),np.log(k1),1,full = True)[1]
print test
plt.plot(time,np.exp(np.log(time)*a+b),'r-')
plt.errorbar(time,k1,yerr = error,fmt='ro')

N = np.arange(50,time[-1]*2)
plt.plot(N, cut_off_theoreticalpref(N,m),'b--')
#plt.title("Cut-off degree with system-size: Preferential" )
plt.xlabel("$N$", fontsize = 22)
plt.ylabel("$k_1$", fontsize = 22)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.xlim([50,time[-1]*2])
plt.yscale('log')
plt.xscale('log')
plt.show()
"""
################################################################################
"""
#TASK 1D: DATACOLLAPSE

time = [100,1000,10000,100000,1000000]
m = 4
runs = [10000,1000,100,50,25]

#DO NOT TOUCH - PLOTTING
plt.figure()

colours = ["red","blue","green", "orange", "purple","black"]# enter same number of colours as m size

for i in range(len(runs)):
    data = join_BA(time[i],runs[i],m+1,m)[0]
    distribution = data_prob(data)
    k1 = cut_off_theoreticalpref(time[i],m)
    x,y = distribution.keys(),distribution.values()
    b,c = log_bin(np.array(data), 1.,1.5, 1.5, debug_mode=True, drop_zeros=False)
    plt.loglog(b/k1,c/preferential_theory_func(m,b),'-', color = colours[i], label = 'N = '+ str(time[i]))
    data = 0
    
    
#plt.title("Preferential attachment Data collapse: m = " +str(m) )
plt.xlabel("$k/k_1$", fontsize = 22)
plt.ylabel("$p(k)/p_\infty(k)$", fontsize = 22)
plt.legend(loc = 3)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.show()
"""
################################################################################
"""
#TASK 2B: RANDOM ATTACHMENT - FIXED N VARYING m

#Initialise parameters
time = 1000000
m = [1,2,4,8,16,32]

#DO NOT TOUCH - PLOTTING
plt.figure()

colours = ["red","blue","green", "orange", "purple","black"]# enter same number of colours as m size

for i in range(len(m)):
    data = random_model(m[i]+1,time,m[i])
    distribution = data_prob(data.values())
    
    x,y = distribution.keys(),distribution.values()
    b,c = log_bin(np.array(data.values()), 1.,1.5, 1.2, debug_mode=True, drop_zeros=False)
    plt.loglog(b,c,'o', color = colours[i], label = 'm = '+ str(m[i]))
    
    k = np.arange(1,x[-1]*10,0.1)
    plt.loglog(k, random_theory_func(m[i],k),'--', color = colours[i])

#plt.title("Random attachment: N = " +str(time) )
plt.xlabel("$k$", fontsize = 22)
plt.ylabel("$p(k)$", fontsize = 22)
plt.legend(loc = 3)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.show()
"""
###############################################################################
"""
#TASK 2B: FIXED m VARYING N

time = [100,1000,10000,100000,1000000]
m = 4
runs = [10000,1000,200,100,25]

#DO NOT TOUCH - PLOTTING
plt.figure()

colours = ["red","blue","green", "orange", "purple","black"]# enter same number of colours as m size

for i in range(len(time)):
    data = join_random(time[i],runs[i],m+1,m)[0]
    distribution = data_prob(data)
    
    x,y = distribution.keys(),distribution.values()
    b,c = log_bin(np.array(data), 1.,1.5, 1.05, debug_mode=True, drop_zeros=False)
    plt.loglog(b,c,'o', color = colours[i], label = 'N = '+ str(time[i]))
    #plt.loglog(x,y,'o',color = colours[i], label = 'N = '+ str(time[i]))
    data = 0
    
    
k = np.arange(1,x[-1],0.5)
plt.loglog(k, random_theory_func(m,k),'-', color = "black")
    
#plt.title("Random attachment: m = " +str(m) )
plt.xlabel("$k$", fontsize = 22)
plt.ylabel("$p(k)$", fontsize = 22)
plt.legend(loc = 1)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.show()
"""
##############################################################################3
"""
#TASK 2B: CUT-OFF RANDOM
#Initialise parameters
time = [100,1000,10000,100000,1000000]
m = 4
runs = [10000,1000,200,100,25]

#DO NOT TOUCH
k1 = np.zeros(len(time))
error = np.zeros(len(time))

for i in range(len(time)):
    cutoff = join_random(time[i], runs[i], m+1,m)[1]
    k1[i] = average_height(cutoff)
    error[i] = stand_dev(cutoff)

N = np.arange(50,time[-1]*2)
plt.plot(N, cut_off_theoreticalrand(N,m),'r--')

plt.errorbar(time,k1,yerr = error,fmt='ro')
#plt.title("Cut-off degree with system-size: Random" )
plt.xlabel("$N$", fontsize = 22)
plt.ylabel("$k_1$", fontsize = 22)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.yscale('log')
plt.xscale('log')
plt.show()
"""
################################################################################
"""
#TASK 2B: FITTING DATA

#Initialise parameters
time = 1000000
m = [1,2,4,8,16]

#DO NOT TOUCH - PLOTTING
start = [0,1,2,3,4]

for i in range(len(start)):
    plt.figure()
    data = random_model(m[i]+1,time,m[i])
    raw = data_prob(np.array(data.values()))
    
    x,y = raw.keys(),raw.values()
    b, c = log_bin(np.array(data.values()), 1.,1.5, 1.5, debug_mode=True, drop_zeros=False)

    k = np.arange(1,b[-1],0.5)
    plt.loglog(x,y,'o',label = 'Raw data')
    plt.loglog(k,random_theory_func(m[i],k),'--', label = 'Theoretical')
    plt.loglog(b,c,'-', label = 'Log binned')
    plt.title("Fitting data Random: N = " +str(time) + ' , m = ' +str(m[i]))
    print chisquare(c[start[i]::],random_theory_func(m[i],b)[start[i]::])
    plt.legend(loc = 1)

plt.xlabel("$k$", fontsize = 22)
plt.ylabel("$p(k)$", fontsize = 22)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.show()
"""
###############################################################################
"""
#TASK 3: RANDOM WALK
#Initialise N and m
time = 200000
m = 4 

#DO NOT TOUCH
diam = int(round(np.log(time)/np.log(np.log(time))))

#Walk Lengths
zero = random_walk(m+1,time,m,0)
one = random_walk(m+1,time,m,1)
diameter = random_walk(m+1,time,m,diam)

                ################ Plotting ###############
#L = 0
plt.figure()
zero_dist = (zero.degree()).values()
zero_dist1 = data_prob(np.array(zero_dist))
a,b = np.array(zero_dist1.keys()),zero_dist1.values()
k0 = np.arange(1,a[-1],0.5)

plt.loglog(a, b, 'ro',label = 'L = 0')
plt.loglog(k0, random_theory_func(m,k0),'b--')
#plt.title("Random walk: L = 0")
plt.xlabel("$k$", fontsize = 22)
plt.ylabel("$p(k)$", fontsize = 22)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.legend(loc=1)

#L = 1
plt.figure()
one_dist = (one.degree()).values()
one_dist1 = data_prob(np.array(one_dist))
c,d = np.array(one_dist1.keys()),one_dist1.values()#plotting data
p,q = log_bin(one_dist, 1.,1.5, 1.5, debug_mode=True, drop_zeros=False) #logbinning

h,i = np.polyfit(np.log(p[5:15]),np.log(q[5:15]),1)#fitting
print h

error1 = np.polyfit(np.log(p[5:15]),np.log(q[5:15]),1,full=True)[1]
print error1

plt.loglog(p, q, 'b--', label = 'L = 1')
plt.loglog(c, d, 'ro')
#plt.title("Random walk: L = 1")
plt.xlabel("$k$", fontsize = 22)
plt.ylabel("$p(k)$", fontsize = 22)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.legend(loc=1)


#L = diameter
plt.figure()
diameter_dist = (diameter.degree()).values()
diameter_dist1 = data_prob(np.array(diameter_dist))
e,f = np.array(diameter_dist1.keys()),diameter_dist1.values()#plotting data
r,s = log_bin(diameter_dist, 1.,1.5, 1.5, debug_mode=True, drop_zeros=False) #logbinning

j,k = np.polyfit(np.log(r[5:15]),np.log(s[5:15]),1)#fitting
print j
error2 = np.polyfit(np.log(r[5:15]),np.log(s[5:15]),1,full=True)[1]
print error2


plt.loglog(r, s, 'b--',label = 'L = '+ str(diam))
plt.loglog(e, f, 'ro')
#plt.title("Random walk: L = diameter")
plt.xlabel("$k$", fontsize = 22)
plt.ylabel("$p(k)$", fontsize = 22)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.legend(loc=1)
plt.show()
"""
