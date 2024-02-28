# simplified: one generator and one load at each node

import pandas as pd
import numpy as np
import os
import datetime
import random

import gurobipy as gp
from gurobipy import GRB, max_
from gurobipy import quicksum
from gurobipy import LinExpr
import scipy.sparse as sp

def read_data_from_files(case):

    gen_file = rf"data/gen_data_case{case}.csv"
    bus_file = rf"data/bus_data_case{case}.csv"
    branch_file = rf"data/branch_data_case{case}.csv"
    gen_cost_file = rf"data/gencost_data_case{case}.csv"

    gen_data = pd.read_csv(gen_file)
    bus_data = pd.read_csv(bus_file)
    branch_data = pd.read_csv(branch_file)
    gencost_data = pd.read_csv(gen_cost_file)

    gen_data_case = gen_data.loc[gen_data['case'] == case]
    bus_data_case = bus_data.loc[bus_data['case'] == case]
    branch_data_case = branch_data.loc[branch_data['case'] == case]
    gencost_data_case = gencost_data.loc[gencost_data['case'] == case]

    nodes = int(bus_data_case['bus'].max())
    P_load = (-bus_data_case.Pd).tolist()

    gen_bus = gen_data_case.bus.tolist()
    gen_upper_bounds = gen_data_case.Pmax.tolist()
    gen_costs = gencost_data_case.c1.tolist()

    gens = len(gen_bus)

    gen_bus_dict = {}

    for (g,n) in enumerate(gen_bus):
        if n in gen_bus_dict:
              gen_bus_dict[n].append(g)
        else:
            gen_bus_dict[n] = [g]
    
    for i in range(1,nodes):
        if i in gen_bus_dict:
            pass
        else:
            gen_bus_dict[i] = []

    branch_limits = np.zeros((nodes,nodes))
    B = np.zeros((nodes,nodes))

    for i,row in branch_data_case.iterrows():
        i = int(row.fbus)-1
        j = int(row.tbus)-1
        branch_limits[i,j] = row.rateA
        branch_limits[j,i] = row.rateA
        B[i,j] = 1 / row.x
        B[j,i] = 1 / row.x

        if branch_limits[i,j] == 0:
            B[i,j] = 0

    

    return (nodes,
            gens,
            P_load,
            gen_bus_dict,
            gen_upper_bounds,
            gen_costs,
            branch_limits,
            B,
            branch_data_case)

def write_results(m, P_load, branch_data_case, r_g=None, w_bar=None, carbon_model=False):

    node_output = pd.DataFrame(index=pd.RangeIndex(1,nodes+1,name='Bus'),
                               columns=['Load','Generation','Carbon Intensity','Voltage']) 
    gen_output = pd.DataFrame(index=pd.RangeIndex(1,gens+1, name='Generator'),
                              columns=['Generation','Cost','Emissions'])
    branch_output = pd.DataFrame(
                        index = pd.MultiIndex(levels=[[],[]],
                                            codes=[[],[]],
                                            names=['from', 'to']),
                        columns = ['Flow', 'Pn', 'Pb'])        
    

    node_output.loc[:,'Load'] = P_load

    if carbon_model:
        node_output.loc[:,'Carbon Intensity'] = w_bar

    for i in range(nodes):
        node_output.loc[i+1,'Generation'] = quicksum(m.getVarByName(f"Gen[{k}]").X for k in gen_bus_dict[i+1]).getValue()
        node_output.loc[i+1,'Voltage'] = m.getVarByName(f"Voltage[{i}]").X

    for i in range(gens):
        gen_output.loc[i+1, 'Generation'] = m.getVarByName(f"Gen[{i}]").X
        gen_output.loc[i+1, 'Cost'] = gen_costs[i] * m.getVarByName(f"Gen[{i}]").X

        if carbon_model:
            gen_output.loc[i+1, 'Emissions'] = r_g[i] * m.getVarByName(f"Gen[{i}]").X

    for pair in branch_data_case.loc[:,['fbus','tbus']].values:
        i = pair[0]-1
        j = pair[1]-1
        if carbon_model:
            branch_output.loc[tuple(pair),:] = \
                [m.getVarByName(f"Flow[{i},{j}]").X,
                m.getVarByName(f"Pn[{i},{j}]").X,
                m.getVarByName(f"Pb[{i},{j}]").X]
        else:
            branch_output.loc[tuple(pair),'Flow'] = m.getVarByName(f"Flow[{i},{j}]").X
            
    time_str = datetime.datetime.now().strftime(
        "%d_%m_%Y_%H_%M_%S"
    )

    dir_name = "_output_" + time_str
    os.mkdir(dir_name)

    node_output.to_csv(os.path.join(dir_name, "node_output.csv"))
    gen_output.to_csv(os.path.join(dir_name, "gen_output.csv"))
    branch_output.to_csv(os.path.join(dir_name, "branch_output.csv"))

def create_opf_model(nodes,gens,P_load,gen_bus_dict,gen_costs,branch_limits,B,gen_upper_bounds,r_g=None,w_bar=None):
    
    m = gp.Model()

    # add variables
    P_gen = m.addMVar((gens),lb=np.zeros(gens),ub = gen_upper_bounds, vtype= GRB.CONTINUOUS,name="Gen")
    voltage_angle = m.addVars(nodes,lb=-1000*np.ones(nodes),vtype= GRB.CONTINUOUS, name="Voltage")
    Flow = m.addVars(nodes,nodes,lb=-branch_limits, ub = branch_limits, vtype = GRB.CONTINUOUS, name="Flow")
    if (not r_g is None) and (not w_bar is None):
        P_n = m.addMVar((nodes,nodes), lb=np.zeros((nodes,nodes)),vtype= GRB.CONTINUOUS,name="Pn")
        P_b = m.addMVar((nodes,nodes), lb=np.zeros((nodes,nodes)),vtype= GRB.CONTINUOUS,name="Pb")
        Rg = sp.diags(r_g)
    # add constraints

    # Net flow = load + gen at each node

    m.addConstrs(quicksum(-Flow[j,i] for j in range(nodes)) 
                == P_load[i] 
                + quicksum(P_gen[k] for k in gen_bus_dict[i+1])
                for i in range(nodes))

    # branch limits 

    m.addConstrs(Flow[i,j] <= branch_limits[i,j] 
                for i in range(nodes)
                for j in range(nodes))
    
    # voltage angle

    m.addConstrs(Flow[i,j] == B[i,j] * (voltage_angle[i] - voltage_angle[j])
                 for i  in range(nodes)
                 for j in range(nodes))

    m.addConstr(voltage_angle[0] == 0)
    
    if (not r_g is None) and (not w_bar is None) :
        for j in range(nodes):
            for i in range(nodes):
                m.addGenConstrMax(P_b[j,i],[Flow[i,j],0])
        
        m.addConstrs((quicksum(P_b[i,j] for j in range(nodes)) \
                      + quicksum(P_gen[k] for k in gen_bus_dict[i+1])
                     == P_n[i,i] for i in range(nodes)))

        m.addConstrs(0 == P_n[i,j] for i in range(nodes) for j in range(i+1,nodes))
        m.addConstrs(0 == P_n[i,j] for i in range(nodes) for j in range(0,i))

        # TODO: fix shape problem
        # only 12 gens, 89 nodes
        # need new decision variable?
        m.addConstr(Rg@P_gen <= (P_n-P_b) @ w_bar)

    # add objective

    m.setObjective(quicksum(gen_costs[i] * P_gen[i] for i in range(gens)))

    return m

def generate_time_series_loads(load_profile, uncertainty, P_load):
    avg_load = sum(load_profile)/len(load_profile)
    normalized_load_profile = np.divide(load_profile, avg_load)
    load_series = np.zeros((len(load_profile),len(P_load)))
    for i in range(len(normalized_load_profile)):
        scalar = normalized_load_profile[i]
        for j in range(len(P_load)):
            load = P_load[j]
            fudge = random.random()*uncertainty*2-uncertainty
            load_series[i][j] = (load*scalar)*(1+fudge)
    return load_series

def run_time_series(load_series,nodes,gens,gen_bus_dict,gen_costs,branch_limits,B,gen_upper_bounds,r_g=None,w_bar=None):
    for i in range(len(load_series)):
        P_load = load_series[i]
        m = create_opf_model(nodes,gens,P_load,gen_bus_dict,gen_costs,branch_limits,B,gen_upper_bounds,r_g,w_bar)
        m.optimize()
        print("time:", i)
        print(m.getVars())

case = 89

(nodes,
gens,
P_load,
gen_bus_dict,
gen_upper_bounds,
gen_costs,
branch_limits,
B,
branch_data_case) = read_data_from_files(case)

# nodes = 3
# P_load = np.array([-1,0,-3])

# branch_limits = np.array([[0,10,0],[10,0,10],[0,10,0]])

generator_carbon = np.ones(gens)

carbon_upper_bounds = np.ones(nodes)

# print(nodes,gens,P_load,gen_bus_dict,gen_costs,branch_limits,
#                     B,
#                     gen_upper_bounds)

# branch_reactance = np.ones((nodes,nodes))
# B = 1/branch_reactance
# for i in range(len(branch_limits)):
#     for j in range(len(branch_limits[i])):
#         if branch_limits[i,j] == 0:
#             B[i,j] = 0

# gen_upper_bounds = np.array([2,0,3]) # generator capacity
# gen_costs = np.array([14,15,12]) # gen cost per ouput

# load_profile = [1,2,3,1,1]
# load_series = generate_time_series_loads(load_profile,0.1,P_load)
# run_time_series(load_series,nodes,
#                     gens,
#                     gen_bus_dict,
#                     gen_costs,
#                     branch_limits,
#                     B,
#                     gen_upper_bounds)

# m = create_opf_model(nodes,
#                     gens,
#                     P_load,
#                     gen_bus_dict,
#                     gen_costs,
#                     branch_limits,
#                     B,
#                     gen_upper_bounds)

# run with carbon

m = create_opf_model(nodes,
                    gens,
                    P_load,
                    gen_bus_dict,
                    gen_costs,
                    branch_limits,
                    B,
                    gen_upper_bounds,
                    generator_carbon,
                    carbon_upper_bounds)

m.optimize()
# # m.computeIIS()
# # m.write("test.ilp")
print(m.getVars())

write_results(m, P_load, branch_data_case)




# dummy parameters

# nodes = 3
# P_load = np.array([-1,0,-3])

# branch_limits = np.array([[0,10,0],[10,0,10],[0,10,0]])

# generator_carbon = np.array([1,0,2])

# carbon_upper_bounds = np.array([1.25,2,2])

# branch_reactance = np.array([[1,1,1],[1,1,1],[1,1,1]])
# B = 1/branch_reactance
# for i in range(len(branch_limits)):
#     for j in range(len(branch_limits[i])):
#         if branch_limits[i,j] == 0:
#             B[i,j] = 0

# gen_upper_bounds = np.array([2,0,2]) # generator capacity
# gen_costs = np.array([10,15,12]) # gen cost per ouput