# simplified: one generator and one load at each node

import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import GRB, max_
from gurobipy import quicksum
from gurobipy import LinExpr

def read_data_from_files(gen_file,bus_file,branch_file,gen_cost_file,case):


    gen_data = pd.read_csv(gen_file)
    bus_data = pd.read_csv(bus_file)
    branch_data = pd.read_csv(branch_file)
    gencost_data = pd.read_csv(gen_cost_file)

    gen_data_case = gen_data.loc[gen_data['case'] == case]
    bus_data_case = bus_data.loc[bus_data['case'] == case]
    branch_data_case = branch_data.loc[branch_data['case'] == case]
    gencost_data_case = gencost_data.loc[gencost_data['case'] == case]

    nodes = bus_data_case['bus'].max()
    P_load = (-1 * bus_data_case.Pd).tolist()

    gen_bus = gen_data_case.bus.tolist()
    gen_upper_bounds = gen_data_case.Pmax.tolist()
    gen_costs = gencost_data_case.c1.tolist()

    branch_limits = np.zeros((nodes,nodes))
    B = np.zeros((nodes,nodes))

    for i,row in branch_data_case.iterrows():
        i = int(row.fbus)-1
        j = int(row.tbus)-1

        branch_limits[i,j] = row.rateA
        B[i,j] = 1 / row.x

        branch_limits[j,i] = row.rateA
        B[j,i] = 1 / row.x

        if branch_limits[i,j] == 0:
                B[i,j] = 0
                B[j,i] = 0

    return (nodes,
            P_load,
            gen_bus,
            gen_upper_bounds,
            gen_costs,
            branch_limits,
            B)


def create_opf_model(nodes,gen_costs,branch_limits,B,gen_upper_bounds,generator_carbon,carbon_upper_bounds):
    
    m = gp.Model()

    # add variables
    P_gen = m.addVars(nodes,lb=np.zeros(nodes),ub = gen_upper_bounds, vtype= GRB.CONTINUOUS)
    voltage_angle = m.addVars(nodes,lb=-1000*np.ones(nodes),vtype= GRB.CONTINUOUS)
    Flow = m.addVars(nodes,nodes,lb=-1000*np.ones((nodes,nodes)), ub = branch_limits, vtype = GRB.CONTINUOUS)
    #P_in = m.addMVar((nodes,nodes), lb=np.zeros((nodes,nodes)),vtype= GRB.CONTINUOUS)
    #P_b = m.addMVar((nodes,nodes), lb=np.zeros((nodes,nodes)),vtype= GRB.CONTINUOUS)

    # add constraints

    # Net flow = load + gen at each node

    m.addConstrs(quicksum(-Flow[j,i] for j in range(nodes)) 
                == P_load[i] + P_gen[i] for i in range(nodes))

    # branch limits 

    m.addConstrs(Flow[i,j] <= branch_limits[i,j] 
                for i in range(nodes)
                for j in range(nodes))
    
    # voltage angle

    m.addConstrs(Flow[i,j] == B[i,j] * (voltage_angle[i] - voltage_angle[j])
                 for i  in range(nodes)
                 for j in range(nodes))

    m.addConstr(voltage_angle[0] == 0)

    # for j in range(nodes):
    #     for i in range(nodes):
    #         m.addGenConstrMax(P_b[j,i],[Flow[j,i],0])
    
    # m.addConstrs((quicksum(P_b[j,i] for j in range(nodes)) 
    #             == P_in[i,i] for i in range(nodes)))

    # m.addConstrs(0 == P_in[i,j] for i in range(nodes) for j in range(i+1,nodes))
    # m.addConstrs(0 == P_in[i,j] for i in range(nodes) for j in range(0,i))

    # # m.addConstrs(np.dot((P_in-P_b,carbon_upper_bounds)) >= generator_carbon)
    # m.addConstr((P_in-P_b) @ carbon_upper_bounds >= generator_carbon)

    # add objective

    m.setObjective(quicksum(gen_costs[i] * P_gen[i] for i in range(nodes)))

    return m


gen_file = r"data/gen_data_case3.csv"
bus_file = r"data/bus_data_case3.csv"
branch_file = r"data/branch_data_case3.csv"
gen_cost_file = r"data/gencost_data_case3.csv"
case = 3

(nodes,
P_load,
gen_bus,
gen_upper_bounds,
gen_costs,
branch_limits,
B) = read_data_from_files(gen_file,bus_file,branch_file,gen_cost_file,case)

generator_carbon = np.array([1,0,2])
carbon_upper_bounds = np.array([1.25,2,2])


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

m = create_opf_model(nodes,
                     gen_costs,
                     branch_limits,
                     B,
                     gen_upper_bounds,
                     generator_carbon,
                     carbon_upper_bounds)

m.optimize()
print(m.getVars())





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