# simplified: one generator and one load at each node

import numpy as np
import gurobipy as gp
from gurobipy import GRB
from gurobipy import quicksum
from gurobipy import LinExpr


# dummy parameters

nodes = 3
P_load = np.array([-1,0,-3])

branch_limits = np.array([[0,10,0],[10,0,10],[0,10,0]])

branch_reactance = np.array([[1,1,1],[1,1,1],[1,1,1]])
B = 1/branch_reactance
for i in range(len(branch_limits)):
    for j in range(len(branch_limits[i])):
        if branch_limits[i,j] == 0:
            B[i,j] = 0

gen_upper_bounds = np.array([2,0,2]) # generator capacity
gen_costs = np.array([10,15,12]) # gen cost per ouput

def create_opf_model(nodes,gen_costs,branch_limits,B,gen_upper_bounds):
    
    m = gp.Model()

    # add variables

    P_gen = m.addVars(nodes,lb=[0,0,0],ub = gen_upper_bounds, vtype= GRB.CONTINUOUS,name=["G1", "G2", "G3"])
    voltage_angle = m.addVars(nodes,lb=[-1000,-1000,-1000],vtype= GRB.CONTINUOUS,name=["v1", "v2", "v3"])
    flow_names = [["F"+str(j+1)+str(i+1) for i in range(3)] for j in range(3)]
    lbflow = [[-1000 for i in range(3)] for j in range(3)]
    Flow = m.addVars(nodes,nodes,lb=lbflow, ub = branch_limits, vtype = GRB.CONTINUOUS,name=flow_names)

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

    # add objective

    m.setObjective(quicksum(gen_costs[i] * P_gen[i] for i in range(nodes)))


    return m


def solve_opf_model(m):

    return m.optimize()

m = create_opf_model(nodes,gen_costs,branch_limits,B,gen_upper_bounds)

m.optimize()
print(m.getVars())