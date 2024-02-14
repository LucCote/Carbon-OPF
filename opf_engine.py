# simplified: one generator and one load at each node

import numpy as np
import gurobipy as gp
from gurobipy import GRB
from gurobipy import quicksum
from gurobipy import LinExpr


# dummy parameters

nodes = 3
P_load = np.array([-1,0,-3])

branch_limits = np.array([[0,10,10],[10,0,10],[10,10,0]])

branch_reactance = np.array([[1,0.7,1],[0.7,1,1],[1,0.7,1]])
B = 1/branch_reactance

gen_upper_bounds = np.array([2,2,3]) # generator capacity
gen_costs = np.array([10,15,12]) # gen cost per ouput

def create_opf_model(nodes,gen_costs,branch_limits,B,gen_upper_bounds):
    
    m = gp.Model()

    # add variables

    P_gen = m.addVars(nodes,ub = gen_upper_bounds, vtype= GRB.CONTINUOUS)
    voltage_angle = m.addVars(nodes,vtype= GRB.CONTINUOUS)
    Flow = m.addVars(nodes,nodes,ub = branch_limits, vtype = GRB.CONTINUOUS)

    # add constraints

    # Net flow = load + gen at each node

    m.addConstrs(quicksum(Flow[j,i] for j in range(nodes)) 
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