# IMPORTS
# External modules
from ETGEMs_function_PMI import *
import pandas as pd
import cobra
import ast
from cobra.io import write_sbml_model
from numpy import *
import copy
import math

# definate the function to simulate metabolism and interaction
# model_list: the model treated by Get_Concretemodel_Need_Data
# biomass_list: the biomass reaction in model
# growth_list: the growth rate of strains
# breed_list: the reciprocal of the number of breeding generations per time interval
# parameter_list: the list of the total number of enzymes, the upper bound for substrate input reaction flux (10 mmol/h/gDW) and the maximum value minus the minimum value of reaction thermodynamic driving force
# public_metabolism: the public metabolites
# vmax: the vmax values
# km: the km values
# deta_t: the time interval of iterations (h)
# deta_s: the grid length (cm)
# micro_distribute_threshold: the threshold to determine when the iteration stops
# length: the length of the whole one dimensional space (cm)
# D: the diffusion coefficient (cm^2/h)

def time_space_metabolism_state(model_list, biomass_list, growth_list, breed_list, parameter_list, public_metabolism, vmax, km, deta_t, deta_s, micro_distribute_threshold, length, D, metabolism_t):
    
    number_cell_side = length // deta_s
    number_cell_side = int(number_cell_side)
    public_metabolism_name_list = []
    public_metabolism_concentration_list = []
    for i in public_metabolism['metabolism']:
        public_metabolism_name_list.append(i)
    for j in public_metabolism['concentration']:
        public_metabolism_concentration_list.append(j)
    public_metabolism = dict(zip(public_metabolism_name_list, public_metabolism_concentration_list))
    
    number_model = len(model_list)
    
    k_m = {}
    v_max = {}
    for i in range(number_model):
        for j in range(len(km['reactions_strain'+str(i)])):
            k_m[(i, km['reactions_strain'+str(i)][j])] = km['km_strain'+str(i)][j]
            v_max[(i, vmax['reactions_strain'+str(i)][j])] = vmax['vmax_strain'+str(i)][j]
    
    
    number_public_metabolism = len(public_metabolism)
    distribute_micro_list = {}
    distribute_public_metabolism_list = {}
    distribute_lb_list = {}
    public_metabolism_r_list = []
    #set the initial distribution of straints
    for i in range(number_model):
        distribute_micro = [10]*number_cell_side
        distribute_micro_list.update({i: distribute_micro})
    print(distribute_micro_list)
    #set the initial distribution of metabolites
    for i in public_metabolism_name_list:
        distribute_public_metabolism = multiply(np.mat(np.ones(number_cell_side)), public_metabolism[i])
        distribute_public_metabolism_list.update({i: distribute_public_metabolism})
    #calculate the upperbounds of uptake reactions for public metabolites
    public_reaction_ub_list = {}
    public_reaction_list = {}
    for i in range(number_model):
        public_reaction_ub = {}
        public_reaction = []
        for rea in model_list[i]['reaction_list']:
            if 'EX_' in rea:
                for j in public_metabolism_name_list:
                    try:
                        model_list[i]['coef_matrix'][(j,rea)]
                    except:
                        pass
                    else:
                        ub = np.mat(np.ones(number_cell_side))
                        if model_list[i]['coef_matrix'][(j,rea)] > 0:
                            for m in range(number_cell_side):
                                ub[0,m] = v_max[(i,rea)]*distribute_public_metabolism_list[j][0,m]/(distribute_public_metabolism_list[j][0,m]+k_m[(i,rea)])
                        else:
                            try:
                                model_list[i]['ub_list'][rea]
                            except:
                                ub = ub * 1000
                            else:
                                ub = ub * model_list[i]['ub_list'][rea]
                        public_reaction_ub.update({rea: ub})
                        public_reaction.append(rea)
                        break
            public_reaction_ub_list[i] = public_reaction_ub
            public_reaction_list[i] = public_reaction

            
    ct = 0
    
    distribute_micro_list_t = {ct: distribute_micro_list}
    distribute_public_metabolism_list_t = {ct: distribute_public_metabolism_list}
    distribute_lb_list_t = {ct: distribute_lb_list}
    r = micro_distribute_threshold + 1
    r_threshold = [r]*5
    slip_r = np.mean(r_threshold[-5:])
    number = {}
    various = {}
    for i in range(number_model):
        number[i] = [np.mean(distribute_micro_list[i])]
        various[i] = [np.std(distribute_micro_list[i])/np.mean(distribute_micro_list[i])]
    
            
    # iterative simulation by slip_r
    while slip_r >= micro_distribute_threshold:
        ct += 1
        print(ct)
        distribute_micro_list_t[ct] = copy.deepcopy(distribute_micro_list_t[ct-1])
        micro_decrease = {}
        micro_increase = {}
        
        #simulate the cell wandering
        #micro_decrease: the decrease number of strains
        #micro_increase: the increase number of strains
        #micro_decrease_cell: the decrease number of strains in the current grid
        #micro_increase_fcell: the increase number of strains in the forward grid
        #micro_increase_bcell: the increase number of strains in the backward grid
        for m in range(number_cell_side):
            micro_decrease_cell = {}
            micro_increase_cell = {}
            for i in range(number_model):
                micro_decrease_cell[i] = 0
                micro_increase_cell[i] = 0
            micro_decrease[m] = micro_decrease_cell
            micro_increase[m] = micro_increase_cell
            
            
        if ct > 0:
            met = 'glc__D_e'
            threshold = 0.3
            for m in range(number_cell_side):
                #calculate the number of strains in the internal grids
                if 0<m<number_cell_side-1:
                    for i in range(number_model):
                        if distribute_micro_list_t[ct][i][m] > 0:
                            if distribute_public_metabolism_list[met][0,m-1] > distribute_public_metabolism_list[met][0,m]:
                                if distribute_public_metabolism_list[met][0,m+1] > distribute_public_metabolism_list[met][0,m]:
                                    micro_decrease_cell_ratio = min(threshold,(distribute_public_metabolism_list[met][0,m+1]/(distribute_micro_list_t[ct][i][m+1]+0.001)+distribute_public_metabolism_list[met][0,m-1]/(distribute_micro_list_t[ct][i][m-1]+0.001)-2*distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001))/(distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001)+0.001))
                                    micro_decrease_cell = micro_decrease_cell_ratio * distribute_micro_list_t[ct][i][m]
                                    micro_decrease_cell = int(micro_decrease_cell)
                                    micro_increase_fcell_ratio = (distribute_public_metabolism_list[met][0,m-1]/(distribute_micro_list_t[ct][i][m-1]+0.001)-distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001))/(distribute_public_metabolism_list[met][0,m+1]/(distribute_micro_list_t[ct][i][m+1]+0.001)+distribute_public_metabolism_list[met][0,m-1]/(distribute_micro_list_t[ct][i][m-1]+0.001)-2*distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001))
                                    micro_increase_fcell = micro_increase_fcell_ratio * micro_decrease_cell
                                    micro_increase_fcell = int(micro_increase_fcell)
                                    micro_increase_bcell = micro_decrease_cell-micro_increase_fcell
                                    micro_decrease[m][i] = micro_decrease[m][i] + micro_decrease_cell
                                    micro_increase[m-1][i] = micro_increase[m-1][i] + micro_increase_fcell
                                    micro_increase[m+1][i] = micro_increase[m+1][i] + micro_increase_bcell
                                else:
                                    micro_decrease_cell_ratio = min(threshold,(distribute_public_metabolism_list[met][0,m-1]/(distribute_micro_list_t[ct][i][m-1]+0.001)-distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001))/(distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001)+0.001))
                                    micro_decrease_cell = micro_decrease_cell_ratio * distribute_micro_list_t[ct][i][m]
                                    micro_decrease_cell = int(micro_decrease_cell)
                                    micro_increase_fcell_ratio = min(threshold,(distribute_public_metabolism_list[met][0,m-1]/(distribute_micro_list_t[ct][i][m-1]+0.001)-distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001))/(distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001)+0.001))
                                    micro_increase_fcell = micro_increase_fcell_ratio * distribute_micro_list_t[ct][i][m]
                                    micro_increase_fcell = int(micro_increase_fcell)
                                    micro_decrease[m][i] = micro_decrease[m][i] + micro_decrease_cell
                                    micro_increase[m-1][i] = micro_increase[m-1][i] + micro_increase_fcell
                            elif distribute_public_metabolism_list[met][0,m+1] > distribute_public_metabolism_list[met][0,m]:
                                micro_decrease_cell_ratio = min(threshold,(distribute_public_metabolism_list[met][0,m+1]/(distribute_micro_list_t[ct][i][m+1]+0.001)-distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001))/(distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001)+0.001))
                                micro_decrease_cell = micro_decrease_cell_ratio * distribute_micro_list_t[ct][i][m]
                                micro_decrease_cell = int(micro_decrease_cell)
                                micro_increase_bcell_ratio = min(threshold,(distribute_public_metabolism_list[met][0,m+1]/(distribute_micro_list_t[ct][i][m+1]+0.001)-distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001))/(distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001)+0.001))
                                micro_increase_bcell = micro_increase_bcell_ratio * distribute_micro_list_t[ct][i][m]
                                micro_increase_bcell = int(micro_increase_bcell)
                                micro_decrease[m][i] = micro_decrease[m][i] + micro_decrease_cell
                                micro_increase[m+1][i] = micro_increase[m+1][i] + micro_increase_bcell
                #calculate the number of strains in the first grid
                elif m == 0:
                    for i in range(number_model):
                        if distribute_micro_list_t[ct][i][m] > 0:
                            if distribute_public_metabolism_list[met][0,m+1] > distribute_public_metabolism_list[met][0,m]:
                                micro_decrease_cell_ratio = min(threshold,(distribute_public_metabolism_list[met][0,m+1]/(distribute_micro_list_t[ct][i][m+1]+0.001)-distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001))/(distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001)+0.001))
                                micro_decrease_cell = micro_decrease_cell_ratio * distribute_micro_list_t[ct][i][m]
                                micro_decrease_cell = int(micro_decrease_cell)
                                micro_increase_bcell_ratio = min(threshold,(distribute_public_metabolism_list[met][0,m+1]/(distribute_micro_list_t[ct][i][m+1]+0.001)-distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001))/(distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001)+0.001))
                                micro_increase_bcell = micro_increase_bcell_ratio * distribute_micro_list_t[ct][i][m]
                                micro_increase_bcell = int(micro_increase_bcell)
                                micro_decrease[m][i] = micro_decrease[m][i] + micro_decrease_cell
                                micro_increase[m+1][i] = micro_increase[m+1][i] + micro_increase_bcell
                #calculate the number of strains in the last grid
                elif m == number_cell_side-1:
                    for i in range(number_model):
                        if distribute_micro_list_t[ct][i][m] > 0:
                            if distribute_public_metabolism_list[met][0,m-1] > distribute_public_metabolism_list[met][0,m]:
                                micro_decrease_cell_ratio = min(threshold,(distribute_public_metabolism_list[met][0,m-1]/(distribute_micro_list_t[ct][i][m-1]+0.001)-distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001))/(distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001)+0.001))
                                micro_decrease_cell = micro_decrease_cell_ratio * distribute_micro_list_t[ct][i][m]
                                micro_decrease_cell = int(micro_decrease_cell)
                                micro_increase_fcell_ratio = min(threshold,(distribute_public_metabolism_list[met][0,m-1]/(distribute_micro_list_t[ct][i][m-1]+0.001)-distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001))/(distribute_public_metabolism_list[met][0,m]/(distribute_micro_list_t[ct][i][m]+0.001)+0.001))
                                micro_increase_fcell = micro_increase_fcell_ratio * distribute_micro_list_t[ct][i][m]
                                micro_increase_fcell = int(micro_increase_fcell)
                                micro_decrease[m][i] = micro_decrease[m][i] + micro_decrease_cell
                                micro_increase[m-1][i] = micro_increase[m-1][i] + micro_increase_fcell
        
            #calculate the number of strains after wandering
            for m in range(number_cell_side):
                for i in range(number_model):
                    distribute_micro_list_t[ct][i][m] = distribute_micro_list_t[ct][i][m] - micro_decrease[m][i]
                    distribute_micro_list_t[ct][i][m] = distribute_micro_list_t[ct][i][m] + micro_increase[m][i]
                    distribute_micro_list_t[ct][i][m] = max(0, distribute_micro_list_t[ct][i][m])
                            
                            
        
        #simulate the substrate diffusion by Fick's second law
        for m in range(number_cell_side): 
            if 0<m<number_cell_side-1:
                for met in public_metabolism_name_list:
                    distribute_public_metabolism_list[met][0,m] = distribute_public_metabolism_list[met][0,m] + ((distribute_public_metabolism_list[met][0,m+1]-distribute_public_metabolism_list[met][0,m])/deta_s-(distribute_public_metabolism_list[met][0,m]-distribute_public_metabolism_list[met][0,m-1])/deta_s)/deta_s*D*deta_t
                    if distribute_public_metabolism_list[met][0,m] < 0:
                        print('Warning: the D is too high!')
                    distribute_public_metabolism_list[met][0,m] = max(0,distribute_public_metabolism_list[met][0,m])
            if m == 0:
                distribute_public_metabolism_list[met][0,m] = distribute_public_metabolism_list[met][0,m] + (distribute_public_metabolism_list[met][0,m+1]-distribute_public_metabolism_list[met][0,m])/deta_s/deta_s*D*deta_t
                distribute_public_metabolism_list[met][0,m] = max(0,distribute_public_metabolism_list[met][0,m])
            if m == number_cell_side-1:
                distribute_public_metabolism_list[met][0,m] = distribute_public_metabolism_list[met][0,m] + (distribute_public_metabolism_list[met][0,m-1]-distribute_public_metabolism_list[met][0,m])/deta_s/deta_s*D*deta_t
                distribute_public_metabolism_list[met][0,m] = max(0,distribute_public_metabolism_list[met][0,m])
        
        #simulate the metabolism by ETMs
        biomass_value_list = {}
        number_model_range = []
        for m in range(number_cell_side):
            if m%2 == 0:
                number_model_range.append([0,1])
            elif m%2 == 1:
                number_model_range.append([1,0])
        for m in range(number_cell_side):
            B_value_list = []
            
            
            biomass_value_micro = {}
            
            
            for i in number_model_range[m]:
                if distribute_micro_list_t[ct][i][m] > 0:
                    public_metabolism_flux_list = {}
                    for j in public_metabolism_name_list:
                        public_metabolism_flux_list.update({j: 0})
                    for j in public_reaction_list[i]:
                        model_list[i]['ub_list'][j] = int(public_reaction_ub_list[i][j][0,m])
            
                    biomass_id = biomass_list[i]
                    E_total=parameter_list[i][0]
                    #set the carbon source to glucose
                    substrate_name='EX_glc__D_e_reverse'
                    substrate_value=parameter_list[i][1]
                    biomass_value=growth_list[i]
                    K_value=parameter_list[i][2]

                    try:
                        MDF_Calculation(model_list[i],biomass_value,biomass_id,substrate_name,substrate_value,K_value,E_total,'gurobi')
                    except:
                        pass
                    else:
                        #calculate the MDF values of metabolic networks
                        biomass_value_micro[i] = biomass_value
                        B_value=MDF_Calculation(model_list[i],biomass_value,biomass_id,substrate_name,substrate_value,K_value,E_total,'gurobi')
                        B_value_list.append(B_value)
                        #calculate the biomass yield under the MDF value
                        obj_name=biomass_list[i]
                        obj_target='maximize'
                        if i == 0:
                            max_biomass_under_mdf=Max_Growth_Rate_Calculation0(model_list[i],obj_name,obj_target,substrate_name,substrate_value,K_value,E_total,B_value,'gurobi')
                        elif i == 1:
                            max_biomass_under_mdf=Max_Growth_Rate_Calculation1(model_list[i],obj_name,obj_target,substrate_name,substrate_value,K_value,E_total,B_value,'gurobi')
                        biomass_value=max_biomass_under_mdf*0.9
                        '''
                        if i == 0:
                            min_E=Min_Enzyme_Cost_Calculation0(model_list[i],biomass_value,biomass_id,substrate_name,substrate_value,K_value,E_total,B_value,'gurobi')
                        elif i == 1:
                            min_E=Min_Enzyme_Cost_Calculation1(model_list[i],biomass_value,biomass_id,substrate_name,substrate_value,K_value,E_total,B_value,'gurobi')
                        '''
                        #calculate the minimum value of the sum of the fluxes
                        if i == 0:
                            [min_V,Concretemodel]=Min_Flux_Sum_Calculation0(model_list[i],biomass_value,biomass_id,substrate_name,substrate_value,K_value,E_total,B_value,'gurobi')
                        elif i == 1:
                            [min_V,Concretemodel]=Min_Flux_Sum_Calculation1(model_list[i],biomass_value,biomass_id,substrate_name,substrate_value,K_value,E_total,B_value,'gurobi')
            
                        #translating the results of ETMs into a usable form
                        model=model_list[i]['model']
                        reaction_kcat_MW=model_list[i]['reaction_kcat_MW']
                        reaction_g0=model_list[i]['reaction_g0']
                        coef_matrix=model_list[i]['coef_matrix']
                        metabolite_list=model_list[i]['metabolite_list']
                        use_result = Get_Results_Thermodynamics(model,Concretemodel,reaction_kcat_MW,reaction_g0,coef_matrix,metabolite_list)
            
                        #simulate the fluxes of the public metabolites
                        for rea in public_reaction_list[i]:
                            for met in public_metabolism_name_list:
                                try:
                                    model_list[i]['coef_matrix'][(met,rea)]  
                                except:
                                    pass
                                else:
                                    public_metabolism_flux_list[met] = public_metabolism_flux_list[met] + model_list[i]['coef_matrix'][(met,rea)]* use_result['flux'][rea]
                        
                        #get the metabolism of strains
                        if ct == metabolism_t:
                            if m == 0:
                                if i == 0:
                                    metabolism_strains = {}
                                metabolism_strains[i] = {}
                                for rea in model_list[i]['reaction_list']:
                                    metabolism_strains[i][rea] = use_result['flux'][rea]
                                print (metabolism_strains)
                        
                        #simulate the distribution of the public metabolites
                        distribute_public_metabolism_ori = copy.deepcopy(distribute_public_metabolism_list)
                        distribute_public_metabolism_re = copy.deepcopy(distribute_public_metabolism_ori)
                        for met in public_metabolism_name_list:
                            distribute_public_metabolism_list[met][0,m] = distribute_public_metabolism_list[met][0,m] - public_metabolism_flux_list[met] * deta_t * distribute_micro_list_t[ct][i][m]
                            distribute_public_metabolism_ori[met][0,m] = copy.deepcopy(distribute_public_metabolism_list[met][0,m])
                            if distribute_public_metabolism_list[met][0,m] < 0:
                                print(met+'_1: ', distribute_public_metabolism_list[met][0,m])
                                
                        #simulate the multiplication and death rates of strains and the distribution of public metabolites after multiplication or death        
                        death_rate = 0
                        birth_rate = 1
                        for rea in public_reaction_list[i]:
                            if 'reverse' not in rea:
                                met = rea[3::]
                                if met in public_metabolism_name_list:
                            
                                    intracellular_c = use_result['met_concentration'][rea].lstrip(';')
                                    intracellular_c = intracellular_c.replace(';',',"')
                                    intracellular_c = intracellular_c.replace(' :','" :')
                                    intracellular_c = '{"' + intracellular_c + '}'
                                    intracellular_c_dict = ast.literal_eval(intracellular_c)
                                
                                    bd_rate_met = (distribute_public_metabolism_list[met][0,m] - 0.1)/intracellular_c_dict[met]/distribute_micro_list_t[ct][i][m]
                                    if bd_rate_met < 0:
                                        death_rate_lb = (public_metabolism_flux_list[met] * deta_t * distribute_micro_list_t[ct][i][m] - distribute_public_metabolism_re[met][0,m] + 0.1)/(public_metabolism_flux_list[met] * deta_t * distribute_micro_list_t[ct][i][m] + intracellular_c_dict[met] * distribute_micro_list_t[ct][i][m])
                                        death_rate = max(death_rate, death_rate_lb)
                                    else:
                                        birth_rate = min(birth_rate, bd_rate_met)
                        death_rate = min(death_rate, 1)
                        birth_rate = min(birth_rate,1)
                        if death_rate > 0:
                            birth_rate = 0
                            for met in public_metabolism_name_list:
                                distribute_public_metabolism_list[met][0,m] = distribute_public_metabolism_list[met][0,m] + public_metabolism_flux_list[met] * deta_t * math.ceil(distribute_micro_list_t[ct][i][m] * death_rate)
                                distribute_public_metabolism_ori[met][0,m] = copy.deepcopy(distribute_public_metabolism_list[met][0,m])
                    
                        for rea in public_reaction_list[i]:
                            if 'reverse' not in rea:
                                met = rea[3::]
                                if met in public_metabolism_name_list:
                            
                                    intracellular_c = use_result['met_concentration'][rea].lstrip(';')
                                    intracellular_c = intracellular_c.replace(';',',"')
                                    intracellular_c = intracellular_c.replace(' :','" :')
                                    intracellular_c = '{"' + intracellular_c + '}'
                                    intracellular_c_dict = ast.literal_eval(intracellular_c)
                                
                                    distribute_public_metabolism_list[met][0,m] = distribute_public_metabolism_ori[met][0,m] - intracellular_c_dict[met]*math.floor(distribute_micro_list_t[ct][i][m]*birth_rate/breed_list[i]) + intracellular_c_dict[met]*math.ceil(death_rate*distribute_micro_list_t[ct][i][m])
                                    if distribute_public_metabolism_list[met][0,m] < 0:
                                        print(met+'_2: ', distribute_public_metabolism_list[met][0,m])
                                        distribute_public_metabolism_list[met][0,m] = 0
                        
                        #simulate the distribution of strains after multiplication or death
                        distribute_micro_increase = math.floor(distribute_micro_list_t[ct][i][m]*birth_rate/breed_list[i])
                        distribute_micro_decrease = math.ceil(distribute_micro_list_t[ct][i][m]*death_rate)
                        distribute_micro_list_t[ct][i][m] = max(0,distribute_micro_list_t[ct][i][m] + distribute_micro_increase - distribute_micro_decrease)
            
            #calculate the upperbounds of nutrient uptake rates after substrate diffusion, cell wandering, metabolism, multiplication and death
            biomass_value_list[m] = biomass_value_micro
            
            for i in range(number_model):
                for rea in model_list[i]['reaction_list']:
                    if 'EX_' in rea:
                        for j in public_metabolism_name_list:
                            try:
                                model_list[i]['coef_matrix'][(j,rea)]
                            except:
                                pass
                            else:
                                if model_list[i]['coef_matrix'][(j,rea)] > 0:
                                    public_reaction_ub_list[i][rea][0, m] = v_max[(i,rea)]*distribute_public_metabolism_list[j][0,m]/(distribute_public_metabolism_list[j][0,m]+k_m[(i,rea)])
        '''
        print('metabolite list: ', distribute_public_metabolism_list['glc__D_e'])
        print('microbial distribution: ', distribute_micro_list_t[ct])
        print('biomass values: ', biomass_value_list)
        
        print(ct, distribute_micro_list_t[ct])
        print(ct-1, distribute_micro_list_t[ct-1])
        '''
        # calculate the mean number and uniformity of distribution of strains at this iteration
        for i in range(number_model):
            strain_number = np.mean(distribute_micro_list_t[ct][i])
            strain_various = np.std(distribute_micro_list_t[ct][i])
            number[i].append(strain_number)
            various[i].append(strain_various)
            print('strain_number: ', i, strain_number)
            print('strain_various: ', i, strain_various)
        
        # calculate the slip_r at this iteration
        if ct > 1:
            r = 0
            for i in range(number_model):
                for m in range(number_cell_side):
                    if distribute_micro_list_t[ct-1][i][m] > 0:
                        r = r + ((distribute_micro_list_t[ct][i][m]-distribute_micro_list_t[ct-1][i][m])/(distribute_micro_list_t[ct-1][i][m]))**2
                    else:
                        r = r + (distribute_micro_list_t[ct][i][m])**2
        '''
        print('micro_distribute_difference: ', r)
        '''
        r_threshold.append(r)
        slip_r = np.mean(r_threshold[-5:])
        fd_r_threshold = r_threshold[5:]
        print('fd_r_threshold: ', fd_r_threshold)
        print('slip_r: ', slip_r)
    return(distribute_micro_list_t, number, various, metabolism_strains)