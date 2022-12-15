############################################################################

# Created by: Prof. Valdecy Pereira, D.Sc.
# UFF - Universidade Federal Fluminense (Brazil)
# email:  valdecy.pereira@gmail.com
# pyMissingAHP

# Citation: 
# PEREIRA, V. (2022). Project: pyMissingAHP, GitHub repository: <https://github.com/Valdecy/pyMissingAHP>

############################################################################

# Required Libraries
import itertools
import math
import matplotlib.pyplot as plt
plt.style.use('bmh')
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import scipy.stats as stats

from fractions import Fraction
from pyMissingAHP.util.ahp import ahp_method
from pyMissingAHP.util.fuzzy_ahp import fuzzy_ahp_method
from pyMissingAHP.util.ga import genetic_algorithm

################################################################################

# pyMissingAHP Class
class load_ahp():
    def __init__(self, original, min_limit = [], max_limit = [], wdm = 'geometric', scale = 'discrete', fuzzy_scale = False, custom_fuzzy_scale = [], custom_rank = []):
      self.alpha   = 0
      self.wd      = wdm 
      self.f_flag  = fuzzy_scale
      self.s_flag  = scale
      self.order   = custom_rank
      if (len(self.order) < len(original)):
          for _ in range(0, len(original) - len(self.order)):
              self.order.append(-1)
      if (len(self.order) > len(original)):
          self.order = self.order[:len(original)]
      if (self.f_flag == False):
        self.dataset = original
        if (len(min_limit) > 0):
          self.dataset_min = min_limit
        else:
          self.dataset_min = np.zeros(( self.dataset.shape))
          for i in range(0, self.dataset.shape[0]):
            for j in range(1, self.dataset.shape[1]):
              if (j > i):
                self.dataset_min[i, j] = 1/9
        if (len(max_limit) > 0):
          self.dataset_max = max_limit
        else:
          self.dataset_max = np.zeros(( self.dataset.shape))
          for i in range(0, self.dataset.shape[0]):
            for j in range(1, self.dataset.shape[1]):
              if (j > i):
                self.dataset_max[i, j] = 9
      elif (self.f_flag == True):
        self.dataset = np.zeros((len(original), len(original[0])), dtype = object)
        for i in range(0, len(original)):
            for j in range(0, len(original[i])):
              self.dataset[i, j] = original[i][j]
        if (len(min_limit) > 0):
          self.dataset_min = np.zeros((len(min_limit), len(min_limit[0])), dtype = object)
          for i in range(0, len(min_limit)):
              for j in range(0, len(min_limit[i])):
                self.dataset_min[i, j] = min_limit[i][j]
        else:
          self.dataset_min = np.zeros((len(original), len(original[0])), dtype = object)
          for i in range(0, self.dataset.shape[0]):
            for j in range(1, self.dataset.shape[1]):
              if (j > i):
                self.dataset_min[i, j] = (1/9, 1/9, 1/9)
        if (len(max_limit) > 0):
          self.dataset_max = np.zeros((len(max_limit), len(max_limit[0])), dtype = object)
          for i in range(0, len(max_limit)):
              for j in range(0, len(max_limit[i])):
                self.dataset_max[i, j] = max_limit[i][j]
        else:
          self.dataset_max = np.zeros((len(original), len(original[0])), dtype = object)
          for i in range(0, self.dataset.shape[0]):
            for j in range(1, self.dataset.shape[1]):
              if (j > i):
                self.dataset_max[i, j] = (9, 9, 9)
      if (self.f_flag == False):
        self.saaty_scale = [-1, 1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        self.saaty_strg  = [str(Fraction(item).limit_denominator()) if item < 1 else str(item) for item in self.saaty_scale]
      elif (self.f_flag == True and len(custom_fuzzy_scale) == 0):
        self.saaty_scale = [ ( -1,  -1,  -1), (1/9, 1/9, 1/9), (1/9, 1/8, 1/7), (1/8, 1/7, 1/6), (1/7, 1/6, 1/5), (1/6, 1/5, 1/4),(1/5, 1/4, 1/3), (1/4, 1/3, 1/2), (1/3, 1/2, 1), (1, 1, 1), (1, 2, 3), (2, 3, 4), (3, 4, 5), (4, 5, 6), (5, 6, 7), (6, 7, 8), (7, 8, 9), (9, 9, 9) ]
        self.saaty_strg  = []
        for fuzzy_number in self.saaty_scale:
          srt = '('
          for item in fuzzy_number:
            if (item < 1):
              srt = srt + str(Fraction(item).limit_denominator()) + ', '
            else:
              srt = srt + str(item) + ', '
          srt = srt[:-2]+')'
          self.saaty_strg.append(srt)
      elif (self.f_flag == True and len(custom_fuzzy_scale) > 0):
        self.saaty_scale = [item for item in custom_fuzzy_scale]
        if (( -1,  -1,  -1) not in self.saaty_scale):
            self.saaty_scale.insert(0, ( -1,  -1,  -1))
        self.saaty_strg  = []
        for fuzzy_number in self.saaty_scale:
          srt = '('
          for item in fuzzy_number:
            if (item < 1):
              srt = srt + str(Fraction(item).limit_denominator()) + ', '
            else:
              srt = srt + str(item) + ', '
          srt = srt[:-2]+')'
          self.saaty_strg.append(srt)
      self.minv      = []
      self.maxv      = []
      self.miss_pos  = []
      self.var_pos   = []
      v              = -1
      for i in range(0, self.dataset.shape[0]):
        for j in range(i, self.dataset.shape[1]):
          #if (j > i):
          if (self.dataset_min[i,j] != self.dataset_max[i,j]):
              v = v + 1
              self.miss_pos.append([i,j])
              self.var_pos.append(v)  
      for pair in self.miss_pos:
          i, j = pair
          if (self.s_flag == 'discrete' and self.f_flag == False):
              self.minv.append(0.0)
              self.maxv.append(1.0)
          elif (self.s_flag != 'discrete' and self.f_flag == False):
              self.minv.append(self.dataset_min[i,j])
              self.maxv.append(self.dataset_max[i,j])
          elif (self.f_flag == True):
              self.minv.append(0.0)
              self.maxv.append(1.0)

    ################################################################################

    # Function: Load AHP Parameters
    def run_ahp(self):
        flag = -1 in self.dataset 
        if (flag == False):
            self.weights, self.rc = ahp_method(self.dataset, self.wd)
            rank                   = self.rank_descending(self.weights)
            for i in range(0, self.weights.shape[0]):
              print('w(g'+str(i+1)+'): ', f'{self.weights[i]:.3f}', '; rank(g'+str(i+1)+'): ', rank[i])
            print('')
            if (self.rc > 0.10):
              print('RC: ' +  f'{self.rc[i]:.4f}', ' The solution is inconsistent, the pairwise comparisons must be reviewed')
            else:
              print('RC: ' + f'{self.rc[i]:.4f}', ' The solution is consistent') 
        self.plot_environment(self.dataset)
        return

    # Function: Load Fuzzy AHP Parameters
    def run_fuzzy_ahp(self):
        flag = False
        for i in range(0, self.dataset.shape[0]):
            for j in range(i, self.dataset.shape[1]):
                a, b, c = self.dataset[i, j]
                if (a == -1 or b == -1 or c == -1):
                    flag = True
                    break
        if (flag == False):
            self.weights, self.rc = fuzzy_ahp_method(self.dataset)
            rank                  = self.rank_descending(self.weights)
            for i in range(0, self.weights.shape[0]):
              print('w(g'+str(i+1)+'): ', f'{self.weights[i]:.3f}', '; rank(g'+str(i+1)+'): ', rank[i])
            print('')
            if (self.rc > 0.10):
              print('RC: ' + f'{self.rc[i]:.4f}', ' The solution is inconsistent, the pairwise comparisons must be reviewed')
            else:
              print('RC: ' + f'{self.rc[i]:.4f}', ' The solution is consistent') 
        self.plot_environment(self.dataset)
        return
    
    ################################################################################
    
    # Function: Load GA Parameters
    def run_ga(self, population_size, mutation_rate, list_of_functions, generations, mu, eta, elite, verbose):
        self.size = population_size
        self.m_r  = mutation_rate
        self.gen  = generations
        self.mu   = mu
        self.eta  = eta
        self.elt  = elite
        self.vbs  = verbose
        self.lof  = []
        self.qidx = [item for item in list_of_functions]
        for i in range(0, len(list_of_functions)):
            if   (list_of_functions[i] == 'f0' and self.f0 not in self.lof):
                self.lof.append(self.f0)
            elif (list_of_functions[i] == 'f1' and self.f1 not in self.lof):
                self.lof.append(self.f1)
        self.run()
        return

    # Function: Load GA_M Parameters
    def run_ga_m(self, population_size, mutation_rate, list_of_functions, generations, mu, eta, elite, verbose, step):
        self.size = population_size
        self.m_r  = mutation_rate
        self.gen  = generations
        self.mu   = mu
        self.eta  = eta
        self.elt  = elite
        self.vbs  = verbose
        self.lof  = []
        self.qidx = [item for item in list_of_functions]
        for i in range(0, len(list_of_functions)):
            if   (list_of_functions[i] == 'f0' and self.f0 not in self.lof):
                self.lof.append(self.f0)
            elif (list_of_functions[i] == 'f1' and self.f1 not in self.lof):
                self.lof.append(self.f1)
            elif (list_of_functions[i] == 'f0_f1' and self.f0_f1 not in self.lof):
                self.lof.append(self.f0_f1)
        cnt_max = int(1/step)
        cnt_tot = -1
        sol_m   = []
        print('')
        print('Starting Optimization...')
        print('')
        while self.alpha <= (1 + step):
            cnt_tot    = cnt_tot + 1
            print('Iteration ', cnt_tot,' of ',cnt_max,' (alpha: ', f'{abs(self.alpha):.3f}',', 1-alpha: ', f'{abs(1-self.alpha):.3f}',')')
            self.run()
            self.alpha = self.alpha + step
            sol_m.append(self.solution)
        self.solution = [item for item in sol_m]
        self.alpha    = 0  
        return
    
    ################################################################################
    
    # Function: Rank Decending (Adapted from: https://stackoverflow.com/questions/39059371/can-numpys-argsort-give-equal-element-the-same-rank)
    def rank_descending(self, x):
        _, inv = np.unique(-x, return_inverse = True, return_counts = False)
        return inv
    
    # Function: All Possible Discrete Combinations
    def brute_force(self, view_report = False, consistent_only = True):
        print('Considering Only Discrete Values')
        rang_elements      = []
        self.bf_solutions  = []
        report_lst         = []
        count              = 0
        for pair in self.miss_pos:
            i, j    = pair
            idx_min = self.saaty_scale.index(self.dataset_min[i,j])
            idx_max = self.saaty_scale.index(self.dataset_max[i,j])
            rang_elements.append( [self.saaty_scale[i] for i in range(0, len(self.saaty_scale)) if (i >= idx_min and i <= idx_max)] )
        combinations    = list(itertools.product(*rang_elements))      
        print('Total Number of Combinations: ',   len(combinations))
        print('Total Number of Missing Values: ', len(combinations[0]))
        dataset_   = np.array(self.dataset, copy = True) 
        rc_bf      = float('+inf') 
        for item in combinations:
            for k in range(0, len(item)):
                i, j           = self.miss_pos[k]
                dataset_[i, j] = item[k]
                if (self.f_flag == False):
                    dataset_[j, i] = 1/item[k]
                else:
                    a, b, c        = dataset_[i, j]
                    dataset_[j, i] = (1/c, 1/b, 1/a)  
            if (self.f_flag == False):
                w_, rc_ = ahp_method(dataset_, self.wd)
            else:
                w_, rc_ = fuzzy_ahp_method(dataset_)
            if (consistent_only == True):
                if (rc_ <= 0.1):
                    order_    = self.rank_descending(w_)
                    self.bf_solutions.append(np.array(dataset_, copy = True))
                    flat_list = [w_.tolist(), order_.tolist(), [rc_]]
                    flat_list = [item for sublist in flat_list for item in sublist]
                    report_lst.append(flat_list)
                    count     = count + 1
                    if (rc_bf > rc_):
                        rc_bf = rc_
            elif (consistent_only == False):
                 order_     = self.rank_descending(w_)
                 self.bf_solutions.append(np.array(dataset_, copy = True))
                 flat_list  = [w_.tolist(), order_.tolist(), [rc_]]
                 flat_list  = [item for sublist in flat_list for item in sublist]
                 report_lst.append(flat_list)
                 if (rc_ <= 0.1):
                     count = count + 1
                     if (rc_bf > rc_):
                         rc_bf = rc_     
        print('Total Number of Consistent Solutions: ', count)
        print('Minimum RC: ', f'{rc_bf:.4f}')
        col_names = []
        for i in range(0, len(w_)):
            col_names.append('w(g'+str(i+1)+')') 
        for i in range(0, len(order_)):
            col_names.append('rank(g'+str(i+1)+')') 
        col_names.append('Consistency')
        self.report_df = pd.DataFrame(report_lst, columns = col_names)
        if (view_report == True):
            print('')
            print(self.report_df)
        return self.bf_solutions, self.report_df
    
    ################################################################################
    
    # Function: Plot Environment
    def plot_environment(self, environment, plot_size = 10):
        keys_j   = self.saaty_scale[1:]
        values_j = self.saaty_strg[1:]
        dict_j   = dict(zip(keys_j, values_j))
        fig      = plt.figure(figsize = [plot_size, plot_size], facecolor = 'w')
        ax       = fig.add_subplot(111, xticks = range(environment.shape[1] + 1), yticks = range(environment.shape[0] + 1), position = [0.1, 0.1, 0.8, 0.8])
        plt.gca().invert_yaxis()
        ax.grid(color = 'k', linestyle = '-', linewidth = 1)
        ax.xaxis.set_tick_params(bottom = 'off', top   = 'off', labelbottom = 'off')
        ax.yaxis.set_tick_params(left   = 'off', right = 'off', labelleft   = 'off')
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        for i in range(0, environment.shape[0]):
            for j in range(0, environment.shape[1]):
                if ( environment[i, j] in dict_j):
                    ax.annotate(dict_j[environment[i, j]] , xy = (0.4 + j, 0.55 + i), fontsize = 10, fontweight = 'bold')
                elif(self.f_flag == True):
                    ax.annotate(environment[i, j] , xy = (0.4 + j, 0.55 + i), fontsize = 10, fontweight = 'bold')
                else:
                    ax.annotate(f'{environment[i, j]:.3f}' , xy = (0.4 + j, 0.55 + i), fontsize = 10, fontweight = 'bold')
                if (i == j):
                  ax.annotate(dict_j[environment[i, j]] , xy = (0.4 + j, 0.55 + i), fontsize = 10, fontweight = 'bold', color = 'w')
                  black_stone = mpatches.Rectangle( (j, i), 1, 1, linewidth = 1, edgecolor = 'k', facecolor = 'k', clip_on = False)
                  ax.add_patch(black_stone)
                if ( ( [i, j] in self.miss_pos ) and j > i ):
                    if ( environment[i, j] in dict_j):
                        ax.annotate(dict_j[environment[i, j]] , xy = (0.4 + j, 0.55 + i), fontsize = 10, fontweight = 'bold', color = 'k')
                    elif(self.f_flag == True):
                        ax.annotate(environment[i, j] , xy = (0.4 + j, 0.55 + i), fontsize = 10, fontweight = 'bold')
                    else:
                        ax.annotate(f'{environment[i, j]:.3f}' , xy = (0.4 + j, 0.55 + i), fontsize = 10, fontweight = 'bold', color = 'k')
                    yellow_stone = mpatches.Rectangle( (j, i), 1, 1, linewidth = 1, edgecolor = 'k', facecolor = 'yellow', clip_on = False)
                    ax.add_patch(yellow_stone)
        return

    # Function: Scatter Plot 
    def plot_scatter(self, view = 'browser', pf = False):
        def pareto_front_points(pts, pf_min = True):
            def pareto_front(pts, pf_min):
                pf = np.zeros(pts.shape[0], dtype = np.bool_)
                for i in range(0, pts.shape[0]):
                    cost = pts[i, :]
                    if (pf_min == True):
                        g_cost = np.logical_not(np.any(pts > cost, axis = 1))
                        b_cost = np.any(pts < cost, axis = 1)
                    else:
                        g_cost = np.logical_not(np.any(pts < cost, axis = 1))
                        b_cost = np.any(pts > cost, axis = 1)
                    dominated = np.logical_and(g_cost, b_cost)
                    if  (np.any(pf) == True):
                        if (np.any(np.all(pts[pf] == cost, axis = 1)) == True):
                            continue
                    if not (np.any(dominated[:i]) == True or np.any(dominated[i + 1 :]) == True):
                        pf[i] = True
                return pf
            idx     = np.argsort(((pts - pts.mean(axis = 0))/(pts.std(axis = 0) + 1e-7)).sum(axis = 1))
            pts     = pts[idx]
            pf      = pareto_front(pts, pf_min)
            pf[idx] = pf.copy()
            return pf
        if (view == 'browser'):
            pio.renderers.default = 'browser'
        data     = []
        inc_list = [ 'Index: '+str(item)+'<br>'+'Solution: Inconsistent' for item in list(self.ind_incon.index)] 
        con_list = [ 'Index: '+str(item)+'<br>'+'Solution: Consistent'   for item in list(self.ind_con.index  )] 
        if (pf == True):
            pts           = self.indicators.iloc[:,:-1]
            pts['f0(MI)'] = pts['f0(MI)'].astype(float)
            pts['f1(KT)'] = pts['f1(KT)'].astype(float)
            pts.drop_duplicates()
            pf            = pareto_front_points(pts.values, pf_min = True)
            front         = pts.iloc[pf, :]
            x, y          = zip(*sorted(zip(front.iloc[:,0], front.iloc[:,1])))
            f_trace       = go.Scatter(x    = x, 
                                       y    = y,
                                       mode = 'lines',
                                       line = dict(color = 'black', width = 1)
                                       )
            data.append(f_trace)
        if (len(con_list) > 0):
          s_trace = go.Scatter(
                              x         = self.ind_con.iloc[:, 0], # ~self.ind_con.index.isin(self.rank_idx)
                              y         = self.ind_con.iloc[:, 1], # ~self.ind_con.index.isin(self.rank_idx)
                              opacity   = 0.85,
                              mode      = 'markers+text',
                              marker    = dict(symbol = 'circle-dot', size = 8, color = 'red'),
                              hovertext = con_list,
                              name      = ''
                              )
          data.append(s_trace)
        if (len(inc_list) > 0):
          n_trace = go.Scatter(
                              x         = self.ind_incon.iloc[:, 0],
                              y         = self.ind_incon.iloc[:, 1],
                              opacity   = 0.5,
                              mode      = 'markers+text',
                              marker    = dict(symbol = 'circle-dot', size = 10, color = 'purple'),
                              hovertext = inc_list,
                              name      = ''
                              )
          data.append(n_trace)
        layout  = go.Layout(showlegend   = False,
                            hovermode    = 'closest',
                            margin       = dict(b = 10, l = 5, r = 5, t = 10),
                            plot_bgcolor = 'white',
                            xaxis        = dict(  showgrid       = True, 
                                                  zeroline       = False, 
                                                  showticklabels = True, 
                                                  title          = self.indicators.columns[0],
                                                  tickmode       = 'array', 
                                                  gridcolor      = 'grey',
                                                  spikedash      = 'solid',
                                                  spikecolor     = 'blue',
                                                  spikethickness = 2
                                              ),
                            yaxis        = dict(  showgrid       = True, 
                                                  zeroline       = False, 
                                                  showticklabels = True,
                                                  title          = self.indicators.columns[1],
                                                  tickmode       = 'array', 
                                                  gridcolor      = 'grey',
                                                  spikedash      = 'solid',
                                                  spikecolor     = 'blue',
                                                  spikethickness = 2
                                                )
                            )
        fig_aut = go.Figure(data = data, layout = layout)
        fig_aut.update_traces(textfont_size = 10, textfont_color = 'white') 
        fig_aut.update_layout(autotypenumbers = 'convert types')
        fig_aut.show() 
        return

    ################################################################################
    
    # Functions: Discretize
    def judgement_discretize(self, variable, i, j):
        min_val = self.dataset_min[i,j]
        max_val = self.dataset_max[i,j]
        if (self.f_flag == False):
            judgs = [i for i in self.saaty_scale[1:] if i >= min_val and i <= max_val]
        else:
            min_val = self.saaty_scale[1:].index(min_val)
            max_val = self.saaty_scale[1:].index(max_val)
            judgs   = [i for i in self.saaty_scale[1:] if self.saaty_scale[1:].index(i) >= min_val and self.saaty_scale[1:].index(i) <= max_val]
        ranges  = [1/(len(judgs))*i for i in range(0, len(judgs))] 
        lower   = []
        upper   = []
        value   = 1
        if (self.f_flag == True):
            value = (1, 1, 1)
        for j in range(0, len(ranges)-1):
            lower.append(ranges[j])
            upper.append(ranges[j+1])
        lower.append(ranges[-1])
        upper.append(1.0)
        for j in range(0, len(judgs)):
            if (variable >= lower[j] and variable < upper[j]):
                value = judgs[j]
        return value
    
    # Functions: Decode Solution
    def convert_solution(self, sol):
        data = np.array(self.dataset, copy = True) 
        for pair in self.miss_pos:
            i, j = pair
            idx  = self.miss_pos.index([i, j])
            if (self.s_flag == 'discrete' and self.f_flag == False):
                data[i, j] = self.judgement_discretize(sol[:-1][self.var_pos[idx]], i, j)
                data[j, i] =  1/data[i, j]
            elif (self.s_flag != 'discrete' and self.f_flag == False):
                idx        = self.miss_pos.index([i, j])
                data[i, j] = sol[:-1][self.var_pos[idx]]
                data[j, i] =  1/data[i, j]
            elif (self.f_flag == True):
                data[i, j] = self.judgement_discretize(sol[:-1][self.var_pos[idx]], i, j)
                a, b, c    = data[i, j]
                data[j, i] = (1/c, 1/b, 1/a) 
        return data

    ################################################################################
    
    # Functions: Objective Function 0 - Consistency Ratio (MI) 
    def f0(self, variables):
        data = np.array(self.dataset, copy = True)  
        for pair in self.miss_pos:
            i, j = pair
            idx  = self.miss_pos.index([i, j])
            if (self.s_flag == 'discrete' and self.f_flag == False):
                data[i, j] = self.judgement_discretize(variables[self.var_pos[idx]], i, j)
                data[j, i] = 1/data[i, j]
            elif (self.s_flag != 'discrete' and self.f_flag == False):
                idx        = self.miss_pos.index([i, j])
                data[i, j] = variables[self.var_pos[idx]]
                data[j, i] = 1/data[i, j]
            elif (self.f_flag == True):
                data[i, j] = self.judgement_discretize(variables[self.var_pos[idx]], i, j)
                a, b, c    = data[i, j]
                data[j, i] = (1/c, 1/b, 1/a)        
        if (self.f_flag == False):
          _, adj_rc = ahp_method(data, self.wd)
        else:
          _, adj_rc = fuzzy_ahp_method(data)
        return adj_rc 
    
    # Functions: Objective Function 1 - Kendall Tau (KT)
    def f1(self, variables):
        data = np.array(self.dataset, copy = True)  
        for pair in self.miss_pos:
            i, j = pair
            idx  = self.miss_pos.index([i, j])
            if (self.s_flag == 'discrete' and self.f_flag == False):
                data[i, j] = self.judgement_discretize(variables[self.var_pos[idx]], i, j)
                data[j, i] = 1/data[i, j]
            elif (self.s_flag != 'discrete' and self.f_flag == False):
                idx        = self.miss_pos.index([i, j])
                data[i, j] = variables[self.var_pos[idx]]
                data[j, i] = 1/data[i, j]
            elif (self.f_flag == True):
                data[i, j] = self.judgement_discretize(variables[self.var_pos[idx]], i, j)
                a, b, c    = data[i, j]
                data[j, i] = (1/c, 1/b, 1/a)       
        if (self.f_flag == False):
          w1, _ = ahp_method(data, self.wd)
        else:
          w1, _ = fuzzy_ahp_method(data)
        w1             = self.rank_descending(w1)
        w2             = []
        for i in range(0, w1.shape[0]):
            if (i <= len(self.order)):
                if (self.order[i] != -1):
                    w2.append(self.order[i])
                else:
                    w2.append(w1[i])
            else:
                w2.append(w1[i])
        w2             = np.array(w2)
        kendall_tau, _ = stats.kendalltau(w1, w2)
        if (math.isnan(kendall_tau)):
            kendall_tau = -1
        return -kendall_tau

    # Function: Objective Function 1 - Kendall Tau (KT) Solutions from Brute Force  
    def f1_bf(self, idx, custom_rank):
        if (self.f_flag == False):
          w1, _ = ahp_method(self.bf_solutions[idx], self.wd)
        else:
          w1, _ = fuzzy_ahp_method(self.bf_solutions[idx])
        w1             = self.rank_descending(w1)
        w2             = custom_rank
        for i in range(0, w1.shape[0]):
            if (i <= len(w2)):
                if (w2[i] == -1):
                    w2.append(w1[i])
            else:
                w2.append(w1[i])
        w2             = np.array(w2)
        kendall_tau, _ = stats.kendalltau(w1, w2)
        if (math.isnan(kendall_tau)):
            kendall_tau = -1
        return kendall_tau    

    # Function: MultiObjective Function - Consistency Ratio (MI) & Kendall Tau (KT)
    def f0_f1(self, variables):
        result = self.alpha*self.f0(variables)/1 + (1 - self.alpha)*self.f1(variables)/10
        return result
        
    ################################################################################

    # Function: GA
    def run(self):
        self.solution = genetic_algorithm(target_function = self.lof[0],
                                          population_size = self.size,
                                          min_values      = self.minv, 
                                          max_values      = self.maxv, 
                                          generations     = self.gen,
                                          elite           = self.elt, 
                                          mutation_rate   = self.m_r,
                                          mu              = self.mu, 
                                          eta             = self.eta,
                                          verbose         = self.vbs
                                          )
        return self.solution

    ################################################################################

    # Function: Get Solutions
    def get_solution(self):
        self.indicators = []
        self.solutions  = []
        category        = 'inconsistent'
        count_con       = 0
        count_inc       = 0
        sol             = self.solution[0:int( ( (self.dataset.shape[0]**2 - self.dataset.shape[0])/2 ) + 1)]
        f0_             = self.f0(sol[:-1])
        if ('f1' in self.qidx):
          f1_ = f'{-self.f1(sol[:-1]):.4f}'
        else:
          f1_ = '-//-'
        if (f0_ > 0.1):
          category  = 'inconsistent'
          count_inc = count_inc + 1
        else:
          category  = 'consistent'
          count_con = count_con + 1
        self.indicators.append((f'{f0_:.4f}', f1_, category))
        self.solutions.append(sol)
        self.indicators = pd.DataFrame(self.indicators, columns = ['f0(MI)', 'f1(KT)', 'Consistency'])
        self.indicators = self.indicators.drop(columns = self.indicators.columns[(self.indicators == '-//-').any()])
        print(self.indicators)
        return self.indicators
    
    # Function: Get Solutions
    def get_solutions(self):
        self.indicators = []
        self.solutions  = []
        category        = 'inconsistent'
        count_con       = 0
        count_inc       = 0
        for idx in range(0, len(self.solution)):
          sol = self.solution[idx][0:int( ( (self.dataset.shape[0]**2 - self.dataset.shape[0])/2 ) + 1)]
          f0_ = self.f0(sol[:-1])
          f1_ = f'{-self.f1(sol[:-1]):.4f}'
          if (f0_ > 0.1):
            category  = 'inconsistent'
            count_inc = count_inc + 1
          else:
            category  = 'consistent'
            count_con = count_con + 1
          self.indicators.append((f'{f0_:.4f}', f1_, category))
          self.solutions.append(sol)
        self.indicators = pd.DataFrame(self.indicators, columns = ['f0(MI)', 'f1(KT)', 'Consistency'])
        self.ind_con    = self.indicators[self.indicators['Consistency'] =='consistent']
        self.ind_incon  = self.indicators[self.indicators['Consistency'] =='inconsistent']
        if (len(self.ind_con) > 0):
          self.ind_con = self.ind_con.drop_duplicates()
          self.ind_con = self.ind_con.drop(columns = self.ind_con.columns[(self.ind_con == '-//-').any()])
        if (len(self.ind_incon) > 0):
          self.ind_incon = self.ind_incon.drop_duplicates()
          self.ind_incon = self.ind_incon.drop(columns = self.ind_incon.columns[(self.ind_incon == '-//-').any()])
        self.indicators = self.indicators.drop(columns = self.indicators.columns[(self.indicators == '-//-').any()])
        print('Total Number of Inconsistent Solutions: ', count_inc)
        print('Total Number of Consistent Solutions: '  , count_con)
        print('Total Number of Unique Consistent Solutions: ', self.ind_con.shape[0])
        print(self.indicators)
        return self.indicators

    # Function: Get Solutions from Brute Force
    def get_solutions_bf(self, custom_rank = []):
        def pareto_front_points(pts, pf_min = True):
            def pareto_front(pts, pf_min):
                pf = np.zeros(pts.shape[0], dtype = np.bool_)
                for i in range(0, pts.shape[0]):
                    cost = pts[i, :]
                    if (pf_min == True):
                        g_cost = np.logical_not(np.any(pts > cost, axis = 1))
                        b_cost = np.any(pts < cost, axis = 1)
                    else:
                        g_cost = np.logical_not(np.any(pts < cost, axis = 1))
                        b_cost = np.any(pts > cost, axis = 1)
                    dominated = np.logical_and(g_cost, b_cost)
                    if  (np.any(pf) == True):
                        if (np.any(np.all(pts[pf] == cost, axis = 1)) == True):
                            continue
                    if not (np.any(dominated[:i]) == True or np.any(dominated[i + 1 :]) == True):
                        pf[i] = True
                return pf
            idx     = np.argsort(((pts - pts.mean(axis = 0))/(pts.std(axis = 0) + 1e-7)).sum(axis = 1))
            pts     = pts[idx]
            pf      = pareto_front(pts, pf_min)
            pf[idx] = pf.copy()
            return pf
        self.indicators = []
        category        = 'inconsistent'
        count_con       = 0
        count_inc       = 0
        for idx in range(0, len(self.report_df)):
          f0_ = self.report_df.iloc[idx,-1]
          f1_ = self.f1_bf(idx, custom_rank)
          if (f0_ > 0.1):
            category  = 'inconsistent'
            count_inc = count_inc + 1
          else:
            category  = 'consistent'
            count_con = count_con + 1
          self.indicators.append((f'{f0_:.4f}', f1_, category))
        self.indicators = pd.DataFrame(self.indicators, columns = ['f0(MI)', 'f1(KT)', 'Consistency'])
        self.ind_con    = self.indicators[self.indicators['Consistency'] =='consistent']
        self.ind_incon  = self.indicators[self.indicators['Consistency'] =='inconsistent']
        if (len(self.ind_con) > 0):
          self.ind_con = self.ind_con.drop(columns = self.ind_con.columns[(self.ind_con == '-//-').any()])
        if (len(self.ind_incon) > 0):
          self.ind_incon = self.ind_incon.drop(columns = self.ind_incon.columns[(self.ind_incon == '-//-').any()])
        self.indicators = self.indicators.drop(columns = self.indicators.columns[(self.indicators == '-//-').any()])
        pts             = self.indicators.copy(deep = True)
        pts['f0(MI)']   = pts['f0(MI)'].astype(float)
        pts['f1(KT)']   = pts['f1(KT)'].astype(float)
        pts             = pts.drop_duplicates()
        pf              = pareto_front_points(pts.iloc[:,:-1].values, pf_min = True)
        self.front      = pts.iloc[pf, :]
        return 

    ################################################################################
    
    # Function: Check PCM
    def check_adj_pcm(self, plot_size = 10, idx = 0):
        self.dataset_adj = self.convert_solution(self.solutions[idx])
        self.plot_environment(self.dataset_adj, plot_size = plot_size)
        if (self.f_flag == False):
          self.weights_adj, self.rc_adj = ahp_method(self.dataset_adj, wd = self.wd)
        elif(self.f_flag == True):
          self.weights_adj, self.rc_adj = fuzzy_ahp_method(self.dataset_adj)
        rank = self.rank_descending(self.weights_adj)
        for i in range(0, self.weights_adj.shape[0]):
          print('w(g'+str(i+1)+'): ', f'{self.weights_adj[i]:.3f}', '; rank(g'+str(i+1)+'): ', rank[i])
        print('')
        if (self.rc_adj > 0.100000):
          print('RC: ' + f'{self.rc_adj:.4f}', ' The solution is inconsistent, the pairwise comparisons must be reviewed')
        else:
          print('RC: ' + f'{self.rc_adj:.4f}', ' The solution is consistent')
        return self.weights_adj, rank, self.rc_adj

    ################################################################################
