import kahypar
import copy

def graphPartitioningKahypar(adj_hash,parts_num,e=0.03):
  num_nodes = len(adj_hash)
  num_nets  = 0
  
  if(parts_num > num_nodes): return -1,list(range(num_nodes))

  hyperedge_indices = [0]
  hyperedges = []
  
  node_weights = [1]*num_nodes
  edge_weights = []

  for s_no,l_no_map in adj_hash.items():
    for l_no,ew in l_no_map.items():
      if(s_no > l_no): continue
      num_nets += 1
      edge_weights.append(ew)
      hyperedges.append(s_no)
      hyperedges.append(l_no)
      hyperedge_indices.append(num_nets*2)

  k=parts_num
  
  if(not hyperedges):
    num_nets = 1
    hyperedges = list(range(num_nodes))
    hyperedge_indices = [0,num_nodes]
    
  hypergraph = kahypar.Hypergraph(num_nodes, num_nets, hyperedge_indices, hyperedges, k, edge_weights, node_weights)

  context = kahypar.Context()
  context.loadINIconfiguration("/home/yjt/gaoHengyue/kahypar-master/config/cut_kKaHyPar_sea20.ini")

  context.setK(k)

  context.setEpsilon(e)

  context.suppressOutput(True)

  ret = kahypar.partition(hypergraph, context)

  mapping_arr = []
  for i in range(num_nodes): mapping_arr.append(hypergraph.blockID(i))

  return -1,mapping_arr

def graphPartitioning(adj_hash,parts_num,cut_way = "",nodes_res_map = None,key_res_no = None):
    recursive_flag = True if "r" in cut_way else False if "k" in cut_way else None
    if("w" in cut_way):
        adj,xadj,w,vw = transformAdjHashTo4List(adj_hash,nodes_res_map,key_res_no)
        edgecuts,parts = pymetis.part_graph(nparts=parts_num,adjncy=adj,xadj=xadj,eweights=w,vweights = vw,recursive=recursive_flag)
    else:
        adj,xadj,w = transformAdjHashTo3List(adj_hash)
        edgecuts,parts = pymetis.part_graph(nparts=parts_num,adjncy=adj,xadj=xadj,eweights=w,recursive=recursive_flag)

    return edgecuts,parts

def transformAdjHashTo3List(adj_hash):
    adj = []
    xadj = [0]
    w = []

    xadj_head = 0
    for node_no in range(len(adj_hash)):
        adj_counts = len(adj_hash[node_no])
        for adj_node_no in adj_hash[node_no].keys():
            adj.append(adj_node_no)
            w.append(adj_hash[node_no][adj_node_no])
        xadj.append(xadj_head + adj_counts)
        xadj_head += adj_counts
    
    return adj,xadj,w

def transformAdjHashTo4List(adj_hash,nodes_res_map,key_res_no):
    adj = []
    xadj = [0]
    w = []
    vw = []

    xadj_head = 0
    for node_no in range(len(adj_hash)):
        adj_counts = len(adj_hash[node_no])
        for adj_node_no in adj_hash[node_no].keys():
            adj.append(adj_node_no)
            w.append(adj_hash[node_no][adj_node_no])
        xadj.append(xadj_head + adj_counts)
        xadj_head += adj_counts
        vw.append(nodes_res_map[node_no][key_res_no])
    return adj,xadj,w,vw

def multiLevelRefine(adj_hash,fpga_res_list,nodes_res_map,metis_mapping_arr):
    ori_nodes_num = len(metis_mapping_arr)

    cluster_res_map = {}
    for node_no,cluster_no in enumerate(metis_mapping_arr):
        if(cluster_no not in cluster_res_map): cluster_res_map[cluster_no] = [0]*9
        for res_no in range(9): cluster_res_map[cluster_no][res_no] += nodes_res_map[node_no][res_no]

    coarse_scale_str = ""
    cut_size_str = ""
    max_cut_size_str = ""

    max_cut_size_trace = []
    clusters_adj_hash = getClustersAdjHash(metis_mapping_arr,adj_hash)
    cut_count = 0
    cut_max = -1
    for cluster_no,adj_clusters_map in clusters_adj_hash.items():
        for adj_cluster,edge_w in adj_clusters_map.items():
            cut_count += edge_w
            if(edge_w > cut_max): cut_max = edge_w 
    cut_count /= 2
    cut_size_str += str(cut_count) + '	'
    max_cut_size_str += str(cut_max) + '	'
    coarse_scale_str += "0" + '	'
    max_cut_size_trace.append(cut_max)
    best_cut_max = cut_max
    best_mapping_arr = copy.deepcopy(metis_mapping_arr)

    last_cut_max = cut_max
    last_mapping_arr = copy.deepcopy(metis_mapping_arr)
    last_cluster_res_map = copy.deepcopy(cluster_res_map)
    last_loop_cut = -1
    while(True):
        coarse_scale_i = 0

        window_len = 10
        window_p = 0
        window_arr = [-1]*window_len

        break_flag = True
        coarsen_info_map = {}
        while(True):
            exp_fig = 1.1
            coarse_scale = int((((exp_fig**coarse_scale_i) - 1)*80)/(exp_fig**20)) + 1

            for i in range(len(metis_mapping_arr)-1,0,-1):
                if(metis_mapping_arr[i] == metis_mapping_arr[0]): continue
                break_flag = False
                break
            if(break_flag): break

            coarsened_mapping_arr,super_res_arr,super_adj_hash,super2oriNodes_mapping = graphCoarsening_old(ori_nodes_num,metis_mapping_arr,coarse_scale,adj_hash,nodes_res_map,coarsen_info_map)
            
            superNodesTransposition(coarsened_mapping_arr,super_res_arr,super_adj_hash,fpga_res_list,cluster_res_map,coarsen_info_map)

            metis_mapping_arr = graphUncoarsening(coarsened_mapping_arr,super2oriNodes_mapping,metis_mapping_arr)
            
            clusters_adj_hash = getClustersAdjHash(metis_mapping_arr,adj_hash)
            cut_count = 0
            cut_max = -1
            for cluster_no,adj_clusters_map in clusters_adj_hash.items():
                for adj_cluster,edge_w in adj_clusters_map.items():
                    cut_count += edge_w
                    if(edge_w > cut_max): cut_max = edge_w 
            cut_count /= 2
            cut_size_str += str(cut_count) + '	'
            max_cut_size_str += str(cut_max) + '	'
            coarse_scale_str += str(coarse_scale) + '	'
            max_cut_size_trace.append(cut_max)
            if(cut_max < best_cut_max):
                best_cut_max = cut_max
                best_mapping_arr = copy.deepcopy(metis_mapping_arr)

            window_arr[window_p%window_len] = last_cut_max
            window_p += 1
            coarse_scale_i += 1

            cut_ave = sum(window_arr) / window_len
            cut_var = sum((cut_-cut_ave)**2 for cut_ in window_arr) / window_len
            if(cut_var < 1): break

        if(break_flag): break
        if(best_cut_max/last_loop_cut > 0.98): break
        last_loop_cut = best_cut_max

    return best_mapping_arr

def multiLevelRefinePartitioning(adj_hash,parts_num,fpga_res_list,nodes_res_map):
    cut_way_str = "r"
    cut_way_str = ""
    edge_cuts_num,metis_mapping_arr = graphPartitioning(adj_hash,parts_num,cut_way_str,nodes_res_map,0)
    return multiLevelRefine(adj_hash,fpga_res_list,nodes_res_map,metis_mapping_arr)

def graphUncoarsening(coarsened_mapping_arr,super2oriNodes_mapping,mapping_arr):
    for super_node,small_node_map in super2oriNodes_mapping.items():
        for small_node in small_node_map.keys():
            mapping_arr[small_node] = coarsened_mapping_arr[super_node]
    return mapping_arr

def superNodesTransposition(coarsened_mapping_arr,super_res_arr,super_adj_hash,fpga_res_list,cluster_res_map,coarsen_info_map):
    candidates_hash = {}
    candidates_L2_hash = {}

    avoid_ring_visited_map = {}

    for super_node in range(len(coarsened_mapping_arr)):
        updateNodeInfoInCandidatesHash(super_node,coarsened_mapping_arr,super_adj_hash,candidates_hash)

    while(candidates_hash):
        success_flag, super_node = tryTasnsferOneNode(candidates_hash,super_res_arr,super_adj_hash,fpga_res_list,coarsened_mapping_arr,cluster_res_map,coarsen_info_map)
        
        if(success_flag):
            avoid_ring_visited_map[super_node] = True
            del candidates_hash[super_node]

            for adj_super_node in super_adj_hash[super_node].keys():
                updateNodeInfoInCandidatesHash(adj_super_node,coarsened_mapping_arr,super_adj_hash,candidates_hash)

            for l2_node in candidates_L2_hash.keys():
                if(l2_node not in candidates_hash): continue
                candidates_L2_hash[l2_node] = candidates_hash[l2_node]
                del candidates_hash[l2_node]

            temp_L2_hash = copy.deepcopy(candidates_L2_hash)
            while(temp_L2_hash):
                success_flag, l2_node = tryTasnsferOneNode(temp_L2_hash,super_res_arr,super_adj_hash,fpga_res_list,coarsened_mapping_arr,cluster_res_map,coarsen_info_map)
                del temp_L2_hash[l2_node]
                if(not success_flag): continue
                avoid_ring_visited_map[super_node] = True
                
                del candidates_L2_hash[l2_node]
                for adj_super_node in super_adj_hash[l2_node].keys():
                    updateNodeInfoInCandidatesHash(adj_super_node,coarsened_mapping_arr,super_adj_hash,candidates_hash)

                for l2_node in candidates_L2_hash.keys():
                    if(l2_node not in candidates_hash): continue
                    candidates_L2_hash[l2_node] = candidates_hash[l2_node]
                    if(l2_node in temp_L2_hash): temp_L2_hash[l2_node] = candidates_hash[l2_node]
                    del candidates_hash[l2_node]
        else:
            candidates_L2_hash[super_node] = candidates_hash[super_node]
            del candidates_hash[super_node]

    while(True):
        len_begin = len(candidates_L2_hash)
        temp_L2_hash = copy.deepcopy(candidates_L2_hash)
        while(temp_L2_hash):
            success_flag, l2_node = tryTasnsferOneNode(temp_L2_hash,super_res_arr,super_adj_hash,fpga_res_list,coarsened_mapping_arr,cluster_res_map,coarsen_info_map)
            del temp_L2_hash[l2_node]
            if(not success_flag): continue
            del candidates_L2_hash[l2_node]
            for adj_super_node in super_adj_hash[l2_node].keys():
                updateNodeInfoInCandidatesHash(adj_super_node,coarsened_mapping_arr,super_adj_hash,candidates_hash)

            for l2_node in candidates_L2_hash.keys():
                if(l2_node not in candidates_hash): continue
                candidates_L2_hash[l2_node] = candidates_hash[l2_node]
                if(l2_node in temp_L2_hash): temp_L2_hash[l2_node] = candidates_hash[l2_node]
                del candidates_hash[l2_node]
        if(len_begin == len(candidates_L2_hash)): break

def tryTasnsferOneNode(candidates_hash,super_res_arr,super_adj_hash,fpga_res_list,coarsened_mapping_arr,cluster_res_map,coarsen_info_map):
    max_gain = -1
    target_node = -1
    for super_node,node_info in candidates_hash.items():
        if(node_info["gain"] <= max_gain): continue
        max_gain = node_info["gain"]
        target_node = super_node

    target_clusters = candidates_hash[target_node]["target_cluster"]
    if(len(target_clusters)>1):
        clusters_degree = []
        for candidate_cluster in target_clusters:
            cluster_degree = 0
            for s_super_node,adj_info in super_adj_hash.items():
                if(coarsened_mapping_arr[s_super_node] != candidate_cluster): continue   
                for adj_node,edge_weight in adj_info.items():
                    if(coarsened_mapping_arr[adj_node] == candidate_cluster): continue   
                    cluster_degree += edge_weight
            clusters_degree.append(cluster_degree)

        for i in range(len(target_clusters)-1):
            for j in range(len(target_clusters)-1-i):
                if(clusters_degree[j] >= clusters_degree[j+1]): continue
                clusters_degree[j],target_clusters[j],clusters_degree[j+1],target_clusters[j+1] = clusters_degree[j+1],target_clusters[j+1],clusters_degree[j],target_clusters[j]  

    target_node_res = super_res_arr[target_node]
    ori_cluster = coarsened_mapping_arr[target_node]
    target_cluster_arr = copy.deepcopy(candidates_hash[target_node]["target_cluster"])
    while(target_cluster_arr):
        target_cluster = target_cluster_arr.pop()

        res_OK = True
        for res_no in range(9):
            if(cluster_res_map[target_cluster][res_no] + super_res_arr[target_node][res_no] <= fpga_res_list[res_no]): continue
            res_OK = False
            break

        if(not res_OK): continue

        for res_no in range(9):
            cluster_res_map[target_cluster][res_no]                     += super_res_arr[target_node][res_no]
            cluster_res_map[coarsened_mapping_arr[target_node]][res_no] -= super_res_arr[target_node][res_no]

        coarsened_mapping_arr[target_node] = target_cluster

        if(ori_cluster in coarsen_info_map): del coarsen_info_map[ori_cluster]      
        if(target_cluster in coarsen_info_map): del coarsen_info_map[target_cluster]   
        return True,target_node
    
    return False,target_node

def updateNodeInfoInCandidatesHash(super_node,coarsened_mapping_arr,super_adj_hash,candidates_hash):
    in_cluster = coarsened_mapping_arr[super_node]
    node_clusters_weight = {}
    for adj_super_node,weight in super_adj_hash[super_node].items():
        adj_cluster = coarsened_mapping_arr[adj_super_node]
        if(adj_cluster not in node_clusters_weight): node_clusters_weight[adj_cluster] = 0
        node_clusters_weight[adj_cluster] += weight

    max_weight = -1
    max_weight_cluster = []
    for adj_cluster,weight in node_clusters_weight.items():
        if(max_weight > weight): continue
        if(max_weight == weight):
            max_weight_cluster.append(adj_cluster)
            continue
        max_weight = weight
        max_weight_cluster = [adj_cluster]

    if(max_weight == -1):
        if(super_node in candidates_hash): del candidates_hash[super_node]                 
        return True           
    if(in_cluster not in node_clusters_weight): node_clusters_weight[in_cluster] = 0
    if(max_weight == node_clusters_weight[in_cluster]):
        if(super_node in candidates_hash): del candidates_hash[super_node]                 
        return True

    candidates_hash[super_node] = {}
    candidates_hash[super_node]["target_cluster"] = max_weight_cluster
    candidates_hash[super_node]["gain"] = max_weight

def graphCoarsening_old(ori_nodes_num,mapping_arr,coarse_scale,adj_hash,nodes_res_map,coarsen_info_map):
    super_res_arr = []                                                  
    coarsened_mapping_arr = []                                          
    super2oriNodes_mapping = {}                                         

    super_nodes_no = 0

    sub_graph_info_map = {}

    for old_node_no,cluster_no in enumerate(mapping_arr):
        if(cluster_no in coarsen_info_map):
            sub_graph_info_map[cluster_no] = {}
            continue
        if(cluster_no not in sub_graph_info_map): sub_graph_info_map[cluster_no] = [{},[],{}] 
        sub_graph_adj_hash, new2old_nodes_arr,old2new_nodes_arr = sub_graph_info_map[cluster_no]
        if(old_node_no not in old2new_nodes_arr):
            old2new_nodes_arr[old_node_no] = len(new2old_nodes_arr)
            sub_graph_adj_hash[len(new2old_nodes_arr)] = {}
            new2old_nodes_arr.append(old_node_no)

    for old_node_no,old_nodes_map in adj_hash.items():
        cluster_no = mapping_arr[old_node_no]
        if(cluster_no in coarsen_info_map): continue
        sub_graph_adj_hash, new2old_nodes_arr,old2new_nodes_arr = sub_graph_info_map[cluster_no]
        s_node = old2new_nodes_arr[old_node_no]
        for old_adj_node,edge_weight in old_nodes_map.items():
            if(mapping_arr[old_adj_node] != cluster_no): continue
            l_node = old2new_nodes_arr[old_adj_node]
            sub_graph_adj_hash[s_node][l_node] = edge_weight
            sub_graph_adj_hash[l_node][s_node] = edge_weight

    cluster_no_visited_map = {}
    for old_node_no,cluster_no in enumerate(mapping_arr):
        if(cluster_no in cluster_no_visited_map): continue
        sub_graph_info = sub_graph_info_map[cluster_no]
        cluster_no_visited_map[cluster_no] = True
        
        temp_super_nodes_mapping = {}

        if(sub_graph_info):
            sub_graph_adj_hash, new2old_nodes_arr = sub_graph_info[0],sub_graph_info[1]
            coarsen_info_map[cluster_no] = (sub_graph_adj_hash,new2old_nodes_arr)
        else:
            sub_graph_adj_hash, new2old_nodes_arr = coarsen_info_map[cluster_no]    
             
        edgecuts,metis_mapping_arr = graphPartitioningKahypar(sub_graph_adj_hash,coarse_scale)
        
        if(cluster_no == 3 and coarse_scale == 2 and False):
            old2cluster_map = {}
            for new_node,cluster_no in enumerate(metis_mapping_arr):
                old2cluster_map[new2old_nodes_arr[new_node]] = cluster_no
            old_node_arr = list(old2cluster_map.keys())
            old_node_arr.sort()

        for sub_node in range(len(metis_mapping_arr)):
            sub_super_node = metis_mapping_arr[sub_node]
            if(sub_super_node not in temp_super_nodes_mapping):
                temp_super_nodes_mapping[sub_super_node] = super_nodes_no
                super_res_arr.append([0]*9)
                coarsened_mapping_arr.append(cluster_no)
                super2oriNodes_mapping[super_nodes_no] = {}
                super_nodes_no += 1

            super_node_now = temp_super_nodes_mapping[sub_super_node]
            old_small_node = new2old_nodes_arr[sub_node]

            for res_no in range(9): super_res_arr[super_node_now][res_no] += nodes_res_map[old_small_node][res_no]

            super2oriNodes_mapping[super_node_now][old_small_node] = True

    super_mapping_arr = [0]*ori_nodes_num
    for super_node,ori_nodes_map in super2oriNodes_mapping.items():
        for ori_node in ori_nodes_map.keys():
            super_mapping_arr[ori_node] = super_node

    super_adj_hash = getClustersAdjHash(super_mapping_arr,adj_hash)     

    return coarsened_mapping_arr,super_res_arr,super_adj_hash,super2oriNodes_mapping

def getClustersAdjHash(mapping_arr,adj_hash):
    clusters_hash = {}
    for node_no,cluster_no in enumerate(mapping_arr):
        if(cluster_no not in clusters_hash): clusters_hash[cluster_no] = {}
        for adj_node,edge_weight in adj_hash[node_no].items():
            adj_cluster = mapping_arr[adj_node]
            if(adj_cluster == cluster_no): continue
            if(adj_cluster not in clusters_hash[cluster_no]): clusters_hash[cluster_no][adj_cluster] = 0
            clusters_hash[cluster_no][adj_cluster] += edge_weight
    return clusters_hash