#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: alexa.baltsen
"""
import sys
import time
import mysql.connector 
from mysql.connector import errorcode
from collections import OrderedDict
from anytree import Node, RenderTree
from anytree.exporter import DotExporter
from anytree.dotexport import RenderTreeGraph


superfam_db = {'user':'MU_staging',
     'password':'HjUVhg7B4kC0t2GjuScydf0fMC04Ah',
     'host':'staging-superfam.cluster-ro-curlgfcoklkj.us-west-2.rds.amazonaws.com',
     'database':'superfam'
     }   #superfam connection


def test_connection(**self):
    try:
        cnx = mysql.connector.connect(**self)
    
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Issue with connection credentials\n")
            
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist\n")
        
        else:
            print(err)
        
    else:
        cnx.close()
        return True


if test_connection(**superfam_db) != True:
    sys.exit("Connection with mysql server not established")


cnx = mysql.connector.connect(**superfam_db)

        
#replace tuple in list to list    
def tuple_to_list(query):
    q = []
    for i in query:
        t = list(i)
        q.append(t)  
    return q


def query_where(obj, tbl, inst):
    '''
    Queries database given connection info and each query field
    Returns a tuple 
    '''
    cursor = cnx.cursor()

    query = "SELECT {} FROM {} WHERE {}".format(obj, tbl, inst)
          
    cursor.execute(query)
    
    return tuple_to_list(cursor.fetchall())
            
    cursor.close()
    cnx.close() 
        
    
def query_select_column(col, tbl):
    '''
    Queries database and returns entire table
    '''
    cursor = cnx.cursor()
    
    query = "SELECT %s FROM %s" % (col, tbl)
    
    cursor.execute(query)
    
    return tuple_to_list(cursor.fetchall())

    cursor.close()
    cnx.close() 


# Splits classification columns into each level by counting '.'s 
def levelfy(col, tbl, inst, i):
    #if o 
    return query_where(col, tbl, '(SELECT LENGTH({}) - LENGTH(REPLACE({}, ".", "")))={} and not({}="")'.format(inst, inst, i, col))    
    

# Handles case where no result is found from a mysql query
# Returns modified result within nested list or empty list
def querify(query):
    try:
        r = query[0][0]
    except IndexError:
        r = query
    return r


def queritify(query):
    try:
        r = query[0]
    except IndexError:
        r = [query, query]
    return r


# fetches 6 scop level ids 
def get_scop_levels(col, lvl_id):
    return querify(query_where('description', '{}'.format(col), 'id="{}"'.format(lvl_id)))


# returns a modified list of only specific index from list of lists
def get_list_items(lst, index):
    new_list = []
    for entry in lst:
        new_list.append(entry[index])
    return new_list


# returns classification information for scop levels given domain ident
def get_scop_classification(pdb, chain, s, e):
    if chain == '':
        t = query_where('ident, description', 'cla_des_align', 'pdb="{}"'.format(pdb))
        level_IDs = get_list_items(t, 0)
        c = get_list_items(t, 1)
    elif s != '' and e != '':
        t = query_where('ident, description', 'cla_des_align', 'pdb="{}" and (SELECT(MIN(substr(pdb_region,3)) <= {})) and (SELECT(MAX(substr(pdb_region,3)) >= {}))'.format(pdb,e,s))
        level_IDs = get_list_items(t, 0)
        c = get_list_items(t, 1)
    else:
        t = query_where('ident, description', 'cla_des_align', 'pdb="{}" and pdb_region="{}%"'.format(pdb, chain))
        level_IDs = get_list_items(t, 0)
        c = get_list_items(t, 1)
        
    return c, level_IDs


# returns pfam clan and child info
def get_pfam_classification(pdb, s, e):
    if s == '' and e == '':    
        pfam = querify(query_where('pfam_id', 'pdb_map','pdb_id="{}"'.format(pdb)))
        parent = querify(query_where('parent','PFAM_hie','child="{}"'.format(pfam)))
        pfam_desc = querify(query_where('description','PFAM_info','level="pfam" and acc="{}"'.format(pfam)))
        parent_desc = querify(query_where('description','PFAM_info','level="clan" and acc="{}"'.format(parent)))
    else:
        pfam = querify(query_where('pfam_id', 'pdb_map','pdb_id="{}" and (SELECT(MIN(position) <= {})) and (SELECT(MAX(position) >= {}))'.format(pdb,e,s)))
        parent = querify(query_where('parent','PFAM_hie','child="{}"'.format(pfam)))
        pfam_desc = querify(query_where('description','PFAM_info','level="pfam" and acc="{}"'.format(pfam)))
        parent_desc = querify(query_where('description','PFAM_info','level="clan" and acc="{}"'.format(parent)))
        
    return parent_desc, pfam_desc


# returns ec and go ids and entry information given a superfam or fam id 
def get_go_ec_classification(lvl):
    go_id = querify(query_where('go', 'GO_mapping', 'id="{}"'.format(lvl)))
    ec_id = querify(query_where('ec', 'EC_mapping', 'id={}'.format(lvl)))

    if go_id and ec_id:
        go_id = str(go_id).zfill(7)
        go_process, go_function = queritify(query_where('namespace, name', 'GO_info', 'go="{}"'.format(go_id)))
        ec_process, ec_function = queritify(query_where('namespace, description', 'EC_info', 'ec="{}"'.format(ec_id)))
    else:
        go_id, ec_id, go_process, go_function, ec_process, ec_function = '','','','','',''
        
    return go_id, ec_id, go_process, go_function, ec_process, ec_function


# Traverses through GO table
# Gets all children of parent GO id and adds to tree    
def go_getter(node, go_id, p):
    node_name = querify(query_where('name', 'GO_info', 'go="{}"'.format(go_id)))
    node = Node(node_name, parent=p)
    children = query_where('child', 'GO_hie', 'parent="{}" and distance="1"'.format(go_id))
    
    for kid in children:
        child = str(kid).strip('[]').zfill(7)
        kids = child
        go_getter(kids, kids, node)
    
    return node


# Queries cath solid database for 4 levels of sequence clusters
# Is modified to return a list of similar cath ids for each cluster
# also returns cluster lists of benchmark similarities between scop and cath for flagging later
def get_cath_solid(cath):
    c, a, t, h, s, s35, s60, s95, s100 = cath
    sf = query_where('domain_identifier, benchmark_set', 'cath_domain_description', 'class="{}" and architecture="{}" and topology="{}" and superfamily="{}" and s35_cluster="{}"'.format(c, a, t, h, s35))
    osf = query_where('domain_identifier, benchmark_set', 'cath_domain_description', 'class="{}" and architecture="{}" and topology="{}" and superfamily="{}" and s35_cluster="{}" and s60_cluster="{}"'.format(c, a, t, h, s35, s60))
    lsf = query_where('domain_identifier, benchmark_set', 'cath_domain_description', 'class="{}" and architecture="{}" and topology="{}" and superfamily="{}" and s35_cluster="{}" and s60_cluster="{}" and s95_cluster="{}"'.format(c, a, t, h, s35, s60, s95))
    i = query_where('domain_identifier, benchmark_set', 'cath_domain_description', 'class="{}" and architecture="{}" and topology="{}" and superfamily="{}" and s35_cluster="{}" and s60_cluster="{}" and s95_cluster="{}" and s100_cluster="{}"'.format(c,a,t,h,s35,s60,s95,s100))          
    
    return sf, osf, lsf, i, get_list_items(sf,0), get_list_items(osf,0),get_list_items(lsf,0),get_list_items(i,0)


def sf_or_f_go(sf_go_id, f_go_id, par):  
    if sf_go_id == f_go_id: 
        sf_f_go_nodes = go_getter('root', sf_go_id, None)
        sf_f_go_nodes.parent = par
    elif sf_go_id != '':
        sf_go_nodes = go_getter('root', sf_go_id, None)
        sf_go_nodes.parent = par
    elif f_go_id != '':
        f_go_nodes = go_getter('root', f_go_id, None)
        f_go_nodes.parent = par
    
    else:
        return 0


def get_death_domain(uniprot):
    if uniprot != '':
        d = query_where('short_name, process_id, pathway_family', 'death_domains', 'uniprot_id={}'.format(uniprot))
        
    else:
        d = ['','','']
    return d
    

def pipe_input(p,c,u,s,e):
    if p != '' and c == '' and u == '' and s == '' and e == '':
        domain_chain = querify(query_where('pdb_region', 'cla_des_align', 'pdb="{}" LIMIT 1'.format(p)))
        cath_id = '1HDJ{}%'.format(domain_chain[0])
        return cath_id
    elif p != '' and c != '' and u == '' and s == '' and e == '':
        cath_id = '{}{}%'.format(p,c)
        return cath_id
    elif p != '' and u != '' and s == '' and e == '':
        pdb, chain = queritify(query_where('pdb_id, domain_chain','uniprot_pdb_pos_overlap', 'uniprot_id="{}" and description="{} {}%" LIMIT 1'.format(u,p,c)))
        cath_id = '{}{}%'.format(pdb,chain)
        return cath_id
    elif p != '' and u == '' and s != '' and e != '':
        #pdb, chain = query_where('pdb_id, domain_chain','uniprot_pdb_pos_overlap','pdb_id="{}" and ({}>=LEAST(pdb_beg_pos, uniprot_start_pos) and {}<=GREATEST(pdb_end_pos, uniprot_stop_pos)) LIMIT 1'.format(p,e,s))[0][0]
        pdb, chain = queritify(query_where('pdb_id, domain_chain','uniprot_pdb_pos_overlap','pdb_id="{}" and ({}>=LEAST(pdb_beg_pos, uniprot_start_pos) and {}<=GREATEST(pdb_end_pos, uniprot_stop_pos)) LIMIT 1'.format(p,e,s)))
        cath_id = '{}{}%'.format(pdb,chain)
        if pdb == '' and chain == '':
            pdb, chain = queritify(query_where('pdb_id, domain_chain','uniprot_pdb_pos_overlap','pdb_id="{}" and ({}>=LEAST(pdb_beg_pos, uniprot_start_pos) and {}<=GREATEST(pdb_end_pos, uniprot_stop_pos)) LIMIT 1'.format(p,e,s)))
        return cath_id
    elif p != '' and u != '' and s != '' and e != '':
        pdb, chain = queritify(query_where('pdb_id, domain_chain','uniprot_pdb_pos_overlap','(pdb_id="{}" and uniprot="{}") and ({}>=LEAST(pdb_beg_pos, uniprot_start_pos) and {}<=GREATEST(pdb_end_pos, uniprot_stop_pos)) LIMIT 1'.format(p,u,e,s)))
        cath_id = '{}{}%'.format(pdb,chain)
        if pdb == '' and chain == '':
            pdb, chain = queritify(query_where('pdb_id, domain_chain','uniprot_pdb_pos_overlap','pdb_id="{}" and ({}>=LEAST(pdb_beg_pos, uniprot_start_pos) and {}<=GREATEST(pdb_end_pos, uniprot_stop_pos)) LIMIT 1'.format(p,e,s)))
        return cath_id
    

# Tests input for correct identifiers 
def input_test(pdb, chain, uniprot, start, stop):
    if pdb== '' and chain == '' and uniprot== '' and start== '' and stop == '':
        print('No input has been given')
    elif len(pdb) != 4 and pdb != '':
        print('Invalid PDB ID')
    elif chain != '' and str(chain).isalpha() == False:
        print('Invalid domain chain')
    elif len(uniprot) != 6 and uniprot != '':
        print('Invalid Uniprot ID')
    elif (str(start).isdigit() == False or str(stop).isdigit() == False) and (start != '' or stop != ''):
        print('Invalid start or end position')
    else:
        return 0


bins = OrderedDict()

# Set up of classification levels 

bins['Class']=levelfy('description', 'cath_names', 'classification', 0) 
bins['SCOP Class'] = query_select_column('DISTINCT(description)', 'cl')

bins['Architecture']=levelfy('description', 'cath_names', 'classification', 1)
bins['SCOP Fold'] = query_select_column('DISTINCT(description)', 'fo')

bins['Topology']=levelfy('description', 'cath_names', 'classification', 2) 

bins['SCOP Homologous Superfamily'] = query_select_column('DISTINCT(description)', 'su')
bins['Homologous Superfamily']=levelfy('description', 'cath_names', 'classification', 3)

bins['Pfam Clan'] = query_where('description','PFAM_info', 'level="clan"')
bins['SCOP Family'] = query_select_column('DISTINCT(description)', 'fa')
bins['Pfam Family'] = query_where('description','PFAM_info', 'level="pfam"')
bins['SCOP Domain'] = query_select_column('DISTINCT(description)', 'do')

bins['Process'] = query_select_column('DISTINCT(namespace)', 'GO_info')
bins['Function'] = query_select_column('DISTINCT(name)', 'GO_info')


bins['EC Root']=levelfy('description','EC_info','ec',0)
bins['EC 1']=levelfy('description','EC_info','ec',1)
bins['EC 2']=levelfy('description','EC_info','ec',2)
bins['EC 3']=levelfy('description','EC_info','ec',3)

'''
# Tree display, try to map entries to each other if possible
bins_root = Node(bins['Class'])
scop_class = Node(bins['SCOP Class'], parent=bins_root) 
architecture = Node(bins['Architecture'], parent=scop_class) 
scop_fold = Node(bins['SCOP Fold'], parent=architecture)
topology = Node(bins['Topology'], parent=scop_fold)
scop_supfam = Node(bins['SCOP Homologous Superfamily'], parent=topology)
superfam = Node(bins['Homologous Superfamily'], parent=scop_supfam)
clan = Node(bins['Pfam Clan'], parent=superfam)
scop_fam = Node(bins['SCOP Family'], parent=clan)
pfam_fam = Node(bins['Pfam Family'], parent=scop_fam)
scop_domain = Node(bins['SCOP Domain'], parent=pfam_fam)
process = Node(bins['Process'], parent=scop_domain)
function = Node(bins['Function'], parent=process)
ec_root = Node(bins['EC Root'], parent=function)
ec_1 = Node(bins['EC 1'], parent=ec_root)
ec_2 = Node(bins['EC 2'], parent=ec_1)
ec_3 = Node(bins['EC 3'], parent=ec_2)

DotExporter(bins_root).to_dotfile("tree.dot")
#RenderTreeGraph(bins_root).to_picture("~/Documents/domain_classification_tree.png")
'''

#hits = OrderedDict()
hits = {}

# Query through databases to classify domain based on pdbid
#test = '1oip A:25-90,1t89 C:87-171,1pls A:,1E8Z,2WR0'

pdb = ['1HDJ','1oip','1t89','1pls','1E8Z','2WR0','barf','foo']
chain = ['','A','C','A','','T', '',8]
uniprot = ['','', '', '','','','yh8943','']
start_pos = ['',25,200,'','',888,9,0 ]
end_pos = ['',90,300,'','',9,-8,'i']
#seq = []



for p, c, u, s, e in zip(pdb, chain, uniprot, start_pos, end_pos):
    if input_test(p, c, u, s, e) == 0:
        
        # CATH classification section 
        try:              
            cath = query_where('*', 'cath_domain_description', 'domain_identifier LIKE "{}"'.format(pipe_input(p,c,u,s,e)))[0]
        except IndexError:
            print('No classification found for {} {} {} {} {}'.format(p,c,u,s,e))
            break
        #need to fix, break out of if statement and go to next test group
            
        seq_fam_benchmark, orth_seq_fam_benchmark, like_seq_fam_benchmark, ident_benchmark, hits['Sequence Family'],hits['Orthogous Sequence Family'], hits['"Like" Sequence Family'], hits['Identical'] = get_cath_solid(cath[1:10])
        
        hits['Class'], hits['Architecture'], hits['Topology'], hits['Homologous Superfamily'], cath_scop_supfam_mapping = cath[1:6] 
        
        
        # SCOP classification section
        scop_hits, level_IDs = get_scop_classification(p, c, s , e)
        
    
        # Pfam 
        hits['Pfam Clan'], hits['Pfam Family'] = get_pfam_classification(p,s,e)
            
        
        # Superfamily GO & EC
        sf_go_id, sf_ec_id, hits['GO Process'], hits['GO Function'], hits['EC Process'], hits['EC Function'] = get_go_ec_classification(level_IDs[3])
        level_IDs.append(sf_go_id) 
        level_IDs.append(sf_ec_id)
        
        
        # Family GO & EC
        f_go_id, f_ec_id, hits['GO Process'], hits['GO Function'], hits['EC Process'], hits['EC Function'] = get_go_ec_classification(level_IDs[4])
        level_IDs.append(f_go_id) 
        level_IDs.append(f_ec_id)
        
        
        # Death domain 
        
        
        death_domain, process, pathway = get_death_domain(u)
     
        
        print('Creating tree\n')
        
        
        hits_root = Node(p)
        
        hits_class = Node(hits['Class'], parent=hits_root)
        
        hits_scop_class = Node(scop_hits[1], parent=hits_class)
        
        hits_arch = Node(hits['Architecture'], parent=hits_scop_class)
        
        hits_fold = Node(scop_hits[2], parent=hits_arch)
        
        hits_top = Node(hits['Topology'], parent=hits_fold)
        
        hits_scop_supfam = Node(scop_hits[3], parent=hits_top)
        
        hits_superfam = Node(hits['Homologous Superfamily'], parent=hits_scop_supfam)
        hits_seq_fam_name = Node('35%S Sequence Family', parent=hits_superfam)
        
        hits_seq_fam = Node(hits['Sequence Family'], parent=hits_seq_fam_name)
        
        hits_orth_seq_fam_name = Node('60%S Orthogous Sequence Family', parent=hits_seq_fam_name)
        hits_orth_seq_fam = Node(hits['Orthogous Sequence Family'], parent=hits_orth_seq_fam_name)
        
        hits_like_seq_fam_name = Node('95%S "Like" Sequence Family', parent=hits_orth_seq_fam_name)
        hits_like_seq_fam = Node(hits['"Like" Sequence Family'], parent=hits_like_seq_fam_name)
        
        hits_ident_name = Node('100%S Identical', parent=hits_like_seq_fam_name)
        hits_ident = Node(hits['Identical'], parent=hits_ident_name)
        
        hits_clan = Node(hits['Pfam Clan'], parent=hits_superfam)
        
        hits_scop_fam = Node(scop_hits[4], parent=hits_clan)
        
        hits_pfam_fam = Node(hits['Pfam Family'], parent=hits_scop_fam)
        
        hits_scop_domain = Node(scop_hits[5], parent=hits_pfam_fam)
        
        hits_go_process = Node(hits['GO Process'], parent=hits_scop_domain)
        
        hits_go_function = Node(hits['GO Function'], parent=hits_go_process)
        
        '''
        print('sf or f go working')
        sf_or_f_go(sf_go_id, f_go_id, hits_go_function) #adds children of biological process
        print('done with that')
        '''
        
        hits_ec_process = Node(hits['EC Process'], parent=hits_scop_domain)
        
        hits_ec_function = Node(hits['EC Function'], parent=hits_ec_process)
        
        hits_death_domain = Node(death_domain, parent=hits_scop_domain)
        hits_death_process = Node(process, parent=hits_death_domain)
        hits_death_pathway = Node(pathway, parent=hits_death_domain)
        
        for pre, fill, node in RenderTree(hits_root):
            print('%s%s' % (pre, node.name))
            
        print('\n')
            
    else:
        #sys.exit('Incorrect input to run domain classifier')
        print('\n{}:Incorrect input to run domain classifier'.format(p))

