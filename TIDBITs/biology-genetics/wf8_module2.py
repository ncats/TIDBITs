import json
import requests
import pandas
import time
from multiprocessing import Pool
import logging
from requests.exceptions import HTTPError

base_url = 'http://biggim.ncats.io/api'

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logging.debug("test")

"""Linking tissues to columns"""
def get_columns(tissues,table='BigGIM_70_v1'):
    
    columns = []
    for t in tissues:
    	columns += get_column(t,table=table)
    columns =  list(set(columns)) # get rid of redundant columns
    
    a = columns
    b = ['TCGA', 'Pvalue']
    columns = [sa for sa in a if not any(sb in sa for sb in b)]
    
    if len(columns) == 0:
        raise Exception("No Big GIM columns related to %s" % (str(tissues)))
    else:
        print("Returned %i Big GIM columns" % (len(columns), ))

    return columns

"""Running BigGIM"""
def call_biggim(genes, tissues, limit=1000000, average_columns=False, query_id2=False, return_genes=False, N=250): 
    """
    Parameters
    ----------
    genes : list of str 

    columns : list of str

    limit : int

    query_id2 : bool
        If query_id2, the genes are duplicated for the ids2 field. 
        This forces biggim to return all interactions within genes.
    """

    columns = get_columns(tissues,table='BigGIM_70_v1')
    print(columns)
    
    gene_list = ','.join(genes)

    example_query = {
        "restriction_gt": ','.join(["%s,-2.0" % (c,) for c in columns]),
        "table": "BigGIM_70_v1",
        "columns": ','.join(columns),
        "ids1": gene_list,
        "limit": limit,
        "average_columns": average_columns,
    }

    if query_id2: 
        example_query["ids2"] = gene_list

    try:
        query_submit = get('biggim/query', data=example_query)
        jprint(query_submit)
    except requests.HTTPError as e:
        print(e)
        
        jprint(e.response.json())

    error_counter = 0
    try:
        while True:
            try:
                query_status = get('biggim/status/%s' % (query_submit['request_id'],))
                jprint(query_status)

                if query_status['status'] != 'running':
                    # query has finished
                    break
                else:
                    time.sleep(5)
                    print("Checking again")
            except requests.HTTPError as e:
                print(e)
                time.sleep(5)
                error_counter += 1
                if error_counter > 3:
                    print("Giving up")
                    raise
                else:
                    print("Trying again.")
    except requests.HTTPError as e:
        print(e)

        jprint(e.response.json())
        
    result = pandas.concat(map(pandas.read_csv, query_status['request_uri']))
    result = result.iloc[:,1:]
    
    if return_genes: 
        df = result
        print(df)
        df  = df.sort_values(by='mean',ascending=False)
        df = df.reset_index()
        print(df)
        del df['index']
        gene_list = genes
        i = len(gene_list)
        for index, row in df.iterrows():
            if i>=(len(genes)+N):
                break
            gene_list.append(row['Gene1'])
            gene_list.append(row['Gene2'])
            gene_list = list(set(gene_list))
            i = len(gene_list)
         
        gene_list = [int(x) for x in gene_list]
        gene_list = [str(x) for x in gene_list]
        gene_list.sort()

        return gene_list 

    return result

"""Helpers"""

def get_column(some_tissue, table='BigGIM_70_v1'):
    columns = []
    try:
		# query
        md = get(f'metadata/tissue/{some_tissue}')
        for st in md['substudies']:
            for column in st.get('columns',[]):
				#filter table
                if column.get('table', {}).get('name', None) == table:
					#filter r(spearman)
#                    if column.get('interactions_type',None) == 'Spearman Rank Correlation Coefficient':
						#get the name of the columns
                        if column.get('name', None) is not None:
                            columns.append(column.get('name'))
    except HTTPError as e:
		#if 404 error, it could not find the tissue, otherwise something bad happend
        if e.args[0].find('404') == -1:
            raise
    return columns







# A few helper functions for posting and getting api requests
#a couple of simple helper functions
def post(endpoint, data={}, base_url=base_url):
    req = requests.post('%s/%s' % (base_url,endpoint), data=data)
    req.raise_for_status()
    return req.json()

def get(endpoint, data={}, base_url=base_url):
    req = requests.get('%s/%s' % (base_url,endpoint), data=data)
    try:
        req.raise_for_status()
    except requests.HTTPError as e:
        print("Sent: GET %s?%s" % (req.request.url,req.request.body))
        if e.response.status_code == 400:
            print(e.response)
            jprint(e.response.json())
        raise
    print("Sent: GET %s?%s" % (req.request.url,req.request.body))
    return req.json()

def jprint(dct):
    print(json.dumps(dct, indent=2))
    
def wrapper(endpoint, data={}, base_url=base_url):
    try:
        response = get(endpoint, data, base_url)
        jprint(response)
    except requests.HTTPError as e:

        print(e)
        if e.response.status_code == 400:
            jprint(e.response.json())
        raise
    try:
        ctr = 1
        while True:
            query_status = get('%s/status/%s'% (endpoint.split('/')[0],response['request_id'],))
            jprint(query_status)
            if query_status['status'] !='running':
                # query has finished
                break
            else:
                time.sleep(ctr)
                ctr += 5
                #linear backoff
                print("Checking again")
    except requests.HTTPError as e:
        print(e)
        if e.response.status_code == 400:
            jprint(e.response.json())
        raise
    return pandas.concat(map(pandas.read_csv, query_status['request_uri']))