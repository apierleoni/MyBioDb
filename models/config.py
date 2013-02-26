class Limits(object):
    '''maximum number of results to report upon a search '''
    max_query_results=500000
    '''maximum number of results to report upon a search '''
    max_dbname_search=1000
    '''maximum number of results to report upon a search '''
    max_qualifier_search=10000
    '''maximum number of entry loadable at one time'''
    max_entry_load = 1000000


class CacheTimes(object):

    '''seconds to store a query in cache '''
    search=300
