__author__ = 'pierleonia'



def index():

    def parse_query_from_form(vars):
        d = dict()
        for key in vars:
            if key.startswith('query[query]'):
                Lkey = key.split('[')
                ID = Lkey[2].split(']')[0]
                try:
                    ID = int(ID)
                    field = Lkey[3].split(']')[0]
                    value = vars[key]
                    d[str(ID)+'-'+field]=value
                except ValueError:
                    pass

        if 'biodatabase' in vars:
            d['biodatabase'] = vars['biodatabase']
        return d

    if request.vars:
        result_table = TABLE(_id='search-results', _class = 'table table-striped table-bordered table-hover')
        if 'query' in request.vars:# simple search
            return dict(content_body=result_table)
        else:
            return dict(content_body= DIV(SPAN('not implemented', _class = 'alert alert-error'), H3('your query:'),P(BEAUTIFY(request.vars)),result_table))
    else:# show advanced query builder


        biodbsel = SELECT('ALL', _id = 'biodatabase', _name = 'biodatabase')
        for row in biodb(biodb.biodatabase.id>0).select(biodb.biodatabase.name):
            biodbsel.append(OPTION(row.name.replace('_',' '), _value = row.name))

        form = FORM(P('Restrict search to entries of type:', biodbsel),
            TABLE(TBODY(TR(TD(SELECT('AND', 'OR', 'NOT', _id = 'logic', _name = 'logic'),
                ' Search by:',SELECT('Name', 'Description',  'Accession',  'Annotation type', 'Annotation value', 'Feature type',
                    'DBXref', 'Sequence',  'Creation date',  'Last modification date',
                    _id = 'field', _name = 'field'), _style = 'text-align: right;'),#hidden available fields: 'ID','Full Text','Keyword','Created by','Last modified by'
                TD(SELECT('Contains','Equal to','Starts with', 'Greater than','Lower than',_id = 'operator', _name = 'operator',),
                    INPUT(_id = 'query', _name = 'query',),
                    DIV(_class = 'validation badge badge-success',  _name = 'validation',),),
                TD(DIV(A(TAG['i'](_class = "icon-minus"),_href = '#', _class = "btn btn-small pull-right",_id ='removequery',),
                    A(TAG['i'](_class = "icon-plus"),_href = '#', _class = "btn btn-small pull-right",_id ='addquery',),_class="btn-group")),
                _id = "duplicablequery",),),
                _id = "querybuilder"),
            P(INPUT(_type = 'submit', _value = 'Search', _class = "btn"), INPUT(_type = 'reset', _value = 'Clear fields', _class = "btn")),
            _name = 'form')

        script = SCRIPT('''

function validate_query(){
    var field = $(this).parent().parent().children("td:nth-child(1)").children("select[id^='field']").val();
    var operator = $(this).parent().children("select[id^='operator']").val();
    var biodatabase = $("#biodatabase").val();
    var query = $(this).val();
    $(this).parent().children('.validation').load("%s",  {field: field,
                                                          operator: operator,
                                                          biodatabase : biodatabase,
                                                          query : query});
    };
function validate_operator(){
    var field = $(this).parent().parent().children("td:nth-child(1)").children("select[id^='field']").val();
    var operator = $(this).val();
    var biodatabase = $("#biodatabase").val();
    var query = $(this).parent().children("input[id^='query']").val();
    $(this).parent().children('.validation').load("%s",  {field: field,
                                                          operator: operator,
                                                          biodatabase : biodatabase,
                                                          query : query});
    };
function validate_field(){
    var field = $(this).val();
    var operator = $(this).parent().parent().children("td:nth-child(2)").children("select[id^='operator']").val();
    var biodatabase = $("#biodatabase").val();
    var query = $(this).parent().parent().children("td:nth-child(2)").children("input[id^='query']").val();
    $(this).parent().parent().children("td:nth-child(2)").children('.validation').load("%s",  {field: field,
                                                                                              operator: operator,
                                                                                              biodatabase : biodatabase,
                                                                                              query : query});
    };


$(document).ready(function(){
    $("#duplicablequery").dynamicForm("#addquery", "#removequery",
        {   formPrefix:"query",
            limit:5,
            createColor: 'navy',
            removeColor: 'red'
        }
    );
    $("#logic0").hide();


    $("input[id^='query']").live('keyup', validate_query);
    $("select[id^='operator']").live('change', validate_operator);
    $("select[id^='field']").live('change', validate_field);


});'''%(URL(r=request, f = "_validate_query"),
        URL(r=request, f = "_validate_query"),
        URL(r=request, f = "_validate_query"),))


        form.append(script)

        '''workaround for bad renaming of jquery-dynamic-form plugin '''
        if 'query[query][_formkey]' in request.vars.keys():
            request.vars['_formkey'] = request.vars['query[query][_formkey]']
        if 'query[query][_formname]' in request.vars.keys():
            request.vars['_formname'] = request.vars['query[query][_formname]']
        if 'query[query][biodatabase]' in request.vars.keys():
            request.vars['biodatabase'] = request.vars['query[query][biodatabase]']

        if form.accepts(request.vars,session):
            vars = parse_query_from_form(request.vars)
            redirect(URL(r=request, f='search', vars=vars))

        return dict(content_body= form)


def typeahead():
    result_count = 0
    maxresult = 8
    if request.vars.query:
        search = BioSQLSearch(biodb_handler)
        query_obj =  search.quick
        result_count = query_obj.contains(request.vars.query)
        if result_count:
            result_sql = query_obj._select()
            return dict(options = get_typehead_results(result_sql, maxresult))
    
    
def search_handler():
    '''initialize search object '''
    session.forget(response)
    from datetime import datetime
    if request.vars.biodatabase:
        search = BioSQLSearch(biodb_handler, biodatabase = biodatabase)
    else:
        search = BioSQLSearch(biodb_handler)
    data = []
    result_ids = []
    starttime = datetime.now()      #debug REMOVE IN PRODUCTION
    if request.vars.query != None:
        query_obj =  search.quick
        result_count = query_obj.contains(request.vars.query)
        if result_count:
            result_sql = query_obj._select()
    if result_count:
        if result_count > Limits.max_query_results:
            response.flash= '''WARNING: you query retrieved %i results, but only the first %i will be displayed. Please narrow your search to see other results'''%(result_count,Limits.max_query_results)
        data.extend(get_search_result_table_from_ids(result_sql))
        print 'populated table returned', datetime.now() - starttime     #debug REMOVE IN PRODUCTION


    return dict(aaData = data)


def get_search_result_table_from_ids(sql_query):
    '''
    Get a list of bioentry ids and return a list of objects to be served as json to a datatable visualizer
    Optimized for speed.
    Bypass web2py parsing and use directly executesql to get results.

    this can be slow for long results in mysql using pymysql because of a very slow call to
    connections.py:892(_read_rowdata_packet)
    something like 10s on pymysql vs 0.4s on MySQLWorkbench client
    if it is an issue try using MySQLdb also named mysql-python

    '''
    data = []
    bioentry_link = URL(r=request, f= 'view.html', vars=dict(bioentry_id = ''))# define just once for speed improvement
    query =  biodb(biodb.bioentry.bioentry_id.belongs(sql_query))._select(biodb.bioentry.accession,
        biodb.bioentry.name,
        biodb.bioentry.description,
        biodb.bioentry.bioentry_id,
        #limitby = (0,Limits.max_query_results),
        #cacheable=True,
    )
    for accession, name, description, bioentry_id in biodb.executesql(query):
        data.append([A(name,
            _href = bioentry_link+str(bioentry_id),
            _class = 'label'),
                     accession,
                     description])

    return data

def get_typehead_results(sql_query, maxresult):
    '''
    

    '''
    data = []
    bioentry_link = URL(r=request, f= 'view.html', vars=dict(bioentry_id = ''))# define just once for speed improvement
    query =  biodb(biodb.bioentry.bioentry_id.belongs(sql_query))._select(biodb.bioentry.accession,
        biodb.bioentry.name,
        biodb.bioentry.description,
        biodb.bioentry.bioentry_id,
        limitby = (0,maxresult),
        #cacheable=True,
    )
    for accession, name, description, bioentry_id in biodb.executesql(query):
        data.append(A(name+" [%s]"%accession,
            _href = bioentry_link+str(bioentry_id),
            _class = 'typeahead-item',
                     ),
                     )

    return data

def _validate_query():


    field = request.vars.field
    operator = request.vars.operator
    query = request.vars.query
    biodatabase = request.vars.biodatabase
    if biodatabase == 'ALL':
        biodatabase = None

    if field and operator and query :

        mapping_dict = {'Name' : 'name',
                        'Description' : 'description',
                        'Full Text' : 'fulltext',
                        'Quick' : 'quick',
                        'Accession' : 'accession',
                        'ID' : 'id',
                        'Annotation type' : 'annotation_type',
                        'Annotation value' : 'annotation_value',
                        'Feature type' : 'feature_type',
                        'DBXref' : 'dbxref',
                        'Sequence' : 'seq',
                        'Keyword' : 'keyword',
                        'Creation date' : 'created_on',
                        'Created by' : 'created_by',
                        'Last modification date' : 'modified_on',
                        'Last modified by' : 'modified_by',}

        '''initialize search object '''
        search = BioSQLSearch(biodb_handler, biodatabase = biodatabase)
        query_obj = search.__dict__[mapping_dict[field]]

        '''query values'''
        if query and (query[-1] == ';'):
            query = query[:-1]
        query_values = query.split(';')

        '''query type'''
        try:
            if operator == 'Equal to':
                query_count = query_obj == query_values
            elif operator == 'Contains':
                query_count = query_obj.contains(query_values)
            elif operator == 'Starts with':
                query_count = query_obj.startswith(query_values)
            elif operator == 'Greater than':
                query_count = query_obj > query_values
            elif operator == 'Lower than':
                query_count = query_obj < query_values

            if field not in ['Name', 'Description', 'Quick', 'Accession', 'Accession', 'Sequence',
                             'Created by', 'Creation date', 'Last modification date', 'Last modified by' ]:#cannot use query_count since there is no groupby in it
                query_count = query_obj.count()
            return query_count
        except NotImplementedError:
            return SPAN('Wrong operator', _style = 'color:#970000;')
    else:
        return 0
