# -*- coding: utf-8 -*-


response.generic_patterns = ['*'] 

### required - do no delete
def user(): return dict(form=auth())
def download(): return response.download(request,db)
def call(): return service()
### end requires
def index():
    return dict()

def error():
    return dict()
    
def _get_uniprot_xml(uniprot_id):
    import urllib
    from Bio.SeqIO.UniprotIO import UniprotIterator 
    seqrec = cache.disk(str('%s-uniprotxml'%uniprot_id), 
                        lambda: UniprotIterator(urllib.urlopen('http://www.ebi.ac.uk/Tools/webservices/rest/dbfetch/uniprotkb/%s/uniprotxml'%uniprot_id),).next(), 
                        time_expire=180000)#cache on disk for two days (more or less)
    return seqrec #returns a StringIO handler
    
def uniprot_remote_load():
    form=FORM('Please input a UniProt ID to load:',
              INPUT(_name = 'uniprot_id'),
              INPUT(_type = 'submit', _value = 'retrieve from UniProt'),)
    if form.accepts(request.vars, session):
        seqrec = _get_uniprot_xml(form.vars.uniprot_id)
        try:
            created_seqrec_id = biodb_handler.load_seqrecord(seqrec, db = 'UniProt')
            session.flash = 'Entry correctly loaded on BioSQL db with id: %i'%created_seqrec_id
            redirect(URL(r=request, f = 'view', vars = dict(bioentry_id = created_seqrec_id)))
        except Exception, e:
            response.flash = 'Error: %s'%e.message
    return dict(form = form, )


def import_entry():

    from Bio import SeqIO

    biodbs = SELECT(_name = 'biodb')
    for row in  biodb(biodb.biodatabase.id>0).select():
        biodbs.append(OPTION(row.name))#, _value = row.biodatabase_id))
    form=FORM('Please select a file to upload and a supported format:',
        INPUT(_type='file',_name = 'file', requires=IS_NOT_EMPTY()),
        SELECT(('uniprot-xml','genbank', 'fasta' ),_name = 'format'),BR(),
        'Please select the database:',
        biodbs,BR(),
        INPUT(_type = 'submit', _value = 'upload'),)
    if form.accepts(request.vars, session):
        try:
            created_seqrec = []
            errors =[]
            try:
                for seqrec in SeqIO.parse(form.vars.file.file, form.vars.format):
                    try:
                        created_seqrec_id = biodb_handler.load_seqrecord(seqrec, db = form.vars.biodb)
                        created_seqrec.append((seqrec.id, created_seqrec_id))
                    except Exception, e:
                        errors.append((seqrec.id, e))
                    response.flash = '%i entries correctly loaded on db %s'%(len(created_seqrec),
                                                                           form.vars.biodb )
                    if errors:
                        response.flash += '| %i errors while loading'%len(errors)
                    returndiv = DIV(H3('Loaded entries:'))
                    for entry in created_seqrec:
                        returndiv.append(LI(A(entry[0], _href = URL(r=request, f = 'view', vars = dict(bioentry_id = entry[1])))))
                    if errors:
                        returndiv.append(H3('Errors:'))
                        for entry in errors:
                            returndiv.append(LI(entry[0],': ', entry[1]))
                    return dict(main_content = returndiv)
            except:
                response.flash = 'PARSING FAILED: please check that all the entries are in the specified format'

        except Exception, e:
            response.flash = 'Error: %s'%e.message
    return dict(form = form, )
    

def view_uniprot_entries():
    returnul = UL()
    for row in biodb(biodb.bioentry.id>0).select():
        returnul.append(LI(A(row.name,_href = URL(r=request, f = view, vars = dict(bioentry_id = row.id )))))
     
    return dict(returnul = returnul)
   
def input_error_redirect(error):
    session.flash = str(error)
    redirect(URL(r=request, f = 'index'))
    
def view():    
    main_content = DIV(_id = 'maincontent')
    if request.vars.bioentry_id:
        bioentry_id = int(request.vars.bioentry_id)
        biodb_row = biodb((biodb.biodatabase.biodatabase_id == biodb.bioentry.biodatabase_id) &
                           (biodb.bioentry.bioentry_id == bioentry_id)).select(biodb.biodatabase.name).first()
        if biodb_row:
            biodb_name = biodb_row.name
        else:
            input_error_redirect('ERROR! invalid bioentry data')
    elif request.args:# reads: view/dbname/accession
        try:
            biodb_name = request.args[0]
            bioentry_id = biodb_handler._get_bioentry_id(accession = request.args[1], dbid = biodb_handler.dbs[biodb_name])
        except:
            input_error_redirect('ERROR! invalid bioentry data')
    else:
        input_error_redirect('ERROR! invalid bioentry data')
    editable = False
    if request.vars.editable == 'True':
        editable = True
    sidebar = UL(_class = "nav nav-list sidenav affix", _id = 'left-sidebar', **{"data-spy":"affix", "data-offset-top":"100"})
    if not request.vars.views:
        request.vars.views =  ['v_general', 'v_taxonomy', 'v_relationships',
                               'v_qualifiers','v_comments', 'v_sequence',
                               'v_features', 'v_dbxrefs','v_references',   ]
    for view in request.vars.views:
        main_content.append(LOAD(c='default', f=view,
                                 vars={'bioentry_id' : bioentry_id},
                                 ajax=False,
                                 ajaxTrap=True,
                                 content = IMG(_src = URL(request.application, 'static', 'ajax-loader.gif')),
                                 #content = '',
                                 ))
        view_key = ''.join(view.split('v_')[1:])

        sidebar.append(LI(A(' '+view_key, TAG.i(_class="icon-chevron-right"),_href='#'+view_key, _class = 'affix-element')))

    if not editable:
        sidebar.append(LI(A('Edit',
                            _href=URL(r = request, f= 'view', vars =dict(bioentry_id =bioentry_id, editable = True))),
                            _class = 'btn btn-primary btn-block'))
    else:
        sidebar.append(LI(A('Done editing',
                            _href=URL(r = request, f= 'view', vars =dict(bioentry_id =bioentry_id)),
                            _class = 'btn btn-primary btn-block')))



    return dict(main_content = main_content, bioentry_id = bioentry_id, editable = editable, sidebar = sidebar)
    
    
def v_general():
    if request.vars.bioentry_id:
        bioentry_id = int(request.vars.bioentry_id)
        returndiv = DIV(_class = 'unit-view-element')
        row = biodb(biodb.bioentry.id == bioentry_id).select().first()
        
        title = SPAN(row.name,
                     _class = 'edit-bioentry-name editable', 
                     _id = 'bioentry-name', 
                     _title = 'Entry name, click to edit')
        
        content = DIV()
        content.append(H3( SPAN(row.description,
                                              _class = 'edit-bioentry-description editable' , 
                                              _id = 'bioentry-description', 
                                              _title = 'Entry description, click to edit' ))
                       )
        content.append(H4('DB: %s | Accession: %s | id: %i'%(biodb.biodatabase[row.biodatabase_id].name,
                                                         row.accession,
                                                         row.bioentry_id))
                       )
        
        '''timestamps and author '''
        year=day=month = ''
        timestamps = biodb_handler._get_bioentry_time_stamp(bioentry_id)
        author = ''
        try:
            if timestamps:
                author= '%(first_name)s %(last_name)s'%db.auth_user[timestamps['created_by']]
        except:
            pass
        footer=P()
        '''
        if 'date' in seqrec.annotations:
            if isinstance(seqrec.annotations['date'], list):
                year, month, day = seqrec.annotations['date'][0].split('-')
            elif isinstance(seqrec.annotations['date'], str):
                year, month, day = seqrec.annotations['date'].split('-')
        elif timestamps:
            #date_footer.append('Record Date Added: ' + timestamps['created_on'].date().isoformat())
            year = timestamps['created_on'].year
            month = timestamps['created_on'].month
            day = timestamps['created_on'].day
        '''
        if timestamps:
            year = timestamps['created_on'].year
            month = timestamps['created_on'].month
            day = timestamps['created_on'].day
            footer.append('Record Last Modified: '+ timestamps['modified_on'].date().isoformat())
            try:
                if timestamps:
                    footer.append(' by %(first_name)s %(last_name)s'%db.auth_user[timestamps['modified_by']])
                    userdata = db.auth_user[timestamps['created_by']]
            except:
                pass

        
        
        return UnitView(title = title,
                        content = content,
                        counts = '',
                        footer = footer,
                        author = '',
                        _id = 'general',
                        _class = '')
    return DIV()
    
    
def v_taxonomy():
    '''if request.vars.bioentry_id:
        content = DIV()
        return UnitView(title = 'Taxonomy',
                        content = content,
                        counts = '',
                        footer = '',
                        author = '',
                        _class = 'collapsible ')'''
    return DIV()
    
def v_relationships():
    '''if request.vars.bioentry_id:
        content = DIV()
        return UnitView(title = 'Relationships',
                        content = content,
                        counts = '',
                        footer = '',
                        author = '',
                        _class = 'collapsible ')'''
    return DIV()
    
def v_qualifiers():
    if request.vars.bioentry_id:
        content = DIV()
        qualifiers = BioSQLBioentryQualifiers(biodb_handler,int(request.vars.bioentry_id))
        if not qualifiers.simple_qualifiers:
            return DIV()
        else:
            qaltable=TABLE(THEAD(TR(TH('Type'),TH('Annotation'),
                              _style = 'font-weight: bold')),
                        _class ='table table-striped table-bordered table-hover  datatable')
            qualtbody=TBODY()
            for q in qualifiers.qualifiers:
                key = q.key
                bullets = len(q.rank2value) > 1
                values = UL()
                if bullets:
                    for rank,value in sorted(q.rank2value.items()):
                        values.append(LI(SPAN(value, _class = 'edit-field-qualifier editable', _rank = rank)))
                    qualtbody.append(TR(key,values))
                else:
                    rank,value = q.rank2value.items()[0]
                    qualtbody.append(TR(TD(key), TD(SPAN(value, _class = 'edit-field-qualifier editable', _rank = rank))))
                   
            #for k,v in sorted(qualifiers.simple_qualifiers.items()):
            #    qualtbody.append(TR(k,biosql_data_beautify(v)))
            qaltable.append(qualtbody)
            return UnitView(title = 'General annotations (Qualifiers)',
                            content = qaltable,
                            counts = len(qualifiers.qualifiers),
                            footer = '',
                            author = '',
                            _id = 'qualifiers',
                            _class = 'collapsible ')
    return DIV()
  
def v_comments():
    '''if request.vars.bioentry_id:
        content = DIV()
        return UnitView(title = 'Comments',
                        content = content,
                        counts = '',
                        footer = '',
                        author = '',
                        _class = 'collapsible ')'''
    return DIV()
    


def v_features():
    if request.vars.bioentry_id:
        features = BioSQLAlterBioentryFeatures(biodb_handler,int(request.vars.bioentry_id)).get_seqfeatures()
        seq_length = biodb(biodb.biosequence.bioentry_id == int(request.vars.bioentry_id)).select(biodb.biosequence.length).first().length
        feattable=TABLE(THEAD(TR(TH('Type'),TH('Position(s)'),TH('Length'),
                                 TH('Description'),TH('Graphical view'),
                                 TH('Qualifiers'), TH('Identifier'),
                              _style = 'font-weight: bold')),
                        _id = 'feature-table',
                        _class ='table table-striped table-bordered table-hover datatable-not-sorted')
        featbody=TBODY()
        for feat_id,feat in features:
            featbody.append(draw_single_feature(feat, seq_length, feat_id,))
        feattable.append(featbody)
        control_btn = DIV()
        add_btn = A('Add new feature',
                         _href = URL(r = request, 
                                     f = 'add_feature', 
                                     vars = dict(bioentry_id =request.vars.bioentry_id)),
                         _class = 'btn')
        control_btn.append(add_btn)
        delete_btn = SPAN('Delete selected feature',
            _id = 'feature-delete',
            _class = 'btn disabled')
        control_btn.append(delete_btn)
        content = DIV(H4('Feature table'),feattable, control_btn)#, view_feature_dialog)
        return UnitView(title = 'Features',
                        content = content,
                        counts = len(features),
                        footer = '',
                        author = '',
                        _class = 'collapsible ',
                        _id = 'features')
    return DIV()

def draw_feature_on_sequence():
    if request.vars.bioentry_id:
        sequence = biodb_handler._retrieve_seq(int(request.vars.bioentry_id))
        if sequence:
            if request.vars.start:
                request.vars.start = int(request.vars.start)-1
            else:
                request.vars.start = 0
            if (not request.vars.end) or (int(request.vars.end)<=request.vars.start):
                request.vars.end = request.vars.start+1
            else:
                request.vars.end=int(request.vars.end)
            return draw_feature_on_sequence_utility(sequence, 
                                                    request.vars.start, 
                                                    request.vars.end)
            
def v_dbxrefs():
    if request.vars.bioentry_id:
        content = DIV()
        dbxrefs = BioSQLAlterBioentryDBXrefs(biodb_handler).read(int(request.vars.bioentry_id))
        if not dbxrefs:
            return DIV()
        else:
            dbxreftable=TABLE(THEAD(TR(TH('database'),TH('Entry'),
                              _style = 'font-weight: bold')),
                              _id = 'dbxref-table',
                              _class ='table table-striped table-bordered table-hover datatable')
            dbxreftbody=TBODY()
            for dbxref in dbxrefs:
                dbxreftbody.append(draw_single_dbxref(*dbxref))
            dbxreftable.append(dbxreftbody)
            control_panel = DIV(_class = 'dbxref-control-panel')
            add_btn = A('Add new dbxref',
                         _href = URL(r = request,
                                     f = 'add_dbxref',
                                     vars = dict(bioentry_id =request.vars.bioentry_id)),
                         _class = 'btn ')
            control_panel.append(add_btn)
            delete_btn = SPAN('Delete selected dbxref',
                              _id = 'dbxref-delete',
                              _class = 'btn disabled')
            control_panel.append(delete_btn)

            content = DIV(dbxreftable,control_panel)
            return UnitView(title = 'Database cross-references (DBXrefs)',
                            content = content,
                            counts = len(dbxrefs),
                            footer = '',
                            author = '',
                            _class = 'collapsible ',
                            _id = 'dbxrefs')
    return DIV()
    
def v_references():
    if request.vars.bioentry_id:
        content = DIV()
        refs = BioSQLReference(biodb_handler,int(request.vars.bioentry_id))
        if not refs.references:
            return DIV()
        else:
            refstable=TABLE(THEAD(TR(TH('year'),TH('reference'),
                              _style = 'font-weight: bold')),
                        _class ='table table-striped table-bordered table-hover datatable-inverted')
            refstbody=TBODY()
            for ref in refs.references:
                try:
                    year = ref.journal.split('(')[-1].split(')')[0]
                except:
                    year = ''
                refstbody.append(TR(year, draw_single_reference(ref)))
            refstable.append(refstbody)
            add_new = P(A('Add new reference',
                         _href = URL(r = request, 
                                     f = 'add_reference', 
                                     vars = dict(bioentry_id =request.vars.bioentry_id)),
                         _class = 'btn') )
            content = DIV(refstable,add_new)
            return UnitView(title = 'References',
                            content = content,
                            counts =  len(refs.references),
                            footer = '',
                            author = '',
                            _id = 'references',
                            _class = 'collapsible ')
    return DIV()
    
def v_sequence():
    if request.vars.bioentry_id:
        returndiv = DIV()
        seq = BioSQLSeq(biodb_handler, int(request.vars.bioentry_id))
        sequence = seq.seq.tostring()
        if seq.alphabet == 'protein':
            returndiv.append(draw_sequence(sequence, mode = 'protparams'))
        else:
            returndiv.append(draw_sequence(sequence))
        
        return UnitView(title = 'Sequence',
                        content = returndiv,
                        counts = 'Length:' +str(len(sequence)),
                        footer = '',
                        author = '',
                        _id = 'sequence',
                        _class = 'collapsible ')
    return DIV()
    

def return_pubmed_data():
    if request.vars.pmid:
        data = fetch_pubmed_data(request.vars.pmid)
        controldiv = DIV(A('hide' ,
                            _class = 'control-pubmed-data btn btn-mini pull-right',
                            _onclick = 'void(0)'),)
        returndiv = DIV(_class= 'jref-additional-data')
        if data:
            abstract = data['AB']
            links = data['AID']
            mesh_terms = data['MH']
            if links:
                for l in links:
                    if '[doi]' in l:
                        doi = l.split('[doi]')[0].strip()
                        returndiv.append(P(A('Go to journal',
                                             _href = 'http://dx.doi.org/%s'%doi,
                                             _target = '_blank', 
                                             _class ='btn btn-primary btn-small')))

            if abstract:
                returndiv.append(DIV(H4('Abstract'), P(abstract, _class ='jref-abstract')))
            if mesh_terms:
                terms = P(_class ='jref-mesh-terms')
                for term in mesh_terms:
                    terms.append(A(term,
                                   _href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=mesh&term=%%22%s%%22'%term,
                                   _target = '_blank',
                                   _class ='label label-info mesh-tags'))

                returndiv.append(DIV(H4('Mesh Terms'), terms))
        
            return DIV(controldiv, returndiv)
        else:
            return 'Unable to fetch data from ncbi'
    
#@auth.requires_login()         
def alter_single_qualifier_value():
    '''test how it works with annotations consisting of multiple lines'''
    if request.vars.bioentry_id:
        print request.vars
        bioentry_id = int(request.vars.bioentry_id)
        key = request.vars.id
        value = request.vars.value
        rank = request.vars.rank
        try:
            ann = BioSQLQualifier(biodb_handler, bioentry_id=bioentry_id, key = key)
            ann.update_single_value(value = value, rank = rank)
            '''while old_value in  ann.value:
                biodb.value.pop(ann.value.index(old_value))
            ann.value.append(value)
            ann.sync()
            biodb_handler._update_bioentry_time_stamp(bioentry_id)'''
            return value
        except:
            response.headers['Content-Type'] = None
            response.status = 400
            return 'Error'

        
#@auth.requires_login() 
def alter_bioentry():
    if request.vars.bioentry_id:
        bioentry_id = int(request.vars.bioentry_id)
        value = request.vars.value
        try:
            if request.vars.type == 'name':
                biodb.bioentry[bioentry_id] = dict(name = request.vars.value)
            elif request.vars.type == 'description':
                biodb.bioentry[bioentry_id] = dict(description = request.vars.value)
            elif request.vars.type == 'division':
                biodb.bioentry[bioentry_id] = dict(division = request.vars.value)
            elif request.vars.type == 'accession':
                biodb.bioentry[bioentry_id] = dict(accession = request.vars.value)
            elif request.vars.type == 'version':
                biodb.bioentry[bioentry_id] = dict(version = int(request.vars.value))
            elif request.vars.type == 'identifier':
                biodb.bioentry[bioentry_id] = dict(identifier = request.vars.value)
            elif request.vars.type == 'taxon_id':
                biodb.bioentry[bioentry_id] = dict(taxon_id = int(request.vars.value))
            elif request.vars.type == 'biodatabase_id':
                biodb.bioentry[bioentry_id] = dict(biodatabase_id = int(request.vars.value))

            biodb_handler._update_bioentry_time_stamp(bioentry_id)
            return value
        except:
            response.headers['Content-Type'] = None
            response.status = 400
            return 'Error'

def alter_single_feature_value():
#    print 'feature alter values', request.vars
    if request.vars.pk:
        try:
            # read the requeste feature, change the requested value and updated it
            # can be made more efficient by writing single sql update for every type
            if request.vars.type == 'type':
                BioSQLAlterFeature(biodb_handler).update_seqfeature_type(
                                int(request.vars.pk),
                                request.vars.value)
            elif request.vars.type == 'start':
                BioSQLAlterFeature(biodb_handler).update_seqfeature_start(
                    int(request.vars.pk),
                    request.vars.value)
            elif request.vars.type == 'end':
                BioSQLAlterFeature(biodb_handler).update_seqfeature_end(
                    int(request.vars.pk),
                    request.vars.value)
            elif request.vars.type == 'identifier':
                BioSQLAlterFeature(biodb_handler).update_seqfeature_identifier(
                    int(request.vars.pk),
                    request.vars.value)
            elif request.vars.type == 'description':
                BioSQLAlterFeature(biodb_handler).update_seqfeature_description(
                    int(request.vars.pk),
                    request.vars.value)
            elif request.vars.type == 'qualifier-value':
                if '|' in request.vars.value:
                    request.vars.value = request.vars.value.split('|')
                BioSQLAlterFeature(biodb_handler)._update_qualifier(
                        int(request.vars.pk),
                        request.vars.key,
                        request.vars.value)

            return
            #biodb_handler._update_bioentry_time_stamp(bioentry_id)
        except:
            response.headers['Content-Type'] = None
            response.status = 400
            return 'Error'

    response.headers['Content-Type'] = None
    response.status = 501
    return 'Not implemented'+str(request.vars)


def alter_single_dbxref():
    if request.vars.pk:
        try:
            'TODO: change this to remove dbxref from bioentry and add/link a new one'
            if request.vars.type == 'dbname':
                BioSQLAlterBioentryDBXref(biodb_handler).update_dbname(int(request.vars.pk),request.vars.type.value);
            elif request.vars.type == 'accession':
                BioSQLAlterBioentryDBXref(biodb_handler).update_dbname(int(request.vars.pk),request.vars.type.value);
        except:
            response.headers['Content-Type'] = None
            response.status = 400
            return 'Error'
    response.headers['Content-Type'] = None
    response.status = 501
    return 'Not implemented'


def delete_single_dbxref():
    if request.vars.bioentry_id and request.vars.dbxref_id:
        try:
            BioSQLAlterBioentryDBXref(biodb_handler).delete(int(request.vars.bioentry_id),int(request.vars.dbxref_id));
            response.flash = 'dbxref succesfully deleted'
            return 'dbxref succesfully deleted'
        except:
            response.headers['Content-Type'] = None
            response.status = 400
            response.flash = 'Error'
            return 'Error'
    response.headers['Content-Type'] = None
    response.status = 501
    response.flash = 'Error'
    return 'Not implemented'

def delete_single_feature():
    if request.vars.seqfeature_id :
#        try:
            BioSQLAlterFeature(biodb_handler).delete(int(request.vars.seqfeature_id));
            return 'feature succesfully deleted'
#        except:
#            response.headers['Content-Type'] = None
#            response.status = 400
#            return 'Error'
    response.headers['Content-Type'] = None
    response.status = 501
    return 'Not implemented'

def add_feature():
    if request.vars.bioentry_id:
        feat_types = return_feature_type_select(biodb)
        bioentry_id = int(request.vars.bioentry_id)
        max_end_position = biodb(biodb.biosequence.bioentry_id == bioentry_id).select(biodb.biosequence.length).first().length
        drawnsequence = LOAD(request.controller,
                       'draw_feature_on_sequence',
                       vars=dict(bioentry_id=bioentry_id, start = 1, mode= 'add feature'),
                       ajax = False,
                       ajaxTrap=True,)
        seqdiv = DIV(drawnsequence,
                   _class = 'seqdiv',
                   _id = 'seqdiv-add-feature')
        #seqdiv = DIV(LOAD(request.controller,'draw_feature_on_sequence',vars=dict(bioentry_id=request.vars.bioentry_id, start = 1),ajax=True), _id='seqdiv')
        content_body=DIV(H4('Sequence'), seqdiv)
        form=FORM(H4('Set positions:'),
                  P('Start position: ',INPUT(_name = 'start', _id = 'startposition', _value = 1, _style= 'width: 50px',
                                         requires= (IS_INT_IN_RANGE(1,max_end_position), IS_NOT_EMPTY()),),
                  ' End position: ',INPUT(_name = 'end', _id = 'endposition', _style= 'width: 50px',
                                         requires= IS_NULL_OR(IS_INT_IN_RANGE(1,max_end_position))),'(specified positions are included)'),
                  H4('Set feature data:'),
                  P('Type: ',SPAN(feat_types,
                              A('new',_id ='shownewvalue', _class = 'btn',_onclick = 'void(0)'), 
                              DIV(INPUT(_name='newvalue', _id="newvalue"),
                              A('add',_id ='addnewvalue', _class = 'btn', _onclick = 'void(0)'),_id = 'newvaluecontainer'))),
                  P('Description: ',INPUT(_name='description')),
                  P('Identifier: ',INPUT(_name='identifier')),
                  INPUT(_type = 'submit', _value = 'Add feature')#,INPUT(_type = 'reset', ),
                  )
        content_body.append(form)
        seqrefresh = SCRIPT('''
$(document).ready(function(){

//handle new feature value
$('#newvaluecontainer').toggle();
$('#shownewvalue').click(function() {
    $("#newvaluecontainer").slideToggle("slow");
    var text = $(this).text() == 'new' ? 'hide' : 'new';
    $(this).text(text)
    });
$('#addnewvalue').click(function() {
    $('#featuretype-select').append( new Option('NEW: '+$('#newvalue').val(),$('#newvalue').val()));
    $("#featuretype-select").val($('#newvalue').val());
    $('#newvaluecontainer').toggle();
    $('#shownewvalue').text('new');
    });

//live draw
$("#startposition").change(update);
$("#startposition").keyup(update);
$("#endposition").change(update);
$("#endposition").keyup(update);
});

function update(){      
    //$('#seqdiv-add-feature').slideDown('fast');
    var startvalue = $("#startposition").val();
    var endvalue = $("#endposition").val();
$('#seqdiv-add-feature').load("%s",  {bioentry_id: %s, start: startvalue, end : endvalue, mode : 'add feature'});
}'''%(URL(r=request, f = draw_feature_on_sequence), bioentry_id))
        content_body.append(seqrefresh)
         
        if form.accepts(request.vars, session):
            #biodb_archiever.archive(request.vars.bioentry_id)#archive old entry
            from Bio import SeqFeature 
            feat = SeqFeature.SeqFeature(type =form.vars.type)
            form.vars.start = int(form.vars.start)
            if form.vars.description:
                feat.qualifiers['description'] = form.vars.description
            if form.vars.identifier:
                #feat.id = form.vars.identifier
                feat.qualifiers['id'] = form.vars.identifier
            if not form.vars.end :
                form.vars.end = form.vars.start
            else:
                form.vars.end = int(form.vars.end)
                if form.vars.end <form.vars.start:
                    form.vars.end = form.vars.start
            feat.location = SeqFeature.FeatureLocation(form.vars.start-1, form.vars.end)
            biosqlfeat = BioSQLFeature(biodb_handler, bioentry_id = bioentry_id, feature = feat)
            biosqlfeat.sync()
            session.flash='Feature successfully added'
            redirect(URL(r=request, f='view', vars = dict(bioentry_id = bioentry_id)))  
        elif form.errors:
            response.flash='FORM ERROR!'   

        return dict(main_content = UnitView(title = 'Add feature',
                        content = content_body,))
         
def add_qualifier():
    content_body=DIV()
    return dict(main_content = UnitView(title = 'Add qualifier',
                        content = content_body,))
def add_dbxref():
    content_body=DIV(TODO)
    return dict(main_content = UnitView(title = 'Add database cross reference',
                        content = content_body,))
def add_reference():
    if request.vars.bioentry_id:
        content_body=DIV()
        bioentry_id = int(request.vars.bioentry_id)
        max_end_position = biodb(biodb.biosequence.bioentry_id == bioentry_id).select(biodb.biosequence.length).first().length
        form = FORM('PubMed ID: ', INPUT(_name = 'pmid', requires = IS_NOT_EMPTY()), BR(),
                    'Start position: ', INPUT(_name = 'start', requires = IS_NULL_OR(IS_INT_IN_RANGE(1, max_end_position)),_style = 'width: 80px;'), 
                    ' End position: ', INPUT(_name = 'end', requires = IS_NULL_OR(IS_INT_IN_RANGE(1, max_end_position)),_style = 'width: 80px;'), BR(),
                    'Comment: ', INPUT(_name = 'comment'), BR(),
                     INPUT(_type='submit', _value = 'Import from PubMed'))
        content_body.append(form)
        
        if form.accepts(request.vars, session):
            from Bio import SeqFeature
            data = fetch_pubmed_data(form.vars.pmid)
            reference = SeqFeature.Reference()
            #If the start/end are missing, reference.location is an empty list
            if form.vars.start:
                form.vars.start -= 1 #python counting
                if form.vars.end is None:
                    form.vars.end = form.vars.start
                reference.location = [SeqFeature.FeatureLocation(form.vars.start, form.vars.end)]
            #Don't replace the default "" with None.
            if data['AU']: reference.authors = ', '.join(data['AU'])
            if data['TI']: reference.title = data['TI']
            REFERENCE_JOURNAL = "%(name)s %(volume)s:%(pages)s(%(pub_date)s)"
            
            reference.journal = REFERENCE_JOURNAL % dict(name=data['TA'],
                            volume=data['VI'], pages=data['PG'], pub_date=data['DP'].split(' ')[0])
            reference.pubmed_id = form.vars.pmid
            reference.comment = form.vars.comment
            '''TEMPORARY, make a class to handle single references'''
            dbreferences =BioSQLReference (biodb_handler, bioentry_id)
            dbreferences._insert(reference)
            session.flash = 'Reference correctly added'
            redirect(URL(r=request, f = 'view', vars = dict(bioentry_id = bioentry_id)))
        return dict(main_content = UnitView(title = 'Add reference from PubMed',
                            content = content_body,))



## Static analyzer import helpers for controllers:
#if 0:
#    import gluon
#    from gluon.contrib.gql import GQLDB
#    from gluon.html import A
#    from gluon.html import B
#    from gluon.html import BEAUTIFY
#    from gluon.html import BODY
#    from gluon.html import BR
#    from gluon.html import CENTER
#    from gluon.html import CLEANUP
#    from gluon.html import CODE
#    from gluon.html import CRYPT
#    from gluon.html import DIV
#    from gluon.html import FORM
#    from gluon.html import I
#    from gluon.html import IFRAME
#    from gluon.html import IMG
#    from gluon.html import INPUT
#    from gluon.html import IS_ALPHANUMERIC
#    from gluon.html import IS_DATE
#    from gluon.html import IS_DATETIME
#    from gluon.html import IS_DATETIME_IN_RANGE
#    from gluon.html import IS_DATE_IN_RANGE
#    from gluon.html import IS_DECIMAL_IN_RANGE
#    from gluon.html import IS_EMAIL
#    from gluon.html import IS_EMPTY_OR
#    from gluon.html import IS_EQUAL_TO
#    from gluon.html import IS_EXPR
#    from gluon.html import IS_FLOAT_IN_RANGE
#    from gluon.html import IS_IMAGE
#    from gluon.html import IS_INT_IN_RANGE
#    from gluon.html import IS_IN_DB
#    from gluon.html import IS_IN_SET
#    from gluon.html import IS_IPV4
#    from gluon.html import IS_LENGTH
#    from gluon.html import IS_LIST_OF
#    from gluon.html import IS_LOWER
#    from gluon.html import IS_MATCH
#    from gluon.html import IS_NOT_EMPTY
#    from gluon.html import IS_NOT_IN_DB
#    from gluon.html import IS_NULL_OR
#    from gluon.html import IS_SLUG
#    from gluon.html import IS_STRONG
#    from gluon.html import IS_TIME
#    from gluon.html import IS_UPLOAD_FILENAME
#    from gluon.html import IS_UPPER
#    from gluon.html import IS_URL
#    from gluon.html import LABEL
#    from gluon.html import LEGEND
#    from gluon.html import LI
#    from gluon.html import LINK
#    from gluon.html import MARKMIN
#    from gluon.html import MENU
#    from gluon.html import META
#    from gluon.html import OBJECT
#    from gluon.html import OL
#    from gluon.html import ON
#    from gluon.html import OPTGROUP
#    from gluon.html import OPTION
#    from gluon.html import P
#    from gluon.html import PRE
#    from gluon.html import SCRIPT
#    from gluon.html import SELECT
#    from gluon.html import SPAN
#    from gluon.html import TABLE
#    from gluon.html import TAG
#    from gluon.html import TBODY
#    from gluon.html import TD
#    from gluon.html import TEXTAREA
#    from gluon.html import TFOOT
#    from gluon.html import TH
#    from gluon.html import THEAD
#    from gluon.html import TITLE
#    from gluon.html import TR
#    from gluon.html import TT
#    from gluon.html import UL
#    from gluon.html import URL
#    from gluon.html import XHTML
#    from gluon.html import XML
#    from gluon.html import embed64
#    from gluon.html import xmlescape
#    from gluon.sql import SQLDB
#    from gluon.sqlhtml import SQLFORM
#    from gluon.sql import SQLField
#    from gluon.sqlhtml import SQLTABLE
#    from gluon.html import STYLE
#    from gluon.http import redirect
#    import gluon.languages.translator as T
#    from gluon.tools import Auth
#    from gluon.tools import Service
#    global auth; auth = gluon.tools.Auth()
#    global cache; cache = gluon.cache.Cache()
#    global crud; crud = gluon.tools.Crud()
#    global db; db = gluon.sql.DAL()
#    import gluon.compileapp.local_import_aux as local_import
#    global request; request = gluon.globals.Request()
#    global response; response = gluon.globals.Response()
#    global service; service = gluon.tools.Service()
#    global session; session = gluon.globals.Session()
#    global DAL; DAL = gluon.dal()
#    global HTTP; HTTP = gluon.http()
#    global LOAD; LOAD = gluon.compileapp.LoadFactory()
#
## Static analyzer import helpers for models:
#if 0:
#    None
# import gluon
# from gluon.contrib.gql import GQLDB
# from gluon.html import A
# from gluon.html import B
# from gluon.html import BEAUTIFY
# from gluon.html import BODY
# from gluon.html import BR
# from gluon.html import CENTER
# from gluon.html import CLEANUP
# from gluon.html import CODE
# from gluon.html import CRYPT
# from gluon.html import DIV
# from gluon.html import FORM
# from gluon.html import I
# from gluon.html import IFRAME
# from gluon.html import IMG
# from gluon.html import INPUT
# from gluon.html import IS_ALPHANUMERIC
# from gluon.html import IS_DATE
# from gluon.html import IS_DATETIME
# from gluon.html import IS_DATETIME_IN_RANGE
# from gluon.html import IS_DATE_IN_RANGE
# from gluon.html import IS_DECIMAL_IN_RANGE
# from gluon.html import IS_EMAIL
# from gluon.html import IS_EMPTY_OR
# from gluon.html import IS_EQUAL_TO
# from gluon.html import IS_EXPR
# from gluon.html import IS_FLOAT_IN_RANGE
# from gluon.html import IS_IMAGE
# from gluon.html import IS_INT_IN_RANGE
# from gluon.html import IS_IN_DB
# from gluon.html import IS_IN_SET
# from gluon.html import IS_IPV4
# from gluon.html import IS_LENGTH
# from gluon.html import IS_LIST_OF
# from gluon.html import IS_LOWER
# from gluon.html import IS_MATCH
# from gluon.html import IS_NOT_EMPTY
# from gluon.html import IS_NOT_IN_DB
# from gluon.html import IS_NULL_OR
# from gluon.html import IS_SLUG
# from gluon.html import IS_STRONG
# from gluon.html import IS_TIME
# from gluon.html import IS_UPLOAD_FILENAME
# from gluon.html import IS_UPPER
# from gluon.html import IS_URL
# from gluon.html import LABEL
# from gluon.html import LEGEND
# from gluon.html import LI
# from gluon.html import LINK
# from gluon.html import MARKMIN
# from gluon.html import MENU
# from gluon.html import META
# from gluon.html import OBJECT
# from gluon.html import OL
# from gluon.html import ON
# from gluon.html import OPTGROUP
# from gluon.html import OPTION
# from gluon.html import P
# from gluon.html import PRE
# from gluon.html import SCRIPT
# from gluon.html import SELECT
# from gluon.html import SPAN
# from gluon.html import TABLE
# from gluon.html import TAG
# from gluon.html import TBODY
# from gluon.html import TD
# from gluon.html import TEXTAREA
# from gluon.html import TFOOT
# from gluon.html import TH
# from gluon.html import THEAD
# from gluon.html import TITLE
# from gluon.html import TR
# from gluon.html import TT
# from gluon.html import UL
# from gluon.html import URL
# from gluon.html import XHTML
# from gluon.html import XML
# from gluon.html import embed64
# from gluon.html import xmlescape
# from gluon.sql import SQLDB
# from gluon.sqlhtml import SQLFORM
# from gluon.sql import SQLField
# from gluon.sqlhtml import SQLTABLE
# from gluon.html import STYLE
# from gluon.http import redirect
# import gluon.languages.translator as T
# global cache; cache = gluon.cache.Cache()
# import gluon.compileapp.local_import_aux as local_import
# global request; request = gluon.globals.Request()
# global response; response = gluon.globals.Response()
# global session; session = gluon.globals.Session()
# global DAL; DAL = gluon.dal()
# global HTTP; HTTP = gluon.http()
# global LOAD; LOAD = gluon.compileapp.LoadFactory()