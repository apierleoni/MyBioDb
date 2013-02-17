'''
Reusable portions of code, to be used by all controllers
'''


def draw_sequence(sequence, mode = 'simple', alphabet = None):
        
    if mode == 'protparams':
        returndiv = DIV()
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        seq_div=DIV(_style='font-family:monospace',_class='raw-sequence')
        spacer=len(str(len(sequence)))+1
        for i,pos in enumerate(sequence):
            if i==0:
                seq_div.append(XML((str(i+1)+' ').rjust(spacer).replace(' ','&nbsp;')))
            if i%10==0 and i!=0:
                seq_div.append(' ')
            if i%60==0 and i!=0:
                seq_div.append(XML((str(i)).ljust(spacer).replace(' ','&nbsp;')))
                seq_div.append(BR())
                seq_div.append(XML((str(i+1)+' ').rjust(spacer).replace(' ','&nbsp;')))
            seq_div.append(SPAN(pos,_class='seq-position',_title = i+1))
        returndiv.append(seq_div)
        returndiv.append(H3('Protein Parameters'))
        params_table = TABLE(_style= "width:200px;")
        
        protpar=ProteinAnalysis(sequence)
        params_table.append(TR(SPAN('Length:',_class = 'line-header'), '%i aa'%len(sequence)))
        try:
            params_table.append(TR(SPAN('MW:',_class = 'line-header'), '%i KDa'%round(protpar.molecular_weight()/1000,0)))
        except KeyError:
            pass
        try:
            params_table.append(TR(SPAN('pI:',_class = 'line-header'), '%1.2f'%protpar.isoelectric_point()))
        except KeyError:
            pass
        returndiv.append(params_table)
        return returndiv
        
    if mode == 'simple':
        seq_div=DIV(_style='font-family:monospace',_class='raw-sequence')
        spacer=len(str(len(sequence)))+1
        for i,pos in enumerate(sequence):
            if i==0:
                seq_div.append(XML((str(i+1)+' ').rjust(spacer).replace(' ','&nbsp;')))
            if i%10==0 and i!=0:
                seq_div.append(' ')
            if i%60==0 and i!=0:
                seq_div.append(XML((str(i)).ljust(spacer).replace(' ','&nbsp;')))
                seq_div.append(BR())
                seq_div.append(XML((str(i+1)+' ').rjust(spacer).replace(' ','&nbsp;')))
            seq_div.append(SPAN(pos,_class='seq-position', _title = i+1))
        return seq_div
      

def draw_graphic_view(start, end, length):
    min_width = 2
    total_width = 100
    margin_left= round(float(start)/length*total_width,0)
    width = round(float(end-start)/length*total_width,0)
    if width <= min_width :
        width = min_width
    whole = DIV(_class = 'graph-feat-ref', _style = 'width: %i%%; '%total_width)
    feat = DIV( _class = 'graph-feat-drawn', _style = 'margin-left: %i%%; width: %i%%; '%(margin_left,width))
    returndiv= DIV(whole, feat)
    
    return returndiv
        
def draw_single_feature(feat,length, feat_id):
    '''takes a biosqlfeature obj and fill a TR for output in a table'''

    
    if feat.location.start.position == feat.location.end.position:
        feat.location.end.position+=1
        
    description = feat.qualifiers.get('description','')
    id = feat.qualifiers.get('id','')
    
    qualifiers = SPAN()
    all_qualifiers= feat.qualifiers
    for i in ['description', 'id', 'type']:
        if i in all_qualifiers:
            all_qualifiers.pop(i)
    len_k = len(all_qualifiers)
    for i, kv in enumerate(all_qualifiers.items()):
        k,v = kv
        qualifiers.append(SPAN(k,': ', _style = 'font-weight: bold; margin-left: 5px;'))
        qvalue = ''
        if isinstance(v,(list, tuple)):
#            for qrank,value in enumerate(v):
#                qualifiers.append(LI(SPAN(v,_class = 'edit-field-feature-qualifier-value', **{'_data-pk': feat_id,'_data-rank': qrank, '_data-key': k})))
            qvalue='|'.join(v)
        qualifiers.append(SPAN(qvalue,_class = 'edit-field-feature-qualifier-value', **{'_data-pk': feat_id, '_data-key': k}))
#    if i != len_k-1:         #TODO: fix this
#            qualifiers.append('; ')

    graphic_view = DIV()
    if length:
        graphic_view = draw_graphic_view(feat.location.start.position, feat.location.end.position, length),#graphical view
        
    return TR(SPAN(feat.type, _class = 'edit-field-feature-type', _style = 'font-weight: bold', **{'_data-pk': feat_id} ),
              TD(SPAN(str(feat.location.start.position+1),_class = 'edit-field-feature-start', **{'_data-type': 'number','_data-pk': feat_id})
                 +'-'+
                 SPAN(str(feat.location.end.position), _class = 'edit-field-feature-end', **{'_data-type': 'number','_data-pk': feat_id})),

              feat.location.end.position - feat.location.start.position,
              SPAN(description,_class = 'edit-field-feature-description',**{'_data-pk': feat_id}),
              TD(graphic_view, **{'_class': 'view-feature','_rel':"clickover", '_data-placement':"bottom",
                                   '_data-original-title':"Feature type: "+feat.type,
                                   '_data-start':feat.location.start.position+1,
                                   '_data-end':feat.location.end.position}),
              qualifiers,
              SPAN(id,_class = 'edit-field-feature-identifier', **{'_data-pk': feat_id}),
              **{'_data-pk': feat_id})


def draw_feature_on_sequence_utility(sequence, start, end):
    return_div=DIV() 
    seq_div=DIV(_style='font-family:monospace',_class='raw-sequence')
    spacer=len(str(len(sequence)))+1
    if not isinstance(start,(list, tuple)) :
        start = [start]
    if not isinstance(end,(list, tuple)) :
        end = [end]
    for i,pos in enumerate(sequence):
        styleclass='not-active-seq-position'
        for j in xrange(len(start)):
            if start[j]<i+1<=end[j]:
                styleclass='active-seq-position'
        if i==0:
            seq_div.append(XML((str(i+1)+' ').rjust(spacer).replace(' ','&nbsp;')))
        if i%10==0 and i!=0 and i%60!=0:
            seq_div.append(SPAN(' ', _class=styleclass))
        if i%60==0 and i!=0:
            seq_div.append(XML((' '+str(i)).ljust(spacer).replace(' ','&nbsp;')))
            seq_div.append(BR())
            seq_div.append(XML((str(i+1)+' ').rjust(spacer).replace(' ','&nbsp;')))
        seq_div.append(SPAN(pos,_class=styleclass, _title = i+1))
    return_div.append(seq_div) 
    
    return return_div

    
def biosql_data_beautify(value):
    if isinstance(value,str):
        return value
    elif isinstance(value,list):
        if len(value) == 1:
            if isinstance(value[0],str):
                return value[0]
            else:
                return BEAUTIFY(value[0])
        listul = UL()
        for v in value:
            if isinstance(v,str):
                listul.append(LI(v))
            else:
                listul.append(LI(BEAUTIFY(v)))
        return listul
    else:#should not happen on data caming from a biosql db
        return BEAUTIFY(value)
        

def create_dbxref_link(dbname, accession, dbxref_id):
    Dlinks ={'ccds': 'http://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi?REQUEST=ALLFIELDS&DATA=%s',
             'cdd': 'http://www.ncbi.nlm.nih.gov/sites/entrez/query.fcgi?db=cdd&term=%s',
             'embl': 'http://www.ebi.ac.uk/cgi-bin/sva/sva.pl/?search=Go!&query=%s',
             'ensembl': 'http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s',
             'ensg': 'http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s',
             'ensp': 'http://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?p=%s',
             'enst': 'http://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=%s',
             'entrezgene': 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&term=%s',
             'geneid': 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&term=%s',
             'germonline': 'http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s',
             'gi': 'http://www.ncbi.nlm.nih.gov/protein/%s',
             'go': 'http://amigo.geneontology.org/cgi-bin/amigo/term_details?term=%s',
             'hgnc': 'http://www.genenames.org/data/hgnc_data.php?hgnc_id=%s',
             'hgnc_id': 'http://www.genenames.org/data/hgnc_data.php?hgnc_id=%s',
             'hgnc_name': 'http://www.genenames.org/data/hgnc_data.php?match=%s',
             'hprd': 'http://www.hprd.org/resultsQuery?multiplefound=&prot_name=&external=Ref_seq&accession_id=&gene_symbol=&chromo_locus=&function=&ptm_type=&localization=&domain=&motif=&expression=&prot_start=&prot_end=&limit=0&mole_start=&mole_end=&disease=&query_submit=Search&hprd=%s',
             'interpro': 'http://www.ebi.ac.uk/interpro/ISearch?query=%s',
             'ipi': 'http://www.ebi.ac.uk/cgi-bin/dbfetch?db=IPI&id=%s',
             'kegg': 'http://www.genome.jp/dbget-bin/www_bget?%s',
             'merops': 'http://merops.sanger.ac.uk/cgi-bin/make_frame_file?id=%s',
             'mim': 'http://omim.org/entry/%s',
             'mim_gene': 'http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=%s',
             'mim_morbid': 'http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=%s',
             'otthumg': 'http://vega.sanger.ac.uk/Homo_sapiens/geneview?gene=%s',
             'ottp': 'http://vega.sanger.ac.uk/Homo_sapiens/protview?peptide=%s',
             'ottt': 'http://vega.sanger.ac.uk/Homo_sapiens/transview?transcript=%s',
             'panther': 'http://www.pantherdb.org/panther/family.do?clsAccession=%s',
             'pdb': 'http://www.rcsb.org/pdb/explore/explore.do?structureId=%s',
             'pdbsum': 'http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?template=main.html&EBI=TRUE&pdbcode=%s',
             'peptideatlas': 'https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/Search?action=GO&build_type_name=Any&all_fields=on&search_key=%s',
             'pfam': 'http://pfam.sanger.ac.uk//family/%s',
             'pharmgkb': 'http://www.pharmgkb.org/do/serve?objId=%s',
             'pir': 'http://pir.georgetown.edu/cgi-bin/textsearch.pl?submit.x=0&submit.y=0&field0=ALL&search=1&query0=%s',
             'pirsf': 'http://pir.georgetown.edu/cgi-bin/ipcSF?id=%s',
             'prints': 'http://www.bioinf.manchester.ac.uk/cgi-bin/dbbrowser/sprint/searchprintss.cgi?display_opts=Prints&category=None&queryform=false&prints_accn=%s',
             'profile': 'http://www.expasy.ch/prosite/%s',
             'prosite': 'http://www.expasy.ch/prosite/%s',
             'prosite_pattern': 'http://www.expasy.ch/prosite/%s',
             'protein_id': 'http://www.ncbi.nlm.nih.gov/protein/%s',
             'pubmed': 'http://www.ncbi.nlm.nih.gov/sites/entrez/%s',
             'refseq': 'http://www.ncbi.nlm.nih.gov/protein/%s',
             'refseq_dna': 'http://www.ncbi.nlm.nih.gov/nuccore/%s',
             'refseq_peptide': 'http://www.ncbi.nlm.nih.gov/protein/%s',
             'shares_cds_with_ottt': 'http://vega.sanger.ac.uk/Homo_sapiens/transview?transcript=%s',
             'smart': 'http://smart.embl-heidelberg.de/smart/do_annotation.pl?BLAST=DUMMY&DOMAIN=%s',
             'superfamily': 'http://supfam.cs.bris.ac.uk/SUPERFAMILY/cgi-bin/scop.cgi?sunid=%s',
             'tigrfam': 'http://cmr.jcvi.org/cgi-bin/CMR/HmmReport.cgi?hmm_acc=%s',
             'ucsc': 'http://genome.ucsc.edu/cgi-bin/hgGene?hgg_prot=Q9HAU5&hgg_chrom=chr10&hgg_start=12002026&hgg_end=12124814&hgg_type=knownGene&db=hg18&hgg_gene=%s',
             'unigene': 'http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=%s&CID=%s',
             'uniprot': 'http://www.uniprot.org/uniprot/%s',
             'uniprot/sptrembl': 'http://www.uniprot.org/uniprot/%s',
             'uniprot/swissprot': 'http://www.uniprot.org/uniprot/%s',
             'swiss-prot': 'http://www.uniprot.org/uniprot/%s',
             'uniprot/varsplic': 'http://www.uniprot.org/uniprot/%s',
             'unists': 'http://www.ncbi.nlm.nih.gov/genome/sts/sts.cgi?uid=%s',
             'mint': 'http://mint.bio.uniroma2.it/mint/search/search.do?queryType=protein&interactorAc=%s',
             'intact': 'http://www.ebi.ac.uk/intact/pages/interactions/interactions.xhtml?query=%s*',
             'string': 'http://string-db.org/newstring_cgi/show_network_section.pl?identifier=%s',
             'phosphosite':'http://www.phosphosite.org/uniprotAccAction.do?id=%s',
             'pride':'http://www.ebi.ac.uk/pride/searchSummary.do?queryTypeSelected=identification%%20accession%%20number&identificationAccessionNumber=%s',
             'hpa':'http://www.proteinatlas.org/tissue_profile.php?antibody_id=%s',
             'genecards':'http://www.genecards.org/cgi-bin/carddisp.pl?gc_id=%s',
             'nextprot':'http://www.nextprot.org/db/entry/%s',
             'arrayexpress':'http://www.ebi.ac.uk/arrayexpress/genes/uniprot/%s',
             'bgee':'http://bgee.unil.ch/bgee/bgee?uniprot_id=%s',
             'nextbio':'http://www.nextbio.com/b/home/home.nb?id=%s&type=feature',
             'dip':'http://dip.doe-mbi.ucla.edu/dip/Browse.cgi?ID=%s',
             'doi':'http://dx.doi.org/%s',
            }
    dbname = dbname.lower()
    if dbname in Dlinks:
        if dbname == 'unigene':
			cid = accession.split('.')[-1]
			org = accession.split('.')[0]
			return A(accession, _href = Dlinks[dbname]%(org,cid), _target = 'blank')
		#if dbname == 'ensemble':
		#	dbxrefdb = key[:4]
		#if dbname == 'hgnc':
		#	key = key.split(':')[-1]
		#	if key.isdigit():  dbxrefdb = 'HGNC_ID'
		#	else:  dbxrefdb = 'HGNC_NAME'
		#if dbname == 'hprd':
		#	key = key.split('HPRD_')[-1]
		#if dbname == 'ipi':
		#	key = key.split('.')[0]
        return A(accession, _href = Dlinks[dbname]%accession, _target = 'blank', _class = "edit-field-dbxref-accession", **{'_data-pk':dbxref_id})
    else:
        return SPAN(accession,_class = "edit-field-dbxref-accession",**{'_data-pk':dbxref_id})
        
def draw_single_dbxref(dbxref_string, dbxref_id):
    '''takes a biosqlfeature obj and fill a TR for output in a table'''

    decomposed = dbxref_string.split(':')
    dbname = decomposed[0]
    accession = ':'.join(decomposed[1:])
    return TR(SPAN(dbname, _class = "edit-field-dbxref-dbname",**{'_data-pk':dbxref_id}),
              create_dbxref_link(dbname,accession,dbxref_id), **{'_data-pk':dbxref_id})
              
              
              
def draw_single_reference(reference):
    '''takes as input a biopython reference'''
    'comment', 'consrtm', 'location', 

    refdiv = DIV(_class = 'jref')
    refdiv.append(H3(reference.title,_class = 'jref-title'))
    journal_P = P(reference.journal, _class = 'jref-journal')
    refdiv.append(journal_P)
    if reference.pubmed_id:
        refdiv.append(A('PubMed',
                                _href = 'http://www.ncbi.nlm.nih.gov/pubmed/%s'%reference.pubmed_id,
                                _target = 'blank',
                              _class = 'btn btn-primary btn-small'))
        
    authors = reference.authors.split(', ')
    author_P = P(_class = 'jref-authors')
    nauthors = len(authors)
    for i, a in enumerate(authors):
        author_P.append(A(a, 
                          _href = 'http://www.ncbi.nlm.nih.gov/pubmed?term=%%22%s%%22%%5BAuthor%%5D'%a.replace('.',''), 
                          _class = 'jref-author',
                          _target = 'blank'))
        if i< (nauthors-2):
            author_P.append(', ')
        elif i == (nauthors-2):
            author_P.append(' and ')
    refdiv.append(author_P)
    if reference.pubmed_id:
        refdiv.append(DIV(A('show me more...' ,
                            _class = 'load-pubmed-data btn btn-mini pull-right',
                            _onclick = 'void(0)'),
                          _id = 'pubmed-data-%s'%reference.pubmed_id))
    if reference.comment:
        refdiv.append(P(SPAN('Comment: ', _style = 'font-weight: bold'), reference.comment ,_class = 'jref-comment'))
    if reference.consrtm:
        refdiv.append(P(SPAN('Consortium: ', _style = 'font-weight: bold'), reference.consrtm ,_class = 'jref-consrtm'))
    if reference.location:
        try:
            refdiv.append(P(SPAN('This reference is annotated on sequence position(s): ', _style = 'font-weight: bold'), 
                       str(reference.location[0].start)+'-'+str(reference.location[0].end), _class = 'jref-location'))
        except:#skip malformed locations
            pass
    
    return refdiv
    
    
def fetch_pubmed_data(pmid):

    from Bio import Medline,Entrez
    
    try:
        ncbiemail= settings.author_email
    except:
        try:
            ncbiemail= settings.author_email
        except:
            raise Exception('Please set an email to use ncbi services')
    
    Entrez.email = ncbiemail
    Entrez.tool = 'web2biopy'

    try:
        entrez_response=Medline.parse( Entrez.efetch(db="pubmed", id=pmid, rettype="medline",retmode="text",)).next()
        if not entrez_response.has_key('PMID'):
             response.flash='pubmed ID error'
        else:
            return entrez_response
    except IOError:
        session.flash='Remote service not available, please try again.'

       
    return


def return_feature_type_select(biosqldb):
    feat_types = SELECT(sorted([row.name for row in biodb((biosqldb.term.ontology_id == biosqldb.ontology.ontology_id) & \
                                                  (biosqldb.ontology.name == 'SeqFeature Keys' )).select(biosqldb.term.name)]),
                        _name = 'type', _id = 'featuretype-select',)
    return feat_types

    
    