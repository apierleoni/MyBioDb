__author__ = 'pierleonia'

def index():
    content_body = DIV()
    content_body.append(H2("DB Index"))
    content_body.append(UL(LI(A("rebuild index (long process)", _href= URL(r= request, f = 'rebuild_index')))))

    content_body.append(H2("Load Test Data"))
    content_body.append(UL(LI(A("load whole uniprot (long process)", _href= URL(r= request, f = 'local_import_huge')))))

    return dict(content_body= content_body)

def rebuild_index():
    biodb_index.rebuild()
    return "Done"


def local_import_huge():
    '''used to quickly load a big db for testing'''
    import traceback,time
    created_seqrec = []
    errors =[]
    i= 0
    from Bio import SeqIO
    for seqrec in SeqIO.parse('/Users/pierleonia/Downloads/uniprot/uniprot_sprot.xml', 'uniprot-xml'):
        i+=1
        print 'uploading entry ', i,  #used for debug remove in production code
        if 437590<i< settings.max_entry_load:
            start_time = time.time()
            try:
                created_seqrec_id = biodb_handler.load_seqrecord(seqrec, db = 'UniProt')
                created_seqrec.append((seqrec.id, created_seqrec_id))
                biodb.commit()
                biodb_index.add_bioentry_id_to_index(created_seqrec_id)
                print ' -> SUCCESS in ', int(round(time.time() - start_time, 3)*1000), 'ms'
            except Exception, e:
                print ' -> ERROR:', e
                errors.append((seqrec.id, e))
                traceback.print_exc()


        else:
            print ' -> SKIPPED'

    print  '%i entries correctly loaded on db'%len(created_seqrec),
    print '| %i errors while loading'%len(errors)


def view_uniprot_entries():
    returnul = UL()
    for row in biodb(biodb.bioentry.id>0).select():
        returnul.append(LI(A(row.name,_href = URL(r=request, f = view, vars = dict(bioentry_id = row.id )))))

    return dict(returnul = returnul)