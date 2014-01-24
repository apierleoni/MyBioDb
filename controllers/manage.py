__author__ = 'pierleonia'

def index():
    content_body = DIV()
    content_body.append(H2("DB Index"))
    content_body.append(UL(LI(A("rebuild index (long process)", _href= URL(r= request, f = 'rebuild_index'))))
                        )

    return dict(content_body= content_body)

def rebuild_index():
    biodb_index.rebuild()
    return "Done"