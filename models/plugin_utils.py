def url(f,args=request.args,vars={}):
    return URL(r=request,f=f,args=args,vars=vars)

def goto(f,args=request.args,vars={},message='error'):
    session.flash=message
    redirect(url(f,args=args,vars=vars))

def error():
    goto('error')

def get(table, i=0, message='error'):
    try:
        id = int(request.args(i))
    except ValueError:
        goto('error',message=message)
    return table[id] or goto('error',message=message)
