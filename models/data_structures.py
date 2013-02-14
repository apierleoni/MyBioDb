# coding: utf8


'''This file contains all the DAta archetype structures needed for the portal '''



class UnitView():
    '''represent an article object in the main content part
    
        <div class="unitview"> 
            <div class="unitview-header"> 
                <div class="unitview-time"> 
                    <div class="unitview-year">2010</div> 
                    <div class="unitview-date">16<span>apr</span></div> 
                </div> 
                <h1 class="unitview-title">Sample Post 1</h1> 
                <div class=class="unitview-counter">38</div> 
            </header> 
            <p>content text here</p> 
            <footer> 
                <em>Created by:</em> <strong>Author Name</strong> 
                <span class="unitview-tags"><em>Tags:</em> <a class=tags href=#>cool</a><a class=tags href=#>modern</a></span> 
            </footer>
        </div> 
    
    '''
    
    def __init__(self,
                 title = '',
                 counts = '',
                 year = '',
                 day = '',
                 month = '',
                 content = '',
                 author = '',
                 tags = '',
                 button = '',
                 footer = '',
                 watermark = False,
                 script = '',
                 _class = '',
                 _id = '', ):
        
        try:
            self.title = title.xml()
        except AttributeError:
            self.title = title
        try:
            self.counts = counts.xml()
        except AttributeError:
            self.counts = str(counts)
        try:
            self.year = year.xml()
        except AttributeError:
            self.year = str(year)
        try:
            self.day = day.xml()
        except AttributeError:
            self.day = str(day)
        try:
            try:
                month = int(month)
                dmonth = {1:'jan',
                          2:'feb',
                          3:'mar',
                          4:'apr',
                          5:'may',
                          6:'jun',
                          7:'jul',
                          8:'aug',
                          9:'sep',
                          10:'oct',
                          11:'nov',
                          12:'dec',}
                self.month = dmonth[month]
            except:
                
                self.month = month.xml().lower()[:3]
        except AttributeError:
            self.month = str(month).lower()[:3]
        try:
            self.content = content.xml()
        except AttributeError:
            self.content = content
        try:
            self.author = author.xml()
        except AttributeError:
            self.author = author
        try:
            self.footer = footer.xml()
        except AttributeError:
            self.footer = footer
        try:
            self.tags = tags.xml()
        except AttributeError:
            if isinstance(tags, list): #sting list of tags
                self.tags =''
                for tag in tags:
                    self.tags+= '''<a class="unitview-tag" href=%s>%s</a> '''%(URL(c= 'news', f='tags', args= tag), tag)
            
            else:
                self.tags = tags
        try:
            self.button = button.xml()
        except AttributeError:
            self.button = button  #example: <a class=button href=#>Continue Reading</a> 
        
        
        self.script = SCRIPT(script).xml() #an empty script must be initialized to make livequery work
        
        self.watermark = watermark
        self._id = _id
        self._class = _class
        
        
    def xml(self):
        
        #build date
        if self.year or self.day or self.month:
            timediv = '''                          <div class="unitview-time"> 
                          <div class="unitview-year">%s</div> 
                          <div class="unitview-date">%s<span>%s</span></div>
                          </div> '''%(self.year, self.day, self.month)
        else:
            timediv = ''''''
        
        if self.author or self.tags or self.button or self.footer: 
            author_line = tags_line = ''
            if self.author:
                author_line = '''<em>Created by:</em> <strong>%s</strong> '''%self.author
            if self.tags:
                tags_line = '''<span class="unitview-tags"><em>Tags:</em> %s </span> '''%self.tags
            footerdiv='''                <div class="unitview-footer"> 
                   %s
                   %s 
                   %s
                   %s
                </div>''' %(author_line,tags_line,self.button, self.footer)
        else:
            footerdiv='<div class="unitview-footer"></div>'
        
        if self.watermark:
            watermark='''<div class="unitview-watermark" style="margin-top:-50px;float:right; z-index:-10; opacity:0.3;"><img src="%s" alt="primm" style="height: 50px;"/></div>'''%URL(r=request, c='static', f='logo_primm.png')
        else:
            watermark='' 
        
        fullxml = '''        <div id="%s" class="unitview row %s">
            <div class="unitview-header" onclick = "void(0)"> 
                %s
                <h2 class="unitview-title">%s</h2> 
                <div class="unitview-counter">%s</div> 
            </div> 
            <div class="unitview-content">%s</div> 
            %s
            %s
        </div> %s ''' %(self._id, self._class, '', self.title, self.counts, self.content,footerdiv,watermark, self.script)#timediv, self.title, self.counts, self.content,footerdiv,watermark, self.script)
        
        
        return fullxml
    
    def toxml(self):#compatibility
        return self.xml()
    

        
        
        

