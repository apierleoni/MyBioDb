response.title = settings.title
response.subtitle = settings.subtitle
response.meta.author = '%(author)s <%(author_email)s>' % settings
response.meta.keywords = settings.keywords
response.meta.description = settings.description
response.menu = [
('Index',URL('default','index')==URL(),URL('default','index'),[]),
]
#left_sidebar_enabled = True