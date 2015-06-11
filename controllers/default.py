# -*- coding: utf-8 -*-
# this file is released under public domain and you can use without limitations

#########################################################################
## This is a sample controller
## - index is the default action of any application
## - user is required for authentication and authorization
## - download is for downloading files uploaded in the db (does streaming)
#########################################################################

import NW

def algorithm():
    
    d = ""

    form=FORM(DIV(LABEL('Sekwencja A:', XML("&nbsp;"), _for="seqA"), INPUT(_id='seqA', _name='seqA', _type='text', requires = IS_ALPHANUMERIC(error_message='Musisz podać AGCT'), _class='form-control', _placeholder = 'np. AAGCT'), _class='form-group'),XML("&nbsp;"),
              
              DIV(LABEL('Sekwencja B:', XML("&nbsp;"), _for="seqB"), INPUT(_id='seqB', _name='seqB', _type='text', requires = [IS_ALPHANUMERIC(error_message='Musisz podać AGCT')], _class='form-control', _placeholder = 'np. AAGCT'), _class='form-group'), XML("&nbsp;"),
              
              DIV(LABEL('Kara za przerwę:', XML("&nbsp;"), _for="break_penalty"), INPUT(_id='break_penalty', _name='break_penalty', requires=IS_INT_IN_RANGE(-10, 11,error_message='Musisz podać liczbę'), _class='form-control', _placeholder = 'np. -2'), _class='form-group'),XML("&nbsp;"),
              
              DIV(LABEL(B('Praca krokowa:'), XML("&nbsp;"), INPUT(_name='step',value=False,_type='checkbox',_class='checkbox')), _class = 'checkbox form-group' ),XML("&nbsp;"),
              
              INPUT(_type='submit', _class='btn btn-success', _value='Start'))
    
    form['_class']='form-inline'

    if form.accepts(request.vars, session, keepvalues=True):
        x = int(request.vars.break_penalty)
        d=NW.needlemanWunsch(request.vars.seqA,request.vars.seqB, x, 0)

    
    
    
    return dict(form=form, d=d)
    
    


def index():
    """
    example action using the internationalization operator T and flash
    rendered by views/default/index.html or views/generic.html

    if you need a simple wiki simply replace the two lines below with:
    return auth.wiki()
    """
    response.flash = T("Hello World")
    return dict(message=T('Welcome to web2py!'))


def user():
    """
    exposes:
    http://..../[app]/default/user/login
    http://..../[app]/default/user/logout
    http://..../[app]/default/user/register
    http://..../[app]/default/user/profile
    http://..../[app]/default/user/retrieve_password
    http://..../[app]/default/user/change_password
    http://..../[app]/default/user/manage_users (requires membership in
    use @auth.requires_login()
        @auth.requires_membership('group name')
        @auth.requires_permission('read','table name',record_id)
    to decorate functions that need access control
    """
    return dict(form=auth())


@cache.action()
def download():
    """
    allows downloading of uploaded files
    http://..../[app]/default/download/[filename]
    """
    return response.download(request, db)


def call():
    """
    exposes services. for example:
    http://..../[app]/default/call/jsonrpc
    decorate with @services.jsonrpc the functions to expose
    supports xml, json, xmlrpc, jsonrpc, amfrpc, rss, csv
    """
    return service()
