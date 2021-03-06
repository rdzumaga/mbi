# -*- coding: utf-8 -*-
# this file is released under public domain and you can use without limitations


import NW
import SW
import FinalFasta

def algorithm():

    d = ""


    form=FORM(DIV(LABEL('Sekwencja A:', _for="seqA", _class='col-sm-2 control-label'), DIV(INPUT(_id='seqA', _name='seqA', _type='text', _style='text-transform: uppercase', _class='form-control', _placeholder = 'np. AAGCT'), _class='col-md-4'), _class='form-group', _id='seqADiv'),

              DIV(LABEL('Sekwencja B:', _for="seqB", _class='col-sm-2 control-label'), DIV(INPUT(_id='seqB', _name='seqB', _type='text', _style='text-transform: uppercase',  _class='form-control', _placeholder = 'np. AAGCT'), _class='col-md-4'), _class='form-group',  _id='seqBDiv'),

              DIV(LABEL('Kara za przerwę:', _for="break_penalty", _class='col-sm-2 control-label'), DIV(INPUT(_id='break_penalty', _name='break_penalty', _class='form-control', _placeholder = 'np. -2'), _class='col-md-4'), _class='form-group', _id='seqCDiv'),

              DIV(DIV(DIV(LABEL(INPUT(_name='step', _id='checkBox' ,value=False,_type='checkbox'), 'Praca krokowa'), _class='checkbox'),_class='col-sm-offset-2 col-sm-4'), _class = 'form-group' ),

              INPUT(_id='iters', _name='iters', value='-1', _type='hidden')


              )

    form['_class']='form-horizontal'
    form['_id']='myForm'


    if form.accepts(request.vars, session, keepvalues=True):
        penalty=int(request.vars.break_penalty)
        iterx=int(request.vars.iters)
        d=NW.needlemanWunsch(iterx, request.vars.seqA,request.vars.seqB, penalty, "applications/mbi/modules/blosum.txt")


    return dict(form=form, d=d)

def smithWaterman():

	    d = ""


            form=FORM(DIV(LABEL('Sekwencja A:', _for="seqA", _class='col-sm-2 control-label'), DIV(INPUT(_id='seqA', _name='seqA', _type='text', _style='text-transform: uppercase',  _class='form-control', _placeholder = 'np. AAGCT'), _class='col-md-4'), _class='form-group', _id='seqADiv'),

            DIV(LABEL('Sekwencja B:', _for="seqB", _class='col-sm-2 control-label'), DIV(INPUT(_id='seqB', _name='seqB', _type='text', _style='text-transform: uppercase', _class='form-control', _placeholder = 'np. AAGCT'), _class='col-md-4'), _class='form-group',  _id='seqBDiv'),

            DIV(LABEL('Kara za przerwę:', _for="break_penalty", _class='col-sm-2 control-label'), DIV(INPUT(_id='break_penalty', _name='break_penalty', _class='form-control', _placeholder = 'np. -2'), _class='col-md-4'), _class='form-group', _id='seqCDiv'),

            DIV(LABEL('Zgodność:', _for="match", _class='col-sm-2 control-label'), DIV(INPUT(_id='match', _name='match', _class='form-control', _placeholder = 'np. -2'), _class='col-md-4'), _class='form-group', _id='seqDDiv'),

            DIV(LABEL('Niezgodność:', _for="mismatch", _class='col-sm-2 control-label'), DIV(INPUT(_id='mismatch', _name='mismatch', _class='form-control', _placeholder = 'np. -2'), _class='col-md-4'), _class='form-group', _id='seqEDiv'),

            DIV(DIV(DIV(LABEL(INPUT(_name='step', _id='checkBox' ,value=False,_type='checkbox'), 'Praca krokowa'), _class='checkbox'),_class='col-sm-offset-2 col-sm-4'), _class = 'form-group' ),

            INPUT(_id='iters', _name='iters', value='-1', _type='hidden'))

            form['_class']='form-horizontal'
            form['_id']='myForm'


            if form.accepts(request.vars, session, keepvalues=True):
                penalty=int(request.vars.break_penalty)
                iterx=int(request.vars.iters)
                matchint = int(request.vars.match)
                mismatchint = int(request.vars.mismatch)
                d=SW.SmithWaterman(iterx, request.vars.seqA,request.vars.seqB, penalty, matchint, mismatchint)



            return dict(form=form, d=d)


def fasta():


    d = ""

    db=FinalFasta.readDb("applications/mbi/modules/db.txt")

    form=FORM(DIV(LABEL('Sekwencja :', _for="seq", _class='col-sm-2 control-label'), DIV(INPUT(_id='seq', _name='seq', _type='text', _style='text-transform: uppercase', _class='form-control', _placeholder = 'np. AAGCT'), _class='col-md-4'), _class='form-group', _id='seqDiv'),

			  DIV(LABEL('Długość słowa: ', _for="length", _class='col-sm-2 control-label'), DIV(INPUT(_id='length', _name='length', _class='form-control', _placeholder = 'np. 2'), _class='col-md-4'), _class='form-group', _id='lengthDiv'),

			  DIV(LABEL('Próg C: ', _for="threshold_c", _class='col-sm-2 control-label'), DIV(INPUT(_id='threshold_c', _name='threshold_c', _class='form-control', _placeholder = 'np. 2'), _class='col-md-4'), _class='form-group', _id='threshold_cDiv'),

			  DIV(LABEL('Dopasowanie: ', _for="match", _class='col-sm-2 control-label'), DIV(INPUT(_id='match', _name='match', _class='form-control', _placeholder = 'np. 2'), _class='col-md-4'), _class='form-group', _id='matchDiv'),

              DIV(LABEL('Kara za przerwę:', _for="break_penalty", _class='col-sm-2 control-label'), DIV(INPUT(_id='break_penalty', _name='break_penalty', _class='form-control', _placeholder = 'np. -2'), _class='col-md-4'), _class='form-group', _id='break_penaltyDiv')






              )

    form['_class']='form-horizontal'
    form['_id']='myForm'


    if form.accepts(request.vars, session, keepvalues=True):
        penalty=int(request.vars.break_penalty)
        length=int(request.vars.length)
        breakpen=int(request.vars.break_penalty)
        treshold=int(request.vars.threshold_c)
        matchval=int(request.vars.match)
        d=FinalFasta.fasta(request.vars.seq,length,breakpen,treshold,matchval)


    return dict(form=form, d=d, db=db)

def index():

    redirect(URL(request.controller, 'algorithm'))

