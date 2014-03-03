from flask import render_template, flash, redirect, session, Response, url_for, g
from flask.ext.login import login_user, logout_user, current_user, login_required
from fpaste import app, db, lm
import forms
from models import FastaEntry
import models 
import ioroutines as ior

from hashlib import sha256
from base64 import urlsafe_b64encode
from datetime import datetime
from random import random

@lm.user_loader
def load_user(id):
    return models.User.query.get(int(id))

#@lm.needs_refresh_handler
#def refresh():
#    flash("A fresh login is needed to access this information")
#    return redirect(url_for('login'))
    
@app.before_request
def before_request():
    userLoggedIn = False
    if 'user_id' in session:
        g.user = models.User.query.get(int(session['user_id']))
        if g.user is not None:
            userLoggedIn = True
    g.userLoggedIn = userLoggedIn
            
@app.route('/login', methods=['GET', 'POST'])
def login():
    form = forms.UserLoginForm()
    if form.validate_on_submit():
        remember_me = form.remember_me.data
        user = models.User.query.filter_by(nickname=form.username.data.strip()).first()
        uploadedPass = form.userkey.data.strip()
        if user is not None and user.validate_password(uploadedPass):
            login_user(user, remember=remember_me)
            return redirect('/my_activity')
        else:
            form.username.errors.append("Invalid username or password")
    return render_template("login.html", form=form)    
        
@app.route('/logout')
@login_required
def logout():
    logout_user()
    flash("you have successfully logged out")
    return redirect("/login")

@app.route('/')
def redir():
    return redirect('/my_activity')

@app.route('/make_list', methods = ['GET', 'POST'])
@login_required
def make_list():
    form = forms.MakeListFromSelf()
    user = g.user
    
    if form.validate_on_submit():
        newFastaList = models.FastaList.query.filter_by(accessCode = form.fastaList.data).first()
        if newFastaList is None:
            idCode = urlsafe_b64encode(sha256(user.nickname+str(random())+str(random())).digest())[0:20]
            newFastaList = models.FastaList(accessCode = idCode, user_id = user.id)
            db.session.add(newFastaList)
            outputMessage = "Fasta list {} has been created".format(newFastaList.accessCode)
        else:
            outputMessage = "Fasta list {} has been modified".format(newFastaList.accessCode)
        #Process lists for added fastas and removed fastas
        for access in form.fastaAdd.data:
            entry = models.FastaEntry.query.filter_by(accessCode = access).first()
            if entry is not None:
                newFastaList.fastas.append(entry)
        for fasta in newFastaList.fastas:
            accessCode = fasta.accessCode
            if accessCode in form.fastaSubtract.data:
                delIndex = newFastaList.fastas.index(fasta)
                del(newFastaList.fastas[delIndex])
        form.fastaList.data = ""
        form.fastaAdd.data = []
        form.fastaSubtract.data = []
        db.session.commit()
        flash(outputMessage)               
    return render_template('make_list.html', form=form)


@app.route('/import_uniprot', methods = ['GET', 'POST'])
@login_required
def import_uniprot():
    form = forms.MakeListFromUniprot()
    errors = []
    if form.validate_on_submit() and g.userLoggedIn:
        user = g.user
        cutSignalSeq = form.cutSignalSeq.data
        ids = form.uniprotIds.data
        idCode = urlsafe_b64encode(sha256(user.nickname+str(random())+str(random())).digest())[0:20]
        newFastaList = models.FastaList(accessCode = idCode, user_id = user.id)
        db.session.add(newFastaList)
        for id in ids:
            uniprotOutput = ior.get_uniprot_protein_info(id, cutSignalSequence=cutSignalSeq) 
            if not uniprotOutput["success"]:
                errors.append(uniprotOutput["error"])
            else:
                meta = ">" + uniprotOutput["header"] + " " + uniprotOutput["url"]
                seq = uniprotOutput["seq"]
                newFasta = add_fasta_entry(meta, seq, type="protein")
                if newFasta["success"]:
                    fastaObj = newFasta["object"]
                    newFastaList.fastas.append(fastaObj)
                    db.session.add(fastaObj)
            db.session.commit()
    return render_template('make_list_uniprot.html', form=form)


    
    
@app.route('/my_activity')
@login_required
def my_activity():
    fastaDicts = []
    fastaListDicts =[]        
    if g.userLoggedIn:
        #pull individual fasta files
        fastas = models.FastaEntry.query.filter_by(user_id=g.user.id)
        for fasta in fastas:
            fastaDicts.append( {"id":fasta.accessCode , "accession":fasta.accession} )
        fastaLists = models.FastaList.query.filter_by(user_id=g.user.id)
        #pull all fasta lists
        for fastaL in fastaLists:
            fastaListDicts.append( {'id':fastaL.accessCode } )
    return render_template( "user_content.html",
                            fastaList = fastaDicts,
                            fastaListList = fastaListDicts)    



@app.route('/paste', methods = ['GET', 'POST'])
@login_required
def single_paste():
    
    form = forms.SingleFastaPaste()

    if form.validate_on_submit() and g.userLoggedIn:
        metaLine = form.description.data
        seq = form.sequence.data
        addition = add_fasta_entry(metaLine, seq, type="protein") #returns {"success":bool, "error":str, "code":str}
        if addition["success"]:
            code = addition["code"]
            flash("new fasta added with access code = {}".format(code))
        else:
            form.sequence.errors.append(addition["error"])
    return render_template('paste.html',
                           form = form)

@app.route("/fasta/<fasta_id>", methods = ['GET', 'POST'])
def fasta_details(fasta_id):
    form = forms.DeleteHidden()
    fastaDetails = {}
    fasta = models.FastaEntry.query.filter_by(accessCode = fasta_id).first()
    if fasta is None:
        flash("no fata matches this id")
        return redirect("/my_activity")
    
    if form.validate_on_submit():
        flash("fasta entry <{}> deleted".format(fasta_id))
        db.session.delete(fasta)
        db.session.commit()
        return redirect("/my_activity")
    else:
        fastaDetails = {}
        fastaDetails["id"] = fasta.accessCode
        fastaDetails["added"] = str(fasta.added)
        fastaDetails["meta"] = ">" + fasta.accession + " " + fasta.metaData
        return render_template('fasta_details.html',
                               fastaDetails = fastaDetails,
                               form = form)
                               
@app.route("/fastalist/<fastalist_id>", methods = ['GET', 'POST'])
def fasta_list_details(fastalist_id):
    form = forms.DeleteHidden()
    fastaDetails = {}
    fastaList = models.FastaList.query.filter_by(accessCode = fastalist_id).first()
    if fastaList is None:
        flash("no fata matches this id")
        return redirect("/my_activity")
    
    if form.validate_on_submit():
        flash("fasta entry <{}> deleted".format(fastalist_id))
        db.session.delete(fastaList)
        db.session.commit()
        return redirect("/my_activity")
    else:
        fastaDetails = {}
        fastaDetails["accessCode"] = fastaList.accessCode
        #fastaDetails["added"] = str(fasta.added)
        fastaDetails["added"] = "DATE GOES HERE"
        fastaDetails["fastas"] = fastaList.fastas
        return render_template('fasta_list_details.html',
                               fastaDetails = fastaDetails,
                               form = form)

@app.route("/fastalist/<fastalist_id>.fasta", methods = ['GET', 'POST'])
def render_fasta_list(fastalist_id):
    fastaList = models.FastaList.query.filter_by(accessCode = fastalist_id).first()
    if fastaList is None:
        fastaReturn = "notfound"
    else:
        fastaReturn = ""
        for fasta in fastaList.fastas:
            fastaReturn += ">" + fasta.accession + " " + fasta.metaData + "\n"
            counter = 0
            for char in fasta.sequence:
                counter += 1
                fastaReturn += char
                if counter % 80 == 0:
                    counter = 0
                    fastaReturn += "\n"
            fastaReturn += "\n"
    return Response(fastaReturn, content_type="text/plain;charset=UTF-8")  
    
    
@app.route("/fasta/<fasta_id>.fasta")
def render_fasta_file(fasta_id):
    fasta = models.FastaEntry.query.filter_by(accessCode = fasta_id).first()
    if fasta is None:
        fastaReturn = ""
    else:
        fastaReturn = ">" + fasta.accession + " " + fasta.metaData + "\n"
        counter = 0
        for char in fasta.sequence:
            counter += 1
            fastaReturn += char
            if counter % 80 == 0:
                counter = 0
                fastaReturn += "\n"
    return Response(fastaReturn, content_type="text/plain;charset=UTF-8")



def add_fasta_entry(meta, seq, type="protein"): #returns {"success":bool, "error":str, "code":str}
    #meta data line reading and grab type/sequence
    sequence = ''.join([aa.strip() for aa in seq])
    
    # check for existing seq
    binHash = sha256(meta+sequence).digest()
    duplicated, lencut = (True, 15) 
    while duplicated and lencut < 20:
        b64hash = urlsafe_b64encode(binHash)[0:lencut]
        fastaInDb = models.FastaEntry.query.filter_by(accessCode=b64hash).first()
        lencut += 1
        if fastaInDb is None:
            duplicated = False
        elif fastaInDb.sequence == sequence:
            return {"success":True,
                    "error":"",
                    "code":b64hash,
                    "object":fastaInDb}
    
    if duplicated:
        return {"success":False, "error":"duplicated", "code":b64hash}
    else:
        newFasta = FastaEntry()
        #append scrubbed sequence and meta line with accession
        newFasta.sequence = sequence
        newFasta.read_in_meta_line(meta)
        #getcode and asign
        newFasta.accessCode = b64hash
        #added datetime
        newFasta.added = datetime.utcnow()
        #add user meta info
        newFasta.user_id = g.user.id
        db.session.add(newFasta)
        db.session.commit()
        return {"object":newFasta, "success":True, "error":'', "code":b64hash}
        #flash("new fasta added with access code = {}".format(b64hash))
