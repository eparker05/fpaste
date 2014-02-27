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
def make_list():
    form = forms.MakeListFromSelf()
    validated = form.validate_on_submit()
        
    #user login????
    loggedIn = True
    if validated:
        user = g.user
            
    if validated and loggedIn:
        pass
        
    return render_template('make_list.html', form=form)





@app.route('/import_uniprot', methods = ['GET', 'POST'])
@login_required
def import_uniprot():
    form = forms.MakeListFromUniprot()
    errors = []
    if form.validate_on_submit() and g.userLoggedIn:
        user = g.user
        cutSignalSeq = form.cutSignalSeq.data            
        ids = form.uniprotIds
        for id in ids:
            out = ior.get_uniprot_protein_info(id) 
            if not out["success"]:
                errors.append({"link":out["meta"], "name":out["header"], "error":out["error"]})
            else:
                pass                
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
            fastaListDicts.append( {'id':fastL.accessCode } )
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
    
    if form.validate_on_submit():
        return "<p>"+form.deleteThis.data+"</p>"
    else:
        fasta = models.FastaEntry.query.filter_by(accessCode = fasta_id).first()
        if fasta is None:
            flash("no fata matches this id")
            return redirect("/my_activity")
        else:
            fastaDetails = {}
            fastaDetails["id"] = fasta.accessCode
            fastaDetails["added"] = str(fasta.added)
            fastaDetails["meta"] = ">" + fasta.accession + " " + fasta.metaData
            return render_template('fasta_details.html',
                                   fastaDetails = fastaDetails,
                                   form = form)
    
    
    
@app.route("/fasta/<fasta_id>.txt")
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
    metaLine = meta
    sequence = ''.join([aa.strip() for aa in seq])
    
    # check for existing seq
    binHash = sha256(metaLine+sequence).digest()
    duplicated, lencut = (True, 15) 
    while duplicated and lencut < 20:
        b64hash = urlsafe_b64encode(binHash)[0:lencut]
        fastaInDb = models.FastaEntry.query.filter_by(accessCode=b64hash).first()
        if fastaInDb is None:
            duplicated = False
        elif fastaInDb.sequence == sequence:
            return {"success":False,
                    "error":"this sequence exists as record {}".format(b64hash),
                    "code":b64hash}
        else:
            lencut += 1
    
    if duplicated:
        return {"success":False, "error":"duplicated", "code":b64hash}
    else:
        newFasta = FastaEntry()
        #append scrubbed sequence and meta line with accession
        newFasta.sequence = sequence
        newFasta.read_in_meta_line(metaLine)
        #getcode and asign
        newFasta.accessCode = b64hash
        #added datetime
        newFasta.added = datetime.utcnow()
        #add user meta info
        newFasta.user_id = g.user.id
        db.session.add(newFasta)
        db.session.commit()
        return {"success":True, "error":'', "code":b64hash}
        #flash("new fasta added with access code = {}".format(b64hash))
