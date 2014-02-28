
from flask.ext.wtf import Form
from wtforms import Field, TextField, TextAreaField, BooleanField, RadioField, HiddenField, widgets
from wtforms.validators import Required, Length, ValidationError, Regexp, EqualTo

import re

"""custom fields used for fasta-paste forms """

class CsvTextAreaField(Field):
    widget = widgets.TextArea()
    
    def _value(self):
        if self.data:
            return u', '.join(self.data)
        else:
            return u''
            
    def process_formdata(self, valuelist):
        if valuelist:
            self.data = [v.strip() for v in valuelist[0].split(',')]
        else:
            self.data = []

class CsvTextField(CsvTextAreaField):
    widget = widgets.TextInput()

#""" actual forms """
    
class ValidateFastaSeq(object):
    """
    checks that a sequence is valid
    """
    def __init__(self, typefield=None, message=None):
        self.fieldname = typefield
        if message is None:
            self.message = "the given sequence does not conform to standard residue characters"
        else:
            self.message = message
        
    def __call__(self, form, field):
        try:
            fastaTypeField = form[self.fieldname]
        except KeyError:
            raise ValidationError("invalid field name '{}'".format(self.fieldname))
            
        type = fastaTypeField.data
        seqRaw = field.data.upper()
        if not self.run_validation(seqRaw, type):
            raise ValidationError(self.message)
        
    def run_validation(self, seqRaw, type):
        """
        removes formatting characters then checks if other 
        characters are valid.
        returns a bool
        """
        sequence0 = ''.join([aa.strip() for aa in seqRaw]) 
        if type == "protein":
            sequence = [aa.strip() for aa in sequence0 if aa.isalpha() or aa in "*-"]
            sequence = ''.join(sequence)
        else:
            sequence = [na.strip() for na in sequence0 if na in "ATKMBVCNSWDGUYRH-"]
            sequence = ''.join(sequence)
        
        if sequence0 != sequence:
            return False     
        else:
            return True
  
class SingleFastaPaste(Form):
    def validate_true(form, field):
        if field.data == False:
            raise ValidationError('you must have permission to paste this sequence')
    
    fastaType = RadioField('fastaType',
                           choices = [('protein',"Protein"),
                                      ('rna', "RNA"),
                                      ('dna', "DNA")],
                           validators=[Required()])
    description = TextField('description', 
                            validators=[Required(), 
                                        Length(max=255)])
    sequence = TextAreaField('sequence', 
                         validators=[Required(),
                                     Length(min=1, max = 10000),
                                     ValidateFastaSeq(typefield="fastaType")])
    permission = BooleanField('permission',
                              default=False,
                              validators=[validate_true])
    

class MakeListFromSelf(Form):
    """Accept lists and fastas to make a list """
    username = TextField('username', validators=[Required(),
                                                 Length(min=1, max=128)])    
    userkey = TextField('userkey', validators=[Required(),
                                               Length(min=1, max=128)])
    fastas = TextAreaField("fastas", validators=[Length(max=1000)])
    fastaLists = TextAreaField("fastaLists", validators=[Length(max=200)])
    subtractList = TextField("subtractList", validators=[Length(max=50)])

class UserLoginForm(Form):
    """Accept lists and fastas to make a list """
    username = TextField('username', validators=[Required(),
                                                 Length(min=1, max=128)])    
    userkey = TextField('userkey', validators=[Required(),
                                               Length(min=1, max=128)])
    remember_me = BooleanField("remember_me", default=False)

class DeleteHidden(Form):
    deleteThis = HiddenField("deleteThis", default="yes", validators=[Required(), Regexp("yes")])
    
    
class MakeListFromUniprot(Form):
    def validate_uniprotList(form, field):
        unipRegexStr = r'\A([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}){1}\Z'
        if field.data is None:
            pass
        else:
            uniRe, badResults = re.compile(unipRegexStr), ""
            for id in field.data:
                if uniRe.search(id) is None:
                    badResults += ", " + id
            if len(badResults) > 0:
                raise ValidationError("Non Uniprot Ids found: {}".format(badResults))
                
    uniprotIds = CsvTextAreaField("uniprotIds", validators=[Required(), 
                                                 validate_uniprotList])
    cutSignalSeq = BooleanField("cutSignalSeq", default=False)

