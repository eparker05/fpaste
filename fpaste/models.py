from fpaste import db

from hashlib import sha256
from string import letters
from random import choice

#this is used for external reference to the declaritive base
declarative_base = db.Model

ROLE_USER = 0
ROLE_ADMIN = 1

class User(db.Model):
    __tablename__ = "user"
    id = db.Column(db.Integer, primary_key = True)
    nickname = db.Column(db.String(64), index = True, unique = True)
    email = db.Column(db.String(120), index = True, unique=True)
    role = db.Column(db.SmallInteger, default=ROLE_USER)
    password = db.Column(db.String(128))
    
    fastas = db.relationship("FastaEntry", backref="user", lazy="dynamic")
    fastaLists = db.relationship("FastaList", backref="user", lazy="dynamic")
    
    def new_password(self, plaintextPw, rounds=10):
        salt = unicode(''.join(choice(letters) for x in range(10)))
        plaintextPw = unicode(plaintextPw)
        digestpw, i = sha256(plaintextPw+salt), 1
        while i < rounds:
            digestpw = sha256(digestpw.digest())
            i += 1
        header = "#{}>>${}>>".format(rounds,salt)
        self.password = header + digestpw.hexdigest()
        
    def validate_password(self, plaintextPw):
        passRounds, salt, realPwHexHash = self.password.split(">>")
        passRounds = int(passRounds[1:])
        salt = unicode(salt[1:])
        digestpw, i = sha256(plaintextPw+salt), 1
        while i < passRounds:
            digestpw = sha256(digestpw.digest())
            i += 1
        if digestpw.hexdigest() == realPwHexHash:
            return True
        else:
            return False
            
    def __repr__(self):
        return('<User {}>'.format(self.nickname))
    
    def is_authenticated(self):
        return True
    def is_active(self):
        return True
    def is_anonymous(self):
        return False
    def get_id(self):
        return unicode(self.id)

list_to_fasta = db.Table("list_to_fasta", db.Model.metadata,
                   db.Column("fasta", db.Integer, db.ForeignKey("fasta_entry.id"), primary_key=True),
                   db.Column("fastaList", db.Integer, db.ForeignKey("fasta_list.id"), primary_key=True)
                        )

class FastaEntry(db.Model):
    __tablename__ = "fasta_entry"
    
    id = db.Column(db.Integer, primary_key = True)
    accessCode = db.Column(db.String(64), index = True, unique=True)
    accession = db.Column(db.String(256), index = True)
    metaData = db.Column(db.String(256), index = True)
    sequence = db.Column(db.Text, index = True)
    added = db.Column(db.DateTime)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))
    
    fastalists = db.relationship("FastaList", secondary=list_to_fasta,
                        backref=db.backref("fastas", lazy=True))
    
    
    
    def read_in_meta_line(self, metaLine):
        #pre-process and strip
        metaLine = metaLine.strip()
        if metaLine[0] == ">":
            metaLine = metaLine[1:].strip()
        
        #take in data, return true for success or false for failure
        if len(metaLine) < 1:
            return False
        else:
            splitIndex = min(254,metaLine.find(' '))
            end = min(256, len(metaLine))
            if splitIndex == "-1":
                self.accession = metaLine
            else:
                self.accession = metaLine[0:splitIndex]
                self.metaData = metaLine[splitIndex:end].strip()
            return True

                   
class FastaList(db.Model):
    __tablename__ = "fasta_list"
    id = db.Column(db.Integer, primary_key = True)
    accessCode = db.Column(db.String(64), index = True)
    user_id = db.Column(db.Integer, db.ForeignKey("user.id"))
    
    # - this is actually populated thanks to the backref command
    #fastaentry = db.relationship("FastaEntry", back_populates="fastalist")
