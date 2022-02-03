
import os

def isSameVariant(old_object, new_object):
    old_pseudo_id = old_object['chr'] + old_object['start'] + old_object['ref'] + old_object['alt']
    new_pseudo_id = new_object['chr'] + new_object['start'] + new_object['ref'] + new_object['alt']

    return old_pseudo_id == new_pseudo_id  

def makeDBObject(db_line):
    line = db_line.split('\t')
    db_object = {
        'chr': line[0],	
        'start': line[1],
        'end': line[2],	
        'ref': line[3] if line[3] != '-' else '0',	
        'alt': line[4] if line[4] != '-' else '0',	
        'frq': line[5],
        'nchrobs': line[6].strip()
    }
    return db_object

# TODO los cromosomas no siempre vienen de manera numerica (e.g. 19), pueden venir de forma string (e.g. 'chr19'). Ademas de esto, existen los cromosomas X e Y.
# Acomodar la funcion para considerar estos escenarios
def isBefore(old_object, new_object):
    before = False
    if int(old_object['chr']) > int(new_object['chr']):
        before =  False
    elif int(old_object['chr']) < int(new_object['chr']):
        before =  True
    else:
        before =  int(old_object['start']) < int(new_object['start'])
    
    return before

def mergeDBobjects(old_object, new_object):
    updated_object = new_object.copy()
    obs = int(old_object['nchrobs']) + int(new_object['nchrobs'])  
    frq = (float(old_object['frq'])*float(old_object['nchrobs']) + float(new_object['frq'])*float(new_object['nchrobs']))/ obs
    updated_object['nchrobs'] = str(obs)
    updated_object['frq'] = str(frq) 

    return updated_object

def makeDBline(db_object):
    db_line = db_object['chr'] + '\t' + db_object['start'] + '\t' + db_object['end'] + '\t'
    db_line +=  db_object['ref'] if db_object['ref'] != '0' else '-'
    db_line += '\t'
    db_line += db_object['alt'] if db_object['alt'] != '0' else '-'
    db_line += '\t' + db_object['frq'] + '\t' + db_object['nchrobs'] + '\n'
    return db_line

def updateDB(old_db_path, new_db_path, updated_db_path):
    with open(old_db_path, 'r') as old_db:
        with open(new_db_path, 'r') as new_db:
            with open(updated_db_path, 'x') as updated_db:
                old_line = old_db.readline()
                new_line = new_db.readline()
                while old_line and new_line:
                    old_object = makeDBObject(old_line)
                    new_object = makeDBObject(new_line)

                    if isBefore(old_object, new_object):
                        updated_line = makeDBline(old_object)
                        updated_db.write(updated_line)
                        old_line = old_db.readline()
                    else:
                        if isSameVariant(old_object, new_object):
                            updated_line = makeDBline(mergeDBobjects(old_object, new_object))
                            updated_db.write(updated_line)
                            old_line = old_db.readline()
                            new_line = new_db.readline()
                        else:
                            updated_line = makeDBline(new_object)
                            updated_db.write(updated_line)
                            new_line = new_db.readline()

                if old_line:
                    while old_line:
                        makeDBline(old_line)
                        old_line = old_db.readline()
                elif new_line:
                    while new_line:
                        makeDBline(new_line)
                        new_line = new_db.readline()


class DBUpdater:
    def __init__(self, log = False, tmp = './tmp'):
        self.tmp = tmp
        self.log = log

    def set_old_db_path(self, old_db):
        self.old_db_path = old_db
    
    def set_new_db_path(self, new_db):
        self.old_db_path = new_db

    def set_updated_db_path(self, updated_db):
        self.updated_db_pat = updated_db


OLD_DB_PATH ='old_db.txt'
NEW_DB_PATH ='new_db.txt'
UPDATED_DB_PATH ='updated_db.txt'

if os.path.exists(UPDATED_DB_PATH):
    os.remove(UPDATED_DB_PATH)

updateDB(OLD_DB_PATH, NEW_DB_PATH, UPDATED_DB_PATH)
