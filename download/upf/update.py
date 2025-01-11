'''
update the database indexing the pseudopotentials locally

Details
-------
The ABACUS-Pseudopot-Nao-Square backend high-throughput workflows
manage all pseudopotentials downloaded locally in a tag-search 
manner. Each pseudopotential will be tagged with a list of strings,
including the element symbol, the source of the pseudopotential,
the exchange-correlation functional, etc. The tags are used to
search for the pseudopotentials in the database.

Usage
-----
The run of tagging on pseudopotential requires the rules for tagging.
The rules should be stored in a JSON file, with the following format:
```json
{
    "description": "***",
    "rules": [
        {
            "re.folder": "sg15_oncv_upf_2020-02-06",
            "re.file": "^([A-Z][a-z]?)(_ONCV.*\\.upf)$",
            "tags": ["sg15", 
                     "norm-conserving", "NC", 
                     "ONCV", "ONCVPSP"]
        },
        {
            "re.folder": "sg15_oncv_upf_2020-02-06",
            "re.file": "^([A-Z][a-z]?)(.*_PBE.*\\.upf)$",
            "tags": ["PBE", "Perdew-Burke-Ernzerhof"]
        }
    ]
}
```
, once there is a file satisfies the regex in `re.folder` and 
`re.file`, the tags in `tags` will be added to the file. The tags are
appended to the file, instead of replacing the original tags.
Then you can safely run the following command to update the database
(if the database file is not found, it will be created. But this cannot
be applied to the non-existing case of rules file because it is
considered as a fatal error):
```bash
python update.py
```

Version
-------
20241230
'''

# in-built modules
import os
import json
import unittest
import re
import logging
import time

# APNS special modules
import apns.pspot.parse as ampp
import apns.pspot.parse_special.GBRV_Vanderbilt as amppsg

def _impl(fn_psp_db, fn_rules, pseudo_dir):
    '''
    initialize the database of pseudopotentials and numerical orbitals
    
    Parameters
    ----------
    fn_psp_db: str
        the path to the database of pseudopotentials
        
    fn_rules: str
        the path to the rules for tagging the pseudopotentials
    
    pseudo_dir: str
        the root directory of the pseudopotential files. Pseudopotential files
        are allowed to be placed in any level of subdirectories.
    
    Returns
    -------
    list: a list containing the paths of pseudopotential files that not recognized by rules
    '''
    unlabelled = []
    
    logging.info(f'Initialize the pseudopotential tag-search system database...')
    
    # load the databses and rules
    with open(fn_psp_db) as f:
        psp_db = json.load(f)
    logging.info(f'Pseudopotential database loaded/opened from {fn_psp_db}.')
    
    with open(fn_rules) as f:
        rules = json.load(f)['rules']
    logging.info(f'Rules loaded from {fn_rules}.')
    
    for root, _, files in os.walk(pseudo_dir): # big triangle code is a bad practice...
        logging.info(f'Walking through {root}...')
        for file in files:
            m = re.match(r'.*\.upf$', os.path.basename(file))
            if m:
                # unlike the ABACUS orbitals which places many information in its file name,
                # generally we cannot extract any information from the file name of pseudopotentials
                logging.info(f'Found pseudopotential file {file}.')
                fn = os.path.abspath(os.path.join(root, file)) # key is directly the full path of file
                
                logging.info(f'Transversing rules to add tags to {fn}...')
                for rule in rules:
                    redir, refn, add = rule['re.folder'], rule['re.file'], rule['tags']
                    mm = re.search(redir, root) and re.search(refn, file)
                    if mm:
                        logging.info(f'Valid rule found for {fn}.')
                    else:
                        continue # if no rule is found, continue to next rule

                    # it is the first time to add tags, so add element tag
                    # the parse of pseudopotential will be slow, so need to reduce this operation
                    # as much as possible
                    if not fn in psp_db:
                        pp = ampp.as_dict(fn)
                        ppgencode = ampp.determine_code(pp)
                        if ppgencode == 'GBRV':
                            parsed = amppsg.PP_HEADER(pp['PP_HEADER']['data'])
                            element = parsed['attrib']['element'].strip().lower().capitalize()
                        else:
                            element = pp['PP_HEADER']['attrib']['element'].strip().lower().capitalize()
                        psp_db[fn] = [element]
                        
                    # check if tags are already in the database by comparing sets
                    if set(psp_db.get(fn, [])) != set(add):
                        logging.info(f'Refresh. Adding tags {add} to {fn}...')
                        psp_db[fn] = list(set(psp_db[fn] + add))
                    else:
                        logging.info(f'No need to refresh, have a good day: tags {add} already in {fn}.')
                        
                unlabelled.append(fn) if fn not in psp_db else None
                
    # save the database
    with open(fn_psp_db, 'w') as f:
        json.dump(psp_db, f, indent=4)
    logging.info(f'Pseudopotential database updated and saved to {fn_psp_db}.')
    logging.info(f'Pseudopotential tag-search system database initialized, {len(unlabelled)} files not recognized by rules.')
    
    logging.info(f'Backward check: add `sr` and `scalar-relativistic` to pseudopotentials without `fr` and `full-relativistic`...')
    _postprocess(if_without=['fr', 'full-relativistic'], 
                 add=['sr', 'scalar-relativistic'], 
                 fdb=fn_psp_db)
    logging.info(f'Backward check done.')
    
    return unlabelled

def _postprocess(if_without: list, add: list, fdb: str):
    '''for all pseudopotential files, if without tags in `if_without`, add tags in `add`
    useful for some pseudopotentials marked as 'fr' and 'full-relativistic', while those
     not fr, add 'sr' and 'scalar-relativistic' '''
    assert os.path.exists(fdb), 'database file not found'
    with open(fdb) as f:
        database = json.load(f)
    
    record = []
    for fpp in database.keys():
        if not set(if_without) <= set(database[fpp]):
            print(f'Appending tags {add} to {fpp}...')
            record.append(fpp)
            database[fpp] = list(set(database[fpp] + add))

    with open(fdb, 'w') as f:
        json.dump(database, f, indent=4)

    return record 

def main(prefix, fn_psp_db, fn_rules, pseudo_dir=None):
    '''
    update the database indexing the pseudpotential files locally
    
    Parameters
    ----------
    fn_psp_db: str
        the path to the database of pseudopotentials
        
    fn_rules: str
        the path to the rules for tagging the pseudopotentials
    
    pseudo_dir: str
        the root directory of the pseudopotential files. Pseudopotential files
        are allowed to be placed in any level of subdirectories. If this is
        not set, will use the dirname of fn_psp_db as the root directory.
    
    Returns
    -------
    list: a list containing the paths of pseudopotential files that not recognized by rules
    '''
    # init log
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    fnlog = f'{prefix}@{timestamp}.log'
    print(f'Logging to {fnlog}...')
    logging.basicConfig(filename=fnlog, level=logging.INFO)
    logging.info(f'Update started at {timestamp}...')
    
    # sanity
    if not os.path.exists(fn_psp_db):
        # create the database if not exists
        with open(fn_psp_db, 'w') as f:
            json.dump({}, f)
        logging.info(f'Pseudopotential database created at {fn_psp_db}.')
        
    # kernel
    pseudo_dir = pseudo_dir if pseudo_dir else os.path.dirname(fn_psp_db)
    files = _impl(fn_psp_db, fn_rules, pseudo_dir)
    
    # close log
    logging.info(f'Update completed at {time.strftime("%Y%m%d-%H%M%S")}.')
    logging.shutdown()
    return files

class TestRegularExpression(unittest.TestCase):
    def test_refolder(self):
        '''test re.folder key in rules'''
        PseudoDojov10 = '/root/abacus-develop/ABACUS-Pseudopot-Nao-Square/download/pseudopotentials/pseudos_ac_she/116_Lv/Lv-6spd_r.upf'
        refolder = r'pseudos_ac_she/\d+_[A-Z][a-z]?'
        refile = r'[A-Z][a-z]?-\d+[a-z]*_r\.upf' 
        # will have tags ['Lv', 
        #                 'PseudoDojo', 'DOJO', 'abinit', 
        #                 'v1.0', '1.0', 
        #                 'PBE', 'Perdew-Burke-Ernzerhof',
        #                 'NC', 'norm-conserving',
        #                 'full-relativistic', 'rel',
        #                 ...]
        self.assertTrue(re.search(refolder, PseudoDojov10))
        fupf = PseudoDojov10.split('/')[-1]
        print(fupf)
        self.assertTrue(re.match(refile, PseudoDojov10.split('/')[-1]))

if __name__ == '__main__':

    # for more information about the rules and tags, please read
    # the documentation of the package.
    fdb_psp = '/path/to/your/database/file'
    frules  = '/path/to/your/rules/file'
    
    # run
    unlabelled = main('local_pslib_update', fdb_psp, frules)
    n = len(unlabelled) # number of unlabelled files
    if n > 0:
        print('The unlabelled files are:')
        for fn in unlabelled:
            print(f'  {fn}')
    else:
        print('All pseudopotential files are labelled.')
