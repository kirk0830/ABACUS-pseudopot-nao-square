import os
import re
import json

#########
# utils #
#########
def read_apnsjob_desc(fdesc: str):
    with open(fdesc, 'r') as f:
        desc = json.load(f)

    print('Read APNS job description from file: ', fdesc)
    atom_species_symbols = [as_['symbol'] for as_ in desc['AtomSpecies']]
    s = ', '.join(atom_species_symbols)
    print(f'Atom species: {s}')
    pps = [as_['pp'] for as_ in desc['AtomSpecies']]
    s = '\n'.join(pps)
    print(f'Pseudopotentials bond with:\n{s}')
    # naos = [as_['nao'] for as_ in desc['AtomSpecies']]
    # s = ', '.join(naos)
    # print(f'Number of atomic orbitals: {s}')
    cellgen = desc['CellGenerator']
    print(f'''CellGenerator (where the cell generated from)
identifier: {cellgen['identifier']}
source: {cellgen['config']}
''')
    return desc['AtomSpecies'], desc['CellGenerator']

def _hartwigsen_goedecker_hutter(fpp: str):
    '''
    self-defined parse on hard-coded psp file name of HGH family.
    Adapt for your own case if necessary

    Parameters
    ----------
    fpp : str
        the file name of pseudopotential file
    
    Returns
    -------
    elem, family, version, appendix : str
    '''
    family, version = 'Hartwigsen-Goedecker-Hutter', None
    m = re.match(r'([A-Z][a-z]?\.pbe\-)(.*)(hgh\.UPF)', os.path.basename(fpp))
    elem = m.group(1).replace('.pbe-', '')
    appendix = m.group(2)
    appendix = appendix if appendix else None
    return elem, family, version, appendix

def _pwmat_ncpp_pd04(fpp: str):
    '''
    self-defined parse on hard-coded psp file name of PD04 family
    Adapt for your own case if necessary
    
    Parameters
    ----------
    fpp : str
        the file name of pseudopotential file
    
    Returns
    -------
    elem, family, version, appendix : str
    '''
    family, version = 'PD04', None
    m = re.match(r'([A-Z][a-z]?)([\d\+\-\_\w]*)(\.PD04\.PBE\.UPF)', os.path.basename(fpp))
    elem = m.group(1)
    appendix = m.group(2)
    if appendix:
        if appendix[0] in ['+', '-']:
            appendix = appendix[1:]
    appendix = appendix if appendix else None
    return elem, family, version, appendix

def _garrity_bennett_rabe_vanderbilt(fpp: str):
    '''convert the file name to element, family, version, appendix.
    Change log:
    In the previous version of this function, the version of pseudopotential is hard-coded to be 1.5.
    However, 1.5 is the publication version rather than version of every single pseudopotential file.
    In change at 27th Aug 2024, this has been changed to extract the version from the file name.
    Correspondingly the ultrasoft tests are changed accordingly, here is the script to change all contents
    of data file:
    ```python
    path = '/root/abacus-develop/apns_toupdate'
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))\
            and os.path.basename(f).startswith('ultrasoft-')\
                and os.path.basename(f).endswith('.json')]

    pppath = '/root/abacus-develop/pseudopotentials/GBRV_pbe_UPF_v1.5'
    pps = [f for f in os.listdir(pppath) if os.path.isfile(os.path.join(pppath, f))]

    # a dict mapping element to pseudopotential file
    gbrvpat = r'([a-z]+)(_pbe)((_[\w\d]+)?_v[\d\.]+)(\.uspp\.F\.UPF)'
    gbrv = {}
    for pp in pps:
        m = re.match(gbrvpat, pp)
        if m:
            elem = m.group(1).capitalize()
            version = m.group(3).replace('_v', '')
            version = '1.0' if version == '1' else version
            gbrv[elem] = f'GBRV v{version}'

    for file in files:
        print(f'Processing {file}!')
        # first save a copy
        shutil.copyfile(os.path.join(path, file), os.path.join(path, f'{file}.bak'))
        with open(os.path.join(path, file), 'r') as f:
            data = json.load(f)
        for system in data: # data is a list
            name = system['name']
            elem = name.split()[0].split('-')[0]
            assert re.match(r'[A-Z][a-z]?', elem)
            if 'eos' in file: # eos and ecutwfc convergence test has different structure
                for i, test in enumerate(system['data']): # loop over pseudopotentials, 'data' is a list
                    if test['pp'].startswith('GBRV'):
                        assert elem in gbrv
                        correct = gbrv[elem] # correct with this
                        test['pp'] = correct
            elif 'ecutwfc' in file:
                if system['pp'].startswith('GBRV'):
                    assert elem in gbrv
                    correct = gbrv[elem]
                    system['pp'] = correct
                    
        with open(os.path.join(path, file), 'w') as f:
            json.dump(data, f, indent=4)
    ```
    '''

    m = re.match(r'([a-z]+)(_pbe)((_[\w\d]+)?_v[\d\.]+)(\.uspp\.F\.UPF)', os.path.basename(fpp))
    elem = m.group(1).capitalize()
    version = m.group(3).replace('_v', '')
    version = '1.0' if version == '1' else version
    family, appendix = 'GBRV', None
    return elem, family, version, appendix

def _pwmat_ncpp_pd03(fpp: str):
    '''handle the psp file from PD03 family
    Now the file name is always in the regular pattern of
    ([A-Z][a-z]?)(\.PD03\.PBE\.UPF)

    the group(1) is element, group(2) is appendix

    Parameters
    ----------
    fpp : str
        the file name of pseudopotential file
    
    Returns
    -------
    elem, family, version, appendix : str
    '''
    m = re.match(r'([A-Z][a-z]?)(\.PD03\.PBE\.UPF)', os.path.basename(fpp))
    elem = m.group(1)
    family, version, appendix = 'PD03', None, None
    return elem, family, version, appendix

def _oncvpsp_schlipf_gygi_15(fpp: str):
    '''handle the psp file from SG15 family
    Now the file name is always in the regular pattern of
    ([A-Z][a-z]?(_ONCV_PBE)(_)?(FR)?(\-)(\d\.\d)(\.upf))
    
    the group(1) is element, group(6) is version, group(4) is appendix
    
    Parameters
    ----------
    fpp : str
        the file name of pseudopotential file
    
    Returns
    -------
    elem, family, version, appendix : str
    '''
    family = 'SG15'
    m = re.match(r'([A-Z][a-z]?(_ONCV_PBE)(_)?(FR)?(\-)(\d\.\d)(\.upf))', os.path.basename(fpp))
    elem = m.group(1).split('_')[0].lower().capitalize()
    version = m.group(6)
    appendix = 'fr' if m.group(4) is not None else 'sr'
    return elem, family, version, appendix

def _atomic_pslibrary(fpp: str):
    '''handle the psp file from PSlibrary
    Now the file name is always in the regular pattern of
    ([A-Z][a-z]?)(\.)(rel-)?(pbe|pz)(-\w+)?(-)(rrkjus|kjpaw|nc)(_psl\.)?([\.\d]+)?(\.UPF)
    , a detailed explanation can be found on the Quantum ESPRESSO website:
    https://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format

    the group(1) is element, group(7) is appendix, group(9) is version
    
    Parameters
    ----------
    fpp : str
        the file name of pseudopotential file
    
    Returns
    -------
    elem, family, version, appendix : str
    '''
    family = 'PSlibrary'
    m = re.match(r'([A-Z][a-z]?)(\.)(rel-)?(pbe|pz)(-\w+)?(-)(rrkjus|kjpaw|nc)(_psl\.)?([\.\d]+)?(\.UPF)', os.path.basename(fpp))
    elem = m.group(1).lower().capitalize()
    version = '0.3.1' if m.group(9) is None else m.group(9)
    apps = []
    if m.group(7): apps.append(m.group(7).upper())
    if m.group(3): apps.append('fr')
    if m.group(5): apps.append(m.group(5)[1:])
    appendix = ', '.join(apps)
    return elem, family, version, appendix

def _oncvpsp_pseudo_dojo(fpp: str):
    '''handle the psp file download from the Pseudo-Dojo website
    Now the file name is always in the regular pattern of
    (.*)(nc-sr-05_pbe_standard_upf|nc-fr-04_pbe_standard|pbe_s_sr|
    nc-sr-04_pbe_standard_upf|nc-sr-04-3plus_pbe_standard_upf|pseudos_ac_she)(.*)?
    
    the group(1) is element, group(2) is family, group(3) is version, group(4) is appendix
    
    Parameters
    ----------
    fpp : str
        the file name of pseudopotential file
    
    Returns
    -------
    elem, family, version, appendix : str
    
    Raises
    ------
    ValueError
        if the pseudopotential file is not recognized
    '''
    m = re.match(r'([A-Z][a-z]?)(.*upf)', os.path.basename(fpp))
    elem = m.group(1)
    if 'nc-sr-05_pbe_standard_upf' in fpp:
        family, version, appendix = 'PseudoDojo', '0.5', 'sr'
    elif 'nc-fr-04_pbe_standard' in fpp:
        family, version, appendix = 'PseudoDojo', '0.4', 'fr'
    elif 'pbe_s_sr' in fpp:
        family, version, appendix = 'PseudoDojo', '0.3', 'sr'
    elif 'nc-sr-04_pbe_standard_upf' in fpp or 'nc-sr-04-3plus_pbe_standard_upf' in fpp:
        family, version, appendix = 'PseudoDojo', '0.4', 'sr'
    elif 'pseudos_ac_she' in fpp:
        family, version = 'PseudoDojo', '1.0'
        appendix = 'fr' if fpp.endswith('_r.upf') else 'sr'
    else:
        raise ValueError(f'Unrecognized pseudopotential file: {fpp}')
    
    return elem, family, version, appendix

def _cp2k_goedecker_teter_hutter(fpp: str):
    '''handle the psp file converted from CP2K collection
    Now the file name is always in the regular pattern of
    POTENTIAL_(.*)_([A-Z][a-z]?)(-(.*))?.upf

    the group(1) is version, group(2) is element, group(4) is appendix

    Parameters
    ----------
    fpp : str
        the file name of pseudopotential file
    
    Returns
    -------
    elem, family, version, appendix : str
    '''
    pat = r'POTENTIAL_(.*)_([A-Z][a-z]?)(\-(.*))?.upf'
    m = re.match(pat, os.path.basename(fpp))
    version, elem, appendix = m.group(1), m.group(2), m.group(4)
    return elem, 'Goedecker-Teter-Hutter', version, appendix

def _opium_rappe(fpp):
    '''handle the psp file from Rappe group
    Now the file name is always in the regular pattern of
    (.*)_Rappe(.*)?.UPF

    the group(1) is element, group(2) is appendix

    Parameters
    ----------
    fpp : str
        the file name of pseudopotential file
    
    Returns
    -------
    elem, family, version, appendix : str
    '''
    family, version, appendix = 'Rappe', None, None
    elem = os.path.basename(fpp).split('.')[0].lower().capitalize()
    return elem, family, version, appendix

def convert_fpp_to_ppid(fpp: str):
    '''Convert pseudopotential file name to pseudopotential identifier, the one
    more human readable. The conversion is based on the pseudopotential family
    and version. The appendix is also included in the identifier if it is not empty.'''
    func_map = {
        'hgh': _hartwigsen_goedecker_hutter,
        'NCPP-PD04-PBE': _pwmat_ncpp_pd04,
        'GBRV_pbe_UPF_v1.5': _garrity_bennett_rabe_vanderbilt,
        'NCPP-PD03-PBE': _pwmat_ncpp_pd03,
        'sg15_oncv_upf_2020-02-06': _oncvpsp_schlipf_gygi_15,
        'nc-sr-05_pbe_standard_upf': _oncvpsp_pseudo_dojo,
        'nc-fr-04_pbe_standard': _oncvpsp_pseudo_dojo,
        'pbe_s_sr': _oncvpsp_pseudo_dojo,
        'nc-sr-04_pbe_standard_upf': _oncvpsp_pseudo_dojo,
        'nc-sr-04-3plus_pbe_standard_upf': _oncvpsp_pseudo_dojo,
        'pseudos_ac_she': _oncvpsp_pseudo_dojo,
        'cp2k-collection': _cp2k_goedecker_teter_hutter,
        'psl': _atomic_pslibrary,
        'HGH_NLCC': _cp2k_goedecker_teter_hutter,
        'Rappe': _opium_rappe
    }
    def psp_name(family, version, appendix):
        out = ''
        if family is not None:
            out += f'{family}'
        if version is not None:
            out += f' v{version}'
        if appendix is not None and len(appendix) > 0:
            out += f' ({appendix})'
        return out.replace('v ', '').replace('()', '')
    
    if not isinstance(fpp, str):
        raise ValueError(f'fpp should be a string: {fpp}')
    print(f'Converting {fpp}')
    for key in func_map:
        if key in fpp:
            elem, family, version, appendix = func_map[key](fpp)
            return elem, psp_name(family, version, appendix)
    RuntimeWarning(f'Unrecognized pseudopotential file: {fpp}')
    return 'unknown', fpp

def convert_forb_to_orbid(forb: str):
    '''convert the forbidden string to orbital identifier'''
    forb = os.path.basename(forb)
    rcut_match = re.search(r'\d+(\.\d+)?au', forb)
    ecut_match = re.search(r'\d+(\.\d+)?Ry', forb)
    if not rcut_match or not ecut_match:
        raise ValueError('Invalid format for forbidden string')
    rcut = rcut_match.group(0)
    ecut = ecut_match.group(0)
    conf = forb.split('_')[-1].split('.')[0]
    return f'{rcut}, {ecut} ({conf})'

def cal_dict_diff(desc1: dict, desc2: dict) -> dict:
    '''calculate diff between two dict. For the same key, if the value is different,
    record the difference in tuple, the first element is from desc1, the second is from desc2.
    NOTE: if parent key has different value, the child key will not be compared.
    '''
    if desc1 == desc2:
        return {}
    if not desc1:
        # if desc1 is None, then desc2 should not have tuple value
        if any([isinstance(v, tuple) for v in desc2.values()]):
            raise ValueError('desc2 should not have tuple value')
        return {k: (None, v) for k, v in desc2.items()} if desc2 else {}
    if not desc2:
        # if desc2 is None, then desc1 should not have tuple value
        if any([isinstance(v, tuple) for v in desc1.values()]):
            raise ValueError('desc1 should not have tuple value')
        return {k: (v, None) for k, v in desc1.items()} if desc1 else {}
    
    # if they are not concurrently None, they should be dict
    if not (isinstance(desc1, dict) and isinstance(desc2, dict)):
        raise ValueError('desc1 and desc2 should be dict')
    # make sure before comparison, the value is not tuple
    if not all([not isinstance(v, tuple) for v in desc1.values()]):
        raise ValueError('desc1 should not have tuple value')
    if not all([not isinstance(v, tuple) for v in desc2.values()]):
        raise ValueError('desc2 should not have tuple value')

    diff = {}
    for k, v in desc1.items():
        v_ = desc2.get(k)
        if isinstance(v, dict):
            _diff = cal_dict_diff(v, v_)
            diff.update({k: _diff}) if _diff else None
        else:
            if v != v_:
                diff[k] = (v, v_)
    for k, v in desc2.items():
        if k not in desc1:
            diff[k] = (None, v)
    return diff

def cal_desc_diff(desc1: dict, desc2: dict) -> dict:
    '''more specifically, calculate the difference between two APNS job
    descriptions. Because the AtomSpecies key has value as list of dict,
    therefore make difference for each element in list and return'''
    keys = ['AtomSpecies', 'Cell', 'DFTParamSet', 'CellGenerator']
    assert set(desc1.keys()) == set(desc2.keys()), f'desc1 and desc2 should have the same keys: {desc1.keys()} != {desc2.keys()}'
    assert set(keys) == set(desc1.keys()), f'desc1 should have keys: {keys} != {desc1.keys()}'
    assert set(keys) == set(desc2.keys()), f'desc2 should have keys: {keys} != {desc2.keys()}'
    # the assertation above is to ensure the correctness of structure
    # it must be satisfied before comparison
    diff = {}
    # first the AtomSpecies
    asdiff = [cal_dict_diff(as1, as2) for as1, as2 in zip(desc1['AtomSpecies'], desc2['AtomSpecies'])]
    diff.update({'AtomSpecies': asdiff}) if not all([not asd for asd in asdiff]) else None
    # other keys
    for key in keys[1:]:
        d = cal_dict_diff(desc1[key], desc2[key])
        diff.update({key: d}) if d else None
    return diff

def stru_rev_map(structure: str, basename: bool = False):
    '''export the reverse map from file name to chemical formula'''
    import json, os
    try:
        with open(structure, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f'File not found: {structure}')
        return {}
    except json.JSONDecodeError:
        print(f'Invalid JSON format in file: {structure}')
        return {}
    rev_map_ = {}
    for k, v in data.items():
        for v_ in v:
            f = os.path.basename(v_['file']) if basename else v_['file']
            rev_map_[f] = k + f' ({f})'
    return rev_map_

import unittest
class APNS2UtilsTest(unittest.TestCase):
    def test_cal_dict_diff(self):
        desc1 = {
            'a': 1,
            'b': {
                'c': 2,
                'd': 3
            },
            'e': 4
        }
        desc2 = {
            'a': 1,
            'b': {
                'c': 2,
                'd': 3
            },
            'e': 5
        }
        diff = cal_dict_diff(desc1, desc2)
        self.assertEqual(diff, {'e': (4, 5)})
        diff = cal_dict_diff(desc1, None)
        self.assertEqual(diff, {k: (v, None) for k, v in desc1.items()})
        diff = cal_dict_diff(None, desc2)
        self.assertEqual(diff, {k: (None, v) for k, v in desc2.items()})
        diff = cal_dict_diff(None, None)
        self.assertEqual(diff, {})
        diff = cal_dict_diff({}, {})
        self.assertEqual(diff, {})
        diff = cal_dict_diff({}, desc2)
        self.assertEqual(diff, {k: (None, v) for k, v in desc2.items()})
        diff = cal_dict_diff(desc1, {})
        self.assertEqual(diff, {k: (v, None) for k, v in desc1.items()})
        diff = cal_dict_diff(desc1, desc1)
        self.assertEqual(diff, {})
        diff = cal_dict_diff(desc2, desc2)
        self.assertEqual(diff, {})
        diff = cal_dict_diff(desc2, desc1)
        self.assertEqual(diff, {'e': (5, 4)})

        desc2 = {
            'a': 1,
            'b': {
                'c': 3,
                'd': 3
            },
            'e': 4
        }
        diff = cal_dict_diff(desc1, desc2)
        self.assertEqual(diff, {'b': {'c': (2, 3)}})

        desc3 = {
            'a': 1,
            'b': {
                'c': 2,
                'd': 3,
                'e': 4
            }
        } # someone write the e key in b by mistake
        diff = cal_dict_diff(desc1, desc3)
        self.assertEqual(diff, {'b': {'e': (None, 4)}, 'e': (4, None)})

        desc4 = {
            'a': 1,
            'b': {
                'c': 2
            },
            'e': 4
        } # someone forget to write the d key in b
        diff = cal_dict_diff(desc1, desc4)
        self.assertEqual(diff, {'b': {'d': (3, None)}})

    def test_cal_desc_diff(self):
        desc1 = {
            'AtomSpecies': [
                {'symbol': 'H', 'pp': 'H.pbe-rrkjus.UPF', 'nao': 1},
                {'symbol': 'O', 'pp': 'O.pbe-rrkjus.UPF', 'nao': 2}
            ],
            'Cell': {
                'a': 1,
                'b': 2,
                'c': 3
            },
            'DFTParamSet': {
                'basis_type': 'pw'
            },
            'CellGenerator': {
                'identifier': 'pwscf',
                'config': 'pw.x'
            }
        }
        desc2 = {
            'AtomSpecies': [
                {'symbol': 'H', 'pp': 'H.pbe-rrkjus.UPF', 'nao': 1},
                {'symbol': 'O', 'pp': 'O.pbe-rrkjus.UPF', 'nao': 3}
            ],
            'Cell': {
                'a': 1,
                'b': 2,
                'c': 3
            },
            'DFTParamSet': {
                'basis_type': 'pw'
            },
            'CellGenerator': {
                'identifier': 'pwscf',
                'config': 'pw.x'
            }
        }
        diff = cal_desc_diff(desc1, desc2)
        self.assertEqual(diff, {
            'AtomSpecies': [{}, {'nao': (2, 3)}]
        })

    def test__hartwigsen_goedecker_hutter(self):
        fpp = 'H.pbe-hgh.UPF'
        elem, family, version, appendix = _hartwigsen_goedecker_hutter(fpp)
        self.assertEqual(elem, 'H')
        self.assertEqual(family, 'Hartwigsen-Goedecker-Hutter')
        self.assertEqual(version, None)
        self.assertEqual(appendix, None)
    
    def test_handle_pd04(self):
        fpp = 'H.PD04.PBE.UPF'
        elem, family, version, appendix = _pwmat_ncpp_pd04(fpp)
        self.assertEqual(elem, 'H')
        self.assertEqual(family, 'PD04')
        self.assertEqual(version, None)
        self.assertEqual(appendix, None)
        fpp = 'H-sp.PD04.PBE.UPF'
        elem, family, version, appendix = _pwmat_ncpp_pd04(fpp)
        self.assertEqual(elem, 'H')
        self.assertEqual(family, 'PD04')
        self.assertEqual(version, None)
        self.assertEqual(appendix, 'sp')
        fpp = 'Gd3+_f--core-icmod1.PD04.PBE.UPF'
        elem, family, version, appendix = _pwmat_ncpp_pd04(fpp)
        self.assertEqual(elem, 'Gd')
        self.assertEqual(family, 'PD04')
        self.assertEqual(version, None)
        self.assertEqual(appendix, '3+_f--core-icmod1')

    def test_handle_gbrv(self):
        fpp = 'as_pbe_v1.uspp.F.UPF'
        elem, family, version, appendix = _garrity_bennett_rabe_vanderbilt(fpp)
        self.assertEqual(elem, 'As')
        self.assertEqual(family, 'GBRV')
        self.assertEqual(version, '1.0') # indeed, it is 1.0
        self.assertEqual(appendix, None)

    def test_handle_pd03(self):
        fpp = 'H.PD03.PBE.UPF'
        elem, family, version, appendix = _pwmat_ncpp_pd03(fpp)
        self.assertEqual(elem, 'H')
        self.assertEqual(family, 'PD03')
        self.assertEqual(version, None)
        self.assertEqual(appendix, None)
    
    def test_handle_sg15(self):
        fpp = 'H_ONCV_PBE-1.0.upf'
        elem, family, version, appendix = _oncvpsp_schlipf_gygi_15(fpp)
        self.assertEqual(elem, 'H')
        self.assertEqual(family, 'SG15')
        self.assertEqual(version, '1.0')
        self.assertEqual(appendix, 'sr')
        fpp = 'H_ONCV_PBE_FR-1.0.upf'
        elem, family, version, appendix = _oncvpsp_schlipf_gygi_15(fpp)
        self.assertEqual(elem, 'H')
        self.assertEqual(family, 'SG15')
        self.assertEqual(version, '1.0')
        self.assertEqual(appendix, 'fr')
    
    def test_handle_psl(self):
        fpp = 'H.pbe-rrkjus_psl.1.0.0.UPF'
        elem, family, version, appendix = _atomic_pslibrary(fpp)
        self.assertEqual(elem, 'H')
        self.assertEqual(family, 'PSlibrary')
        self.assertEqual(version, '1.0.0')
        self.assertEqual(appendix, 'RRKJUS')
        fpp = 'H.rel-pbe-kjpaw_psl.0.1.UPF'
        elem, family, version, appendix = _atomic_pslibrary(fpp)
        self.assertEqual(elem, 'H')
        self.assertEqual(family, 'PSlibrary')
        self.assertEqual(version, '0.1')
        self.assertEqual(appendix, 'KJPAW, fr')
        fpp = 'Ac.rel-pbe-spfn-rrkjus_psl.1.0.0.UPF'
        elem, family, version, appendix = _atomic_pslibrary(fpp)
        self.assertEqual(elem, 'Ac')
        self.assertEqual(family, 'PSlibrary')
        self.assertEqual(version, '1.0.0')
        self.assertEqual(appendix, 'RRKJUS, fr, spfn')

    def test_handle_gth(self):
        fpp = 'POTENTIAL_LnPP1_Ce.upf'
        elem, family, version, appendix = _cp2k_goedecker_teter_hutter(fpp)
        self.assertEqual(elem, 'Ce')
        self.assertEqual(family, 'Goedecker-Teter-Hutter')
        self.assertEqual(version, 'LnPP1')
        self.assertEqual(appendix, None)
        
    def test_handle_rappe(self):
        fpp = 'W.Rappe.PBE.UPF'
        elem, family, version, appendix = _opium_rappe(fpp)
        self.assertEqual(elem, 'W')
        self.assertEqual(family, 'Rappe')
        self.assertEqual(version, None)
        self.assertEqual(appendix, None)

if __name__ == '__main__':
    unittest.main()