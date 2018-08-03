from lmfdb.base import getDBConnection
from collections import defaultdict
import time
import re
realre = re.compile(r'([0-9]+)(\.[0-9]*)?')
intre = re.compile(r'[0-9]+')

# detect backslashes
"""
integer
integer stored as string
real
real stored as string
string
tate and or time
boolean
comma separated list of integer stored as string
comma separated list of mixed types stored as string
collection of integer
collection of integer stored as string
collection of real
collection of real stored as string
collection of string
collection of mixed types
collection of comma separated list of mixed types stored as string
Record cannot be found containing key
non-primitive type (<class 'bson.objectid.ObjectId'>)
non-primitive type (<class 'bson.binary.Binary'>) # MaassWaveForms.maassform_plots.plot

collection of comma separated list of integer stored as string
non-primitive type (<type 'NoneType'>) # Lfunctions.full_collection
collection of unidentifiable types # knowledge.knowls.history

"""
skip_collections = ['Lfunctions.full_collection',
                    'transitivegroups.groups-old',
                    'knowledge.deleted_knowls',
                    'knowledge.meta',
                    'Lattices.lat_new',
                    'Lattices.lat_nf',  # is_lat_in = [[1334, -10279], [-10279, 84566]], should be string
                    'Lfunction.FarmerTest',
                    'Lfunction.LemurellMaassHighDegree',
                    'Lfunction.LemurellTest',
                    'Lfunctions.new_Lfunctions',
                    'Lfunctions.new_Linstances',
                    'Lfunctions.ecqd1_Lfunctions',
                    'Lfunctions.ecqd1_instances',
                    'MaassWaveForms.HT',
                    'MaassWaveForms.maassform',
                    'MaassWaveForms.metadata',
                    'MaassWaveForms.DirChars',
                    'MaassWaveForms.DirCharsConrey',
                    'MaassWaveForms.DirCharsSage',
                    'MaassWaveForms.CoeffProgress',
                    'SL2Zsubgroups.groups2',
                    'characters.test_collection',
                    'characters.Dirichlet_char_modl',
                    'characters.Dirichlet_gp_modl',
                    'elliptic_curves.test_collection',
                    'embedded_mfs.hecke_orbits',
                    'embedded_mfs.mfs',
                    'hecke_algebras.stats',
                    'hecke_algebras.test_collection',
                    'hgm.motives_copy',
                    'itest.itest',
                    'localfields.back4',
                    'mod_l_eigenvalues.Dirichlet_char_modl',
                    'mod_l_eigenvalues.Dirichlet_gp_modl',
                    'mod_l_eigenvalues.finite_fields',
                    'mod_l_eigenvalues.hecke_algebras',
                    'mod_l_eigenvalues.hecke_algebras_l_adic',
                    'mod_l_eigenvalues.hecke_algebras_orbits',
                    'siegel_modular_forms.paramodular_lifts',

                    'modularforms2.computations',
                    'modularforms2.dimension_table',
                    'modularforms2.webchar',
                    'modularforms2.webmodformspace',
                    'modularforms2.webnewforms',

                    # 'hgm.families',
                    # 'hgm.motives',
                    'hgm.newfamilies',
                    'hgm.newmotives',

                    'artin.representations',
                    'artin.representations_new',
                    'artin.field_data',
                    'artin.field_data_new',
                    'knowledge.knowls',
                    'knowledge.history',
                    'elliptic_curves.curves',
                    'numberfields.fields',
                    'abvar.fq_isog',
                    'MaassWaveForms.Coefficients',
                    'MaassWaveForms.FS',
                    'MaassWaveForms.maassform_plots',
                    'siegel_modular_forms.samples',
                    #'SL2Zsubgroups.groups',

                    #'elliptic_curves.padic_db',
                    #'sato_tate_groups.small_groups',
                    #'sato_tate_groups.st_groups',
                    #'sato_tate_groups.st0_groups',
                    #'transitivegroups.groups',
                    #'transitivegroups.Gmodules',

                    # 'HTPicard.picard',
                    # 'Lattices.lat',
                    #'Lfunctions.Lfunctions',
                    #'Lfunctions.instances',
                    # 'bmfs.dimensions',
                    # 'bmfs.forms',
                    # 'curve_automorphisms.passports',
                    # 'elliptic_curves.IQF_labels',
                    # 'elliptic_curves.nfcurves',
                    # 'finite_fields.finite_fields',
                    # 'genus2_curves.curves',
                    # 'genus2_curves.endomorphisms',
                    # 'genus2_curves.ratpts',
                    # 'genus2_curves.tamagawa_numbers',
                    # 'halfintegralmf.forms',
                    # 'hecke_algebras.hecke_algebras',
                    # 'hecke_algebras.hecke_algebras_l_adic',
                    # 'hecke_algebras.hecke_algebras_orbits',
                    # 'hmfs.fields',
                    # 'hmfs.forms',
                    # 'hmfs.hecke',
                    # 'localfields.fields',
                    # 'mod_l_eigenvalues.modlmf',
                    # 'mod_l_galois.reps',
                    # 'siegel_modular_forms.dimensions',
                    # 'siegel_modular_forms.families',

]

id_ordered = ['artin_reps', 'av_fqisog', 'bmf_forms', 'ec_curves', 'ec_nfcurves',
              'hgcwa_passports', 'lfunc_lfunctions', 'mwf_forms', 'nf_fields', 'smf_fc']

sorts = {'HTPicard.picard': ['ev'],
         'Lattices.lat': ['dim', 'det', 'level', 'class_number', 'label', 'minimum', 'aut'],
         'Lfunctions.Lfunctions': ['degree','group','conductor','Lhash'],
         'Lfunctions.instances': None,
         'MaassWaveForms.FS': ['Weight','Level','Character','Eigenvalue'],
         'MaassWaveForms.maassform_plots': None,
         'MaassWaveForms.Coefficients': None,
         'SL2Zsubgroups.groups': None,
         'abvar.fq_isog': ['g', 'q', 'poly'],
         'artin.representations': ['Dim','Conductor'],
         'artin.field_data': None,
         'belyi.galmaps': ['deg', 'group_num', 'g', 'label'],
         'belyi.passports': ['deg', 'group', 'g', 'plabel'],
         'bmfs.forms': ['level_norm','label'],
         'bmfs.dimensions': ['level_norm','label'],
         'curve_automorphisms.passports': ['genus','dim','cc'],
         'elliptic_curves.curves': ['conductor', 'iso_nlabel', 'lmfdb_number'],
         'elliptic_curves.IQF_labels': None,
         'elliptic_curves.nfcurves': ['field_label','conductor_norm','conductor_label','iso_nlabel','number'],
         'elliptic_curves.padic_db': None,
         'finite_fields.finite_fields': None,
         'genus2_curves.curves': ['cond','class','abs_disc', 'disc_sign','label'],
         'genus2_curves.endomorphisms': None,
         'genus2_curves.ratpts': None,
         'genus2_curves.tamagawa_numbers': ['label','p'],
         'halfintegralmf.forms': ['level','label'],
         'hecke_algebras.hecke_algebras': ['level','weight','num_orbits'],
         'hecke_algebras.hecke_algebras_l_adic': ['level','weight','orbit_label'],
         'hecke_algebras.hecke_algebras_orbits': ['parent_label','orbit'],
         #'hgm.newmotives': ['cond','label'],
         #'hgm.newfamilies': ['label'],
         'hgm.motives': ['cond','label'],
         'hgm.families': ['label'],
         'hmfs.fields': None,
         'hmfs.forms': ['deg','disc','level_norm','level_label','label_nsuffix'],
         'hmfs.hecke': None,
         'localfields.fields': ['p','n','c','label'],
         'mod_l_eigenvalues.modlmf': ['characteristic','deg','level','weight_grading'],
         'mod_l_galois.reps': ['dim','field_order', 'conductor', 'label'],
         'modularforms2.computations': None,
         'modularforms2.dimension_table': None,
         'modularforms2.webchar': None,
         'modularforms2.webmodformspace': None,
         'modularforms2.webnewforms': None,
         'numberfields.fields': ['degree', 'disc_abs', 'disc_sign', 'iso_number'],
         'sato_tate_groups.small_groups': ['order','label'],
         'sato_tate_groups.st_groups': ["weight", "degree", "real_dimension", "identity_component", "name"],
         'sato_tate_groups.st0_groups': ["name"],
         'siegel_modular_forms.dimensions': None,
         'siegel_modular_forms.families': ['ord'],
         'siegel_modular_forms.samples': ['name'],
         'transitivegroups.groups': ["n", "t"],
         'transitivegroups.Gmodules': ["n", "t", "index"],
}

fixes = {'Lfunctions.Lfunctions':{'gamma_factors':'collection of integer', ## should be json, but some as strings
                                  'root_number':'string',
                                  'leading_term':'string',
                                  'central_character':'string',
                                  'sign_arg':'real'},
         'Lfunctions.ecqd1_Lfunctions':{'gamma_factors':'collection of integer'}, ## should be json
         'abvar.fq_isog':{'gal':False,
                          'sort':False,
                          'C_cnts':'space separated list of integer stored as string',
                          'A_cnts':'space separated list of integer stored as string',
                          'poly':'space separated list of integer stored as string'},
         'artin.representations':{'Conductor_key':False},
         'belyi.galmaps': {
             'embeddings': 'collection of real',
             'label': 'string',
             'plabel': 'string',
             },
         'belyi.galmaps': {
             'plabel' : 'string',
             },
         'bmfs.dimensions':{'level_label':'string'},
         'bmfs.forms':{'Lratio':'string',
                       'level_label':'string'},
         'genus2_curves.curves':{'disc_key':False,
                                 'Lhash':'string'},
         'elliptic_curves.nfcurves':{'reg':'real',
                                     'conductor_label':'string'},
         'hecke_algebras.hecke_algebras_l_adic':{'num_gen_l':False,
                                                 'gen_l':False,
                                                 'rel_l':False,
                                                 'charpoly_ql':False},
         'hecke_algebras.hecke_algebras_orbits':{'gen':False,
                                                 'num_gen':False,
                                                 'rel':False},
         'hgm.motives':{"A":'comma separated list of integer stored as string',
                        "A2":'comma separated list of integer stored as string',
                        "A3":'comma separated list of integer stored as string',
                        "A5":'comma separated list of integer stored as string',
                        "A7":'comma separated list of integer stored as string',
                        "A3rev":'comma separated list of integer stored as string',
                        "A5rev":'comma separated list of integer stored as string',
                        "A7rev":'comma separated list of integer stored as string',
                        "Au2":'comma separated list of integer stored as string',
                        "Au3":'comma separated list of integer stored as string',
                        "Au5":'comma separated list of integer stored as string',
                        "Au7":'comma separated list of integer stored as string',
                        "Au2rev":'comma separated list of integer stored as string',
                        "B3":'comma separated list of integer stored as string',
                        "B5":'comma separated list of integer stored as string',
                        "B7":'comma separated list of integer stored as string',
                        "Brev":'comma separated list of integer stored as string',
                        "B2rev":'comma separated list of integer stored as string',
                        "B3rev":'comma separated list of integer stored as string',
                        "B5rev":'comma separated list of integer stored as string',
                        "B7rev":'comma separated list of integer stored as string',
                        "Bu2":'comma separated list of integer stored as string',
                        "Bu2rev":'comma separated list of integer stored as string',
                        "Bu3rev":'comma separated list of integer stored as string',
                        "Bu5rev":'comma separated list of integer stored as string',
                        "Bu7rev":'comma separated list of integer stored as string',
                        "C2":'comma separated list of integer stored as string',
                        "Cu3":'comma separated list of integer stored as string',
                        "Cu5":'comma separated list of integer stored as string',
                        "Cu7":'comma separated list of integer stored as string',
                        "t":'string',
                        "locinfo":'ijson',
         },
         'hgm.families':{"leader":'string',
                         "A":'comma separated list of integer stored as string',
                         "Arev":'comma separated list of integer stored as string',
                         "A2":'comma separated list of integer stored as string',
                         "A3":'comma separated list of integer stored as string',
                         "A5":'comma separated list of integer stored as string',
                         "A7":'comma separated list of integer stored as string',
                         "A2rev":'comma separated list of integer stored as string',
                         "A7rev":'comma separated list of integer stored as string',
                         "Au2":'comma separated list of integer stored as string',
                         "Au3":'comma separated list of integer stored as string',
                         "Au5":'comma separated list of integer stored as string',
                         "Au7":'comma separated list of integer stored as string',
                         "Au3rev":'comma separated list of integer stored as string',
                         "Au5rev":'comma separated list of integer stored as string',
                         "Au7rev":'comma separated list of integer stored as string',
                         "B":'comma separated list of integer stored as string',
                         "Brev":'comma separated list of integer stored as string',
                         "B2":'comma separated list of integer stored as string',
                         "B2rev":'comma separated list of integer stored as string',
                         "B3rev":'comma separated list of integer stored as string',
                         "B5rev":'comma separated list of integer stored as string',
                         "B7rev":'comma separated list of integer stored as string',
                         "Bu3":'comma separated list of integer stored as string',
                         "Bu5":'comma separated list of integer stored as string',
                         "Bu7":'comma separated list of integer stored as string',
                         "Bu2rev":'comma separated list of integer stored as string',
                         "Bu3rev":'comma separated list of integer stored as string',
                         "Bu5rev":'comma separated list of integer stored as string',
                         "Bu7rev":'comma separated list of integer stored as string',
                         "C2":'comma separated list of integer stored as string',
                         "C3":'comma separated list of integer stored as string',
                         "C5":'comma separated list of integer stored as string',
                         "C7":'comma separated list of integer stored as string',
                         "Cu2":'comma separated list of integer stored as string',
                         "Cu3":'comma separated list of integer stored as string',
                         "Cu5":'comma separated list of integer stored as string',
                         "Cu7":'comma separated list of integer stored as string',
                         "snf":'comma separated list of integer stored as string',
                         "famhodge":'comma separated list of integer stored as string',
         },
         'hmfs.forms':{'level_label':'string'},
         'localfields.fields':{'gal':'galois pair',
                               'gms':'string',
                               'top_slope':'string'}, # '03.83333333323/6' 1040 times??
         'mod_l_galois.reps':{'related_objects':False},
         'modularforms2.webmodformspace':{'zetas':False},
         'modularforms2.webnewforms':{'sage_version':False},
         'modularforms2.webchar':{'label':'string'},
         'modularforms2.dimension_table':{'space_labe;':False,
                                          'gamma1_label':'string'},
         'sato_tate_groups.st_groups':{'trace_zero_density':'string',
                                       'component_group':'string'},
         'sato_tate_groups.small_groups':{'derived_group':'string',
                                          'center':'string',
                                          'abelian_quotient':'string',
                                          'label':'string'},
         'siegel_modular_forms.dimensions':{'311':False},
         'transitivegroups.groups':{'orderkey':False},
}

drop_indexes = {'Lfunctions.Lfunctions': ['a1_2d', 'a2_2d', 'a3_2d', 'a4_2d', 'a5_2d', 'a6_2d', 'a7_2d', 'a8_2d', 'a9_2d', 'gamma2_1', 'instances_1', 'gamma1_1', 'gamma3_1'],
                'MaassWaveForms.FS': ["[ ('Weight', 1), ('Level', 1), ('Character', 1), ('Eigenvalue', 1) ]_1", "_cls_1"],
                'abvar.fq_isog': ["principally_polarizable_1_g_1_q_1_label_1","slopes_1_g_1_q_1_label_1","slopes_1","A_counts_1","known_jacobian_1","C_counts_1","p_rank_1_polynomial_1","known_jacobian_1_g_1_q_1_label_1","C_counts_1_g_1_q_1_label_1","decomposition_1","polynomial_1_g_1_q_1_label_1","decomposition_1_g_1_q_1_label_1","polynomial_1","p_rank_1_polynomial_1_g_1_q_1_label_1","angle_ranks_1_g_1_q_1_label_1","principally_polarizable_1","number_field_1_g_1_q_1_label_1", "gal_1_sort_1", "sort_1", "p_rank_1"],
                'elliptic_curves.nfcurves': ["torsion_orders_1"],
                'mod_l_galois.reps': ["type_1"],
                'modularforms2.webmodformspace': ["level_1_weight_1_chi_1"],
                'modularforms2.webnewforms': ["level_1_weight_1_chi_1"],
                'sato_tate_groups.st0_groups': ["degree_1", "weight_1", "name_1", "label_1"],
		'transitivegroups.groups': ["orderkey_1_n_1_t_1"],
}

mod_indexes = {'abvar.fq_isog': {"gal_1": ['galois_n', 'galois_t', 'id'],
                                 "q_1_gal_1_sort_1": ['q', 'galois_n', 'galois_t', 'id']},
               'artin.representations': {'Hide_1_Conductor_key_1': ['Hide', 'Conductor'],
                                         'Dim_1_Hide_1_Conductor_key_1': ['Dim', 'Hide', 'Conductor']},
               'genus2_curves.curves': {'disc_key_1': ['abs_disc'],
                                        'cond_1_class_1_disc_key_1_label_1': ['cond', 'class', 'abs_disc', 'label']},
               'hgm.families': {'hodge_1': ['famhodge']},
               'hgm.motives': {'a2_1': ['A2'],
                               'a3_1': ['A3'],
                               'a5_1': ['A5'],
                               'a7_1': ['A7'],
                               'b2_1': ['B2'],
                               'b3_1': ['B3'],
                               'b5_1': ['B5'],
                               'b7_1': ['B7'],},
}
add_indexes = {'elliptic_curves.nfcurves':{'conductor_norm': ['conductor_norm'],
                                            'base_change_gin': ['base_change']}}

renames = {'HTPicard.picard': 'mwfp_forms',
           'Lattices.lat': 'lat_lattices',
           #'Lattices.lat_nf': 'lat_nflattices',
           'Lfunctions.Lfunctions': 'lfunc_lfunctions',
           'Lfunctions.instances': 'lfunc_instances',
           'MaassWaveForms.FS': 'mwf_forms',
           'MaassWaveForms.maassform_plots': 'mwf_plots',
           'MaassWaveForms.Coefficients': 'mwf_coeffs',
           'SL2Zsubgroups.groups': 'sl2z_subgroups',
           'abvar.fq_isog': 'av_fqisog',
           'artin.representations': 'artin_reps',
           'artin.field_data': 'artin_field_data',
           'belyi.galmaps': 'belyi_galmaps',
           'belyi.passports': 'belyi_passports',
           'bmfs.forms': 'bmf_forms',
           'bmfs.dimensions': 'bmf_dims',
           'curve_automorphisms.passports': 'hgcwa_passports',
           'elliptic_curves.IQF_labels': 'ec_iqf_labels',
           'elliptic_curves.nfcurves': 'ec_nfcurves',
           'finite_fields.finite_fields': 'fq_fields',
           'genus2_curves.curves': 'g2c_curves',
           'genus2_curves.endomorphisms': 'g2c_endomorphisms',
           'genus2_curves.ratpts': 'g2c_ratpts',
           'genus2_curves.tamagawa_numbers': 'g2c_tamagawa',
           'halfintegralmf.forms': 'halfmf_forms',
           'hecke_algebras.hecke_algebras': 'hecke_algebras',
           'hecke_algebras.hecke_algebras_l_adic': 'hecke_ladic',
           'hecke_algebras.hecke_algebras_orbits': 'hecke_orbits',
           'hgm.motives': 'hgm_motives',
           'hgm.families': 'hgm_families',
           'hmfs.fields': 'hmf_fields',
           'hmfs.forms': 'hmf_forms',
           'hmfs.hecke': 'hmf_hecke',
           'localfields.fields': 'lf_fields',
           'mod_l_eigenvalues.modlmf': 'modlmf_forms',
           'mod_l_galois.reps': 'modlgal_reps',
           #'modularforms2.computations': 'mf_progress',
           #'modularforms2.dimension_table': 'mf_dims',
           #'modularforms2.webchar': 'mf_chars',
           #'modularforms2.webmodformspace': 'mf_spaces',
           #'modularforms2.webnewforms': 'mf_newforms',
           'siegel_modular_forms.dimensions': 'smf_dims',
           'siegel_modular_forms.families': 'smf_families',
           'siegel_modular_forms.samples': 'smf_samples',
           'knowledge.knowls': 'kwl_knowls',
           'knowledge.history': 'kwl_history',
           'elliptic_curves.curves': 'ec_curves',
           'elliptic_curves.padic_db': 'ec_padic',
           'sato_tate_groups.small_groups': 'gps_small',
           'sato_tate_groups.st_groups': 'gps_sato_tate',
           'sato_tate_groups.st0_groups': 'gps_sato_tate0',
           'transitivegroups.groups': 'gps_transitive',
           'transitivegroups.Gmodules': 'gps_gmodules',
           'numberfields.fields': 'nf_fields',
}

label_cols = {'HTPicard.picard': 'maass_id',
              'Lattices.lat': 'label',
              'Lfunctions.Lfunctions': 'Lhash',
              'Lfunctions.instances': 'url',
              'SL2Zsubgroups.groups': None,
              'abvar.fq_isog': 'label',
              'belyi.galmaps': 'label',
              'belyi.passports': 'plabel',
              'bmfs.forms': 'label',
              'bmfs.dimensions': 'label',
              'curve_automorphisms.passports': 'total_label',
              'elliptic_curves.IQF_labels': None,
              'elliptic_curves.nfcurves': 'label',
              'finite_fields.finite_fields': None,
              'genus2_curves.curves': 'label',
              'genus2_curves.endomorphisms': 'label',
              'genus2_curves.ratpts': 'label',
              'genus2_curves.tamagawa_numbers': 'label',
              'halfintegralmf.forms': 'label',
              'hecke_algebras.hecke_algebras': 'label',
              'hecke_algebras.hecke_algebras_l_adic': 'label_l',
              'hecke_algebras.hecke_algebras_orbits': 'orbit_label',
              'hgm.motives': 'label',
              'hgm.families': 'label',
              'hmfs.fields': 'label',
              'hmfs.forms': 'label',
              'hmfs.hecke': 'label',
              'localfields.fields': 'label',
              'mod_l_eigenvalues.modlmf': 'label',
              'mod_l_galois.reps': 'label',
              'siegel_modular_forms.dimensions': None,
              'siegel_modular_forms.families': None,
              'siegel_modular_forms.samples': None,
              'elliptic_curves.padic_db': None,
              'sato_tate_groups.small_groups': 'label',
              'sato_tate_groups.st_groups': 'label',
              'sato_tate_groups.st0_groups': 'label',
              'transitivegroups.groups': 'label',
              'transitivegroups.Gmodules': None,

              'elliptic_curves.curves': 'lmfdb_label',
              'numberfields.fields': 'label',
              'MaassWaveForms.FS': 'maass_id',
              'MaassWaveForms.maassform_plots': 'maass_id',
              'MaassWaveForms.Coefficients': 'label',
              'mwf_tables': None, # new table
              'artin.representations': 'Baselabel',
              'artin.field_data': None,
}

col_renames = {}

def _build_typecol(coltypes, included_types, normal, wrapper, postgres_type, collection=False):
    if collection:
        return ', '.join('"{}"'.format(col) for col in coltypes if coltypes[col].startswith('collection of')), normal, wrapper, postgres_type
    else:
        return ', '.join('"{}"'.format(col) for col in coltypes if coltypes[col] in included_types), normal, wrapper, postgres_type

def find_records():
    inv = getDBConnection().inventory
    databases = {db['_id']:db['name'] for db in inv.DB_ids.find()}
    collection_names = {coll['_id']:"%s.%s"%(databases[coll['db_id']], coll['name']) for coll in inv.collection_ids.find({'status':0})}
    records = defaultdict(list)
    for rec in inv.records.find():
        collid = rec['coll_id']
        try:
            cname = collection_names[collid]
        except KeyError:
            continue
        if cname in skip_collections:
            continue
        records[cname].append(rec['schema'])
    for cname in sorted(records.keys()):
        if len(records[cname]) > 1:
            print cname, [len(schema) for schema in records[cname]]
        else:
            del records[cname]
    return records

def rough(exportfile, importfile):
    inv = getDBConnection().inventory
    databases = {db['_id']:db['name'] for db in inv.DB_ids.find()}
    #collections = {coll['_id']:getDBConnection()[databases[coll['db_id']]][coll['name']] for coll in inv.collection_ids.find({'status':status})}
    collection_names = {coll['_id']:"%s.%s"%(databases[coll['db_id']], coll['name']) for coll in inv.collection_ids.find() if coll['status'] in [0,2]}
    # remove duplicates that just differ on the presence of "sort"
    sort_dup = defaultdict(set)
    for index in inv.indexes.find():
        name = index['name']
        cname = collection_names.get(index['coll_id'])
        if cname is not None and name.endswith('_sort_1'):
            sort_dup[cname].add(name[:-7])
    indexes = defaultdict(list)
    for index in inv.indexes.find():
        name = index['name']
        cname = collection_names.get(index['coll_id'])
        if cname is not None and name not in drop_indexes.get(cname,[]) and name not in sort_dup.get(cname,[]):
            mod = mod_indexes.get(cname,{}).get(name)
            if mod is not None:
                indexes[index['coll_id']].append((index['name'], [[key,1] for key in mod]))
            else:
                indexes[index['coll_id']].append((index['name'], index['keys']))
    for cname in add_indexes:
        for collid, _cname in collection_names.iteritems():
            if cname == _cname:
                break
        else:
            raise RuntimeError("Collection not found")
        for iname, sort_keys in add_indexes[cname].iteritems():
            indexes[collid].append((iname, [[key,1] for key in sort_keys]))
    for ilist in indexes.itervalues():
        ilist.sort()
    examples = defaultdict(list)
    types = defaultdict(dict)
    errors = 0
    for fld in inv.fields_auto.find():
        collid = fld['coll_id']
        try:
            cname = collection_names[collid]
        except KeyError:
            # wrong status
            continue
        if cname in skip_collections:
            continue
        typ = fld['data']['type']
        colname = fld['name']
        fix = fixes.get(cname, {}).get(colname)
        if fix is False:
            continue
        elif fix is not None:
            typ = fix
        newname = col_renames.get(cname, {}).get(colname)
        if newname is not None:
            colname = newname
        examples[typ].append(fld['data']['example'])
        types[cname][colname] = typ
    unames = []
    with open(exportfile, 'w') as Fexport:
        if True: # don't want to fix indentation level
            Fexport.write(r"""
from collections import defaultdict
from datetime import datetime, date, time
import json, re, traceback
import lmfdb.base
from sage.functions.other import ceil
from sage.functions.trig import arctan2
from sage.rings.real_mpfr import RealLiteral, RealField, RealNumber
from sage.rings.complex_number import ComplexNumber
from sage.rings.complex_field import ComplexField
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import infinity
from sage.rings.number_field.number_field_element import NumberFieldElement
from sage.rings.number_field.number_field import NumberField, CyclotomicField, NumberField_generic, NumberField_absolute, NumberField_cyclotomic
from sage.rings.number_field.number_field_rel import NumberField_relative
from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_dense_field
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.power_series_poly import PowerSeries_poly
from sage.modules.free_module_element import vector, FreeModuleElement
stripzero_re = re.compile(r'^([0-9-]+)(\.0*)?$')
def intwrap(n):
    if n is None:
        return r'\N'
    if not isinstance(n, basestring):
        n = str(n)
    if not n:
        return r'\N'
    imatch = stripzero_re.match(n)
    if imatch:
        return imatch.group(1)
    else:
        raise ValueError(n)
RR = RealField(53)
rational_re = re.compile(r'^\d+/\d+$')
decimal_re = re.compile(r'^-?\d*\.?\d*(e[+-]?\d+)?$')
def realwrap(n):
    if n is None:
        return r'\N'
    if not isinstance(n, basestring):
        n = str(n)
    if not n:
        return r'\N'
    if rational_re.match(n):
        return str(RR(QQ(n))) # is it okay to use 53-bit precision here?
    if decimal_re.match(n):
        return n
    else:
        raise ValueError(n)
def datetimewrap(n):
    return r'\N' if n is None else str(n)
def binarywrap(n):
    return r'\N' if n is None else r'\\x' + ''.join(c.encode('hex') for c in n)
def asiswrap(n):
    if n is None:
        return r'\N'
    elif isinstance(n, basestring):
        return n
    else: # sometimes inventory lies to us
        return str(n)
def strwrap(n):
    if n is None:
        return r'\N'
    if not isinstance(n, basestring):
        n = str(n)
    return n.replace('\\', '\\\\').replace("\r", r"\r").replace("\n", r"\n").replace("\t", r"\t").replace('"',r'\"')
def boolwrap(b):
    if b is None:
        return r'\N'
    elif b:
        return 't'
    else:
        return 'f'
all_digits = re.compile(r'^-?[0-9]+$')
def _int_or_fix(x):
    if all_digits.match(x):
        return int(x)
    return x.replace('\\','\\\\\\\\')
def splitwrap(n):
    if n is None:
        return r'\N'
    elif not n: # empty string
        n = []
    elif ',' in n:
        n = [_int_or_fix(x.strip()) for x in n.split(',')]
    else:
        n = [_int_or_fix(x.strip()) for x in n.split(' ')]
    return json.dumps(n)
def jsonwrap(n):
    return r'\N' if n is None else Json.dumps(n)
def ijsonwrap(n):
    return r'\N' if n is None else Json.dumps(n, True)
def galpairwrap(n):
    return r'\N' if n is None else str(n[1])
one_re = re.compile(r'1(\.0*)?')
minusone_re = re.compile(r'-1(\.0*)?')
complex_re = re.compile(r'(-?[0-9]?(\.[0-9e\-]*)?)\s*([+-]\s*[0-9]?(\.[0-9e\-]*))\*[Ii]')
def anglewrap(n):
    if n is None:
        return r'\N'
    if isinstance(n, int):
        if n == 1:
            return '0'
        elif n == -1:
            return '0.5'
    if isinstance(n, basestring):
        if one_re.match(n):
            return '0'
        elif minusone_re.match(n):
            return '0.5'
        else:
            cmatch = complex_re.match(n)
            if cmatch:
                gps = cmatch.groups()
                real, imag = gps[0], gps[2]
                prec = ceil(max(len(gps[1]), len(gps[3])) * 3.322) # 3.322 = log_2(10)
                R = RealField(prec)
                pi = R.pi()
                angle = arctan2(R(imag),R(real))/(2*pi)
                if angle < 0:
                    angle += 1
                return str(angle)
    raise ValueError("Unable to parse %s (%s)"%(n, type(n)))
def integer_size(n):
    if n < 2**7:
        return 'smallint'
    elif n < 2**15:
        return 'integer'
    elif n < 2**31:
        return 'bigint'
    else:
        return 'numeric'
def maxval_update(flds, maxvals, fld, rec):
    val = rec.get(fld)
    if val in ['','?']: # another way of indicating null for some tables
        flds[fld] = r'\N'
        return
    if isinstance(val, basestring) and '.' in val: # sometimes the inventory lies a lot
        val = int(round(float(val)))
    flds[fld] = intwrap(val)
    if val is not None:
        val = abs(int(val))
        if val > maxvals[fld]:
            maxvals[fld] = val
def report_time(i, total, name, t0):
    if i and i % 10000 == 0:
        t = datetime.now()
        print "%s: %s/%s ~ %.2f  in %s"%(name, i, total, i/(1.0 * total),  t-t0)
def sort_collection(coll, sort, name):
    t0 = datetime.now()
    total = coll.count()
    if sort is None:
        for i, rec in enumerate(coll.find()):
            report_time(i, total, name, t0)
            yield i, rec
    else:
        keys = []
        sort.append("_id")
        cur = coll.find({}, sort, no_cursor_timeout=True)
        try:
            for i, rec in enumerate(cur):
                report_time(i, total, name + " sort", t0)
                keys.append(tuple(rec.get(k) for k in sort))
        finally:
            cur.close()
        keys.sort()
        id_dict = {}
        for i, key in enumerate(keys):
            id_dict[key[-1]] = i + 1
        cur = coll.find({}, no_cursor_timeout=True)
        try:
            for i, rec in enumerate(cur):
                report_time(i, total, name, t0)
                yield id_dict[rec["_id"]], rec
        finally:
            cur.close()
class Json(object):
    @classmethod
    def dumps(cls, obj, make_int=False):
        return json.dumps(cls.prep(obj, make_int=make_int))

    @classmethod
    def prep(cls, obj, make_int=False, escape_backslashes=True):
        # For now we just hard code the encoding.
        # It would be nice to have something more abstracted/systematic eventually
        if isinstance(obj, tuple):
            return cls.prep(list(obj), make_int, escape_backslashes)
        elif isinstance(obj, list):
            # Lists of complex numbers occur, and we'd like to save space
            # We currently only support Python's complex numbers
            # but support for Sage complex numbers would be easy to add
            if obj and all(isinstance(z, complex) for z in obj):
                return {'__ComplexList__': 0, # encoding version
                        'data': [[z.real, z.imag] for z in obj]}
            elif obj and all(isinstance(z, Rational) for z in obj):
                return {'__QQList__': 0, # encoding version
                        'data': [[int(z.numerator()), int(z.denominator())] for z in obj]}
            elif obj and all(isinstance(z, NumberFieldElement) for z in obj) and all(z.parent() is obj[0].parent() for z in obj[1:]):
                K = obj[0].parent()
                base = cls.prep(K, make_int, escape_backslashes)
                return {'__NFList__': 0, # encoding version
                        'base': base,
                        'data': [cls.prep(z, make_int, escape_backslashes)['data'] for z in obj]}
            else:
                return [cls.prep(x, make_int, escape_backslashes) for x in obj]
        elif isinstance(obj, dict):
            if obj and all(isinstance(k, (int, Integer)) for k in obj):
                return {'__IntDict__': 0, # encoding version
                        'data': [[int(k), cls.prep(v, make_int, escape_backslashes)] for k, v in obj.items()]}
            elif all(isinstance(k, basestring) for k in obj):
                return {k:cls.prep(v, make_int, escape_backslashes) for k,v in obj.iteritems()}
            else:
                raise TypeError("keys must be strings or integers")
        elif isinstance(obj, FreeModuleElement):
            return {'__Vector__': 0, # encoding version
                    'base': cls.prep(obj.base_ring(), make_int, escape_backslashes),
                    'data': [cls.prep(c, make_int, escape_backslashes)['data'] for c in obj]}
        elif isinstance(obj, Integer):
            return int(obj)
        elif isinstance(obj, Rational):
            return {'__Rational__': 0, # encoding version
                    'data': [int(obj.numerator()), int(obj.denominator())]}
        elif isinstance(obj, RealNumber):
            return {'__RealLiteral__': 0, # encoding version
                    'data': obj.literal if isinstance(obj, RealLiteral) else str(obj),
                    'prec': int(obj.parent().precision())}
        elif isinstance(obj, complex):
            # As noted above, support for Sage complex numbers
            # would be easy to add
            return {'__complex__': 0, # encoding version
                    'data': [obj.real, obj.imag]}
        elif isinstance(obj, ComplexNumber):
            return {'__Complex__': 0, # encoding version
                    'prec': int(obj.prec()),
                    'data': [str(obj.real()), str(obj.imag())]}
        elif isinstance(obj, NumberFieldElement):
            return {'__NFElt__': 0, # encoding version
                    'parent': cls.prep(obj.parent(), make_int, escape_backslashes),
                    'data': [cls.prep(c, make_int, escape_backslashes)['data'] for c in obj.list()]}
        elif isinstance(obj, NumberField_generic):
            if isinstance(obj, NumberField_relative):
                return {'__NFRelative__': 0, # encoding version
                        'vname': obj.variable_name(),
                        'data': cls.prep(obj.relative_polynomial(), make_int, escape_backslashes)}
            elif isinstance(obj, NumberField_cyclotomic):
                return {'__NFCyclotomic__': 0, # encoding version
                        'data': int(obj._n())}
            else:
                return {'__NFAbsolute__': 0, # encoding version
                        'vname': obj.variable_name(),
                        'data': cls.prep(obj.absolute_polynomial(), make_int, escape_backslashes)}
        elif obj is ZZ:
            return {'__IntegerRing__': 0, # encoding version
                    'data': 0} # must be present for decoding
        elif obj is QQ:
            return {'__RationalField__': 0, # encoding version
                    'data': 0} # must be present for decoding
        elif isinstance(obj, Polynomial):
            return {'__Poly__': 0, # encoding version
                    'vname': obj.variable_name(),
                    'base': cls.prep(obj.base_ring(), make_int, escape_backslashes),
                    'data': [cls.prep(c, make_int, escape_backslashes)['data'] for c in obj.list()]}
        elif isinstance(obj, PowerSeries_poly):
            if obj.base_ring() is ZZ:
                data = [int(c) for c in obj.list()]
            else:
                data = [cls.prep(c, make_int, escape_backslashes)['data'] for c in obj.list()]
            return {'__PowerSeries__': 0, # encoding version
                    'vname': obj.variable(),
                    'base': cls.prep(obj.base_ring(), make_int, escape_backslashes),
                    'prec': 'inf' if obj.prec() is infinity else int(obj.prec()),
                    'data': data}
        elif escape_backslashes and isinstance(obj, basestring):
            # For use in copy_dumps below
            if all_digits.match(obj):
                return int(obj)
            else:
                return obj.replace('\\','\\\\\\\\').replace("\r", r"\r").replace("\n", r"\n").replace("\t", r"\t").replace('"',r'\"')
        elif obj is None:
            return None
        elif isinstance(obj, date):
            return {'__date__': 0,
                    'data': "%s"%(obj)}
        elif isinstance(obj, time):
            return {'__time__': 0,
                    'data': "%s"%(obj)}
        elif isinstance(obj, datetime):
            return {'__datetime__': 0,
                    'data': "%s"%(obj)}
        elif isinstance(obj, basestring):
            if all_digits.match(obj):
                return int(obj)
            return obj
        elif isinstance(obj, (basestring, int, long, bool, float)):
            return obj
        else:
            raise ValueError("Unsupported type: %s"%(type(obj)))
""")
            for cid in sorted(collection_names.keys(), key=lambda x: collection_names[x]):
                cname = collection_names[cid]
                if cname in skip_collections:
                    continue
                assert cname in sorts, "cname not in sorts:\ncname = %s\nsorts = %s" % (cname, sorts)
                uname = renames[cname]
                unames.append(uname)
                cols = types[cname]
                typecols = {}
                # The keys correspond roughly to postgres types, and the values are pairs, the first being a string listing the columns of that type and the second a code indicating how we'll be handling it.
                typecols['int'] = _build_typecol(cols, ['integer'], False, None, None)
                typecols['real'] = _build_typecol(cols, ['real'], True, 'real', 'numeric')
                typecols['int_as_str'] = _build_typecol(cols, ['integer stored as string'], True, 'int', 'numeric')
                typecols['real_as_str'] = _build_typecol(cols, ['real stored as string'], True, 'real', 'numeric')
                typecols['bool'] = _build_typecol(cols, ['boolean'], True, 'bool', 'boolean')
                typecols['str'] = _build_typecol(cols, ['string'], True, 'str', 'text')
                # "tate" typo is in inventory database
                typecols['datetime'] = _build_typecol(cols, ['tate and or time'], True, 'datetime', 'timestamp')
                typecols['list_as_str'] = _build_typecol(cols, ['comma separated list of integer stored as string'], True, 'split', 'jsonb')
                typecols['spaced_str'] = _build_typecol(cols, ['space separated list of integer stored as string'], True, 'split', 'jsonb')
                typecols['gal_pair'] = _build_typecol(cols, ['galois pair'], True, 'galpair', 'integer')
                typecols['unit_complex'] = _build_typecol(cols, ['unit complex'], True, 'angle', 'numeric')
                typecols['mixlist_as_str'] = _build_typecol(cols, ['comma separated list of mixed types stored as string'], True, 'asis', 'text')
                typecols['json'] = _build_typecol(cols, None, True, 'json', 'jsonb', collection=True)
                typecols['ijson'] = _build_typecol(cols, ['ijson'], True, 'ijson', 'jsonb')
                typecols['binary'] = _build_typecol(cols, ["non-primitive type (<class 'bson.binary.Binary'>)"], True, 'binary', 'bytea')
                typecols['link'] = _build_typecol(cols, ["non-primitive type (<class 'bson.objectid.ObjectId'>)"], False, "Create link to another table", None)
                typecols['missing'] = _build_typecol(cols, ["Record cannot be found containing key"], False, "Record not found so cannot determine type", None)
                wraps = []
                innertypes = []
                for typ, coldata in typecols.iteritems():
                    columns, normal, wrapper, postgres_type = coldata
                    if normal and columns:
                        wraps.append("""                for fld in [{0}]: # {1} columns
                    flds[fld] = {2}wrap(rec.get(fld))""".format(columns, typ, wrapper))
                        innertypes.append("""
        for fld in [{0}]:
            types['{1}'].append(fld)""".format(columns, postgres_type))
                    elif not normal and columns and wrapper:
                        # wrapper gives text for a NotImplementedError
                        wraps.append("""                for fld in [{0}]: # {1} columns: fix
                    raise NotImplementedError("{2}")""".format(columns, typ, wrapper))
                        innertypes.append("""
        for fld in [{}]:
            types['FIX'].append(fld)""".format(columns))
                    elif not normal and columns: # integer type
                        wraps.append("""                for fld in [{}]: # integer columns
                    maxval_update(flds, maxvals, fld, rec)""".format(columns))
                        innertypes.append("""
        for fld in [{}]:
            types[integer_size(maxvals[fld])].append(fld)""".format(columns))
                these_indexes = []
                if sorts[cname] is not None:
                    these_indexes.append("db.{collection}.restore_pkeys()".format(collection=uname))
                for index_name, index_keys in indexes[cid]:
                    if ['sort',1] in index_keys:
                        sort_spot = index_keys.index(['sort',1])
                        index_keys[sort_spot] = ['id',1]
                    if index_keys[0][0] not in ["_id", "metadata"]:
                        for key, direction in index_keys:
                            if key != 'id' and key not in cols:
                                raise RuntimeError("%s (%s): %s not in %s"%(index_name, ', '.join(k for k,d in index_keys), key, cname))
                        index_name = index_name.replace("_1","").replace("_-1","").replace("_sort","").replace("_key","")
                        if index_keys and index_keys[-1][0] == "gin":
                            typ = "gin"
                            index_keys = index_keys[:-1]
                            has_2d = has_json = False
                        else:
                            typ = "btree"
                            has_2d = any(direction == '2d' for col, direction in index_keys)
                            has_json = any(cols[col].startswith("collection of") for col, direction in index_keys if col in cols)
                        s = "db.{collection}.create_index({columns}, '{type}')".format(collection=uname, columns = [col for col, dummy in index_keys], type=typ)
                        if has_2d:
                            s += " # 2d key"
                        if has_json:
                            s += " # json col: may want to use GIN"
                        these_indexes.append(s)
                Fexport.write(r"""
def export_{collection}():
    conn = lmfdb.base.getDBConnection()
    {collection} = conn.{dotted_coll}
    maxvals = defaultdict(int)
    try:
        ordered_cols = [{ordered_cols}]
        with open('exports/{collection}.txt', 'w') as Fout:
            for i, rec in sort_collection({collection}, {export_sort}, '{dotted_coll}'):
                flds = {{}}
{wraps}
                flds["id"] = str(i)
                Fout.write("\t".join(flds[fld] for fld in ["id"] + ordered_cols) + "\n")
        types = defaultdict(list)
{typesets}
        with open('{importfile}', 'a') as Fimp:
            Fimp.write('''
def import_{collection}():
    try:
        db.create_table('{collection}', {{0}}, {label_col}, {import_sort}, {id_ordered}, search_order={{1}})
    except Exception:
        print "Failure in creating {collection}"
        traceback.print_exc()
        return
    try:
        db.{collection}.copy_from('/scratch/importing/{collection}.txt', search_cols={{1}}, includes_ids=True, resort=False)
    except Exception:
        print "Failure loading data into {collection}"
        traceback.print_exc()
        return
    print "Successfully imported {collection}"
    index_{collection}()
def index_{collection}():
    try:
        print "Indexing {collection}"
        {indexes}
    except Exception:
        print "Failure in indexing {collection}"
        traceback.print_exc()
    else:
        print "Successfully indexed {collection}"
'''.format(dict(types), ordered_cols))
    except Exception:
        print "Failure in exporting {collection}"
        traceback.print_exc()
        with open("export_failures.txt",'a') as F:
            F.write("{collection}\n")
    else:
        print "Successfully exported {collection}"
""".format(collection = uname,
           dotted_coll = cname,
           label_col = "None" if label_cols[cname] is None else "'" + label_cols[cname] + "'",
           importfile = importfile,
           import_sort = str(sorts[cname]),
           export_sort = str(sorts[cname]),
           wraps = "\n".join(wraps),
           typesets = "".join(innertypes),
           indexes = "\n        ".join(these_indexes),
           ordered_cols = ', '.join('"{}"'.format(col) for col in cols),
           json_cols = typecols['json'][0],
           id_ordered = (uname in id_ordered),
           ))

            Fexport.write(r"""
def export_all():
    with open("export_failures.txt",'w') as F:
        pass
    with open('{importfile}', 'w') as Fimp:
        Fimp.write('''
from psycopg2 import DatabaseError, connect
import traceback
from collections import defaultdict
from lmfdb.db_backend import db
def import_all():
{importlist}
''')
{exportlist}
""".format(importfile = importfile,
           importlist = "\n".join("    import_{}()".format(uname) for uname in unames),
           exportlist = "\n".join("    export_{}()".format(uname) for uname in unames)))


    # t0 = time.time()
    # maxsize = defaultdict(lambda: defaultdict(int))
    # mantissa = defaultdict(lambda: defaultdict(int))
    # floor = defaultdict(lambda: defaultdict(int))
    # realtypes = set()
    # bads = []
    # badD = defaultdict(set)
    # for cid, collection in collections.iteritems():
    #     cname = collection_names[cid]
    #     print cname, time.time() - t0
    #     cols = types[cname]
    #     intcols = [col for col in cols if cols[col] in ["integer", "integer stored as string"]]
    #     realcols = [col for col in cols if cols[col] == "real"]
    #     realstrcols = [col for col in cols if cols[col] == "real stored as string"]
    #     if intcols or realcols:
    #         cursor = collection.find(no_cursor_timeout=True)
    #         try:
    #             for rec in cursor:
    #                 for col in intcols:
    #                     if col in rec:
    #                         try:
    #                             if abs(int(rec[col])) > maxsize[cname][col]:
    #                                 maxsize[cname][col] = abs(int(rec[col]))
    #                         except ValueError:
    #                             bads.append((cname, col, rec[col], "istring"))
    #                 for col in realcols:
    #                     if col in rec:
    #                         realtypes.add(type(rec[col]))
    #                 for col in realstrcols:
    #                     if col in rec:
    #                         if not isinstance(rec[col], basestring):
    #                             bads.append((cname, col, rec[col], "type"))
    #                             badD[cname].add(col)
    #                             continue
    #                         realmatch = realre.match(rec[col])
    #                         if not realmatch:
    #                             bads.append((cname, col, rec[col], "rstring"))
    #                             badD[cname].add(col)
    #                             continue
    #                         floor[cname][col] = max(floor[cname][col],len(realmatch.group(1)))
    #                         mantissa[cname][col] = max(mantissa[cname][col],len(realmatch.group(2))-1)
    #         finally:
    #             cursor.close()
    #         for col in maxsize[cname]:
    #             print "IntMax {0}: {1}".format(col, maxsize[cname][col])
    #         for col in floor[cname]:
    #             print "RealSize {0}: {1}.{2}".format(col, floor[cname][col], mantissa[cname][col])
    #         for col in badD[cname]:
    #             print "Bad {0}".format(col)
    # return maxsize, floor, mantissa, realtypes, bads, types, examples
