def addstats_all():
    failures = []
    for f in [lf_fields_addstats,
              ec_curves_addstats,
              hmf_forms_addstats,
              lat_lattices_addstats,
              hecke_algebras_addstats,
              hgcwa_passports_addstats,
              g2c_curves_addstats,
              bmf_forms_addstats,
              ec_nfcurves_addstats,
              av_fqisog_addstats,
              artin_reps_addstats,
              hgm_families_addstats,
              nf_fields_addstats]:
        try:
            f()
        except Exception:
            failures.append(f)
    if failures:
        print "%s failures: %s"%(len(failures), failures)
    else:
        print "No failures!"

def ec_curves_addstats():
    cols = ['jinv','torsion','rank','sha','torsion_structure']
    constraints = [None,
                   {'nonmax_primes':{'$contains':[2]}},
                   {'nonmax_primes':{'$contains':[3]}},
                   {'nonmax_primes':{'$contains':[5]}},
                   {'nonmax_primes':{'$contains':[7]}},
                   {'isogeny_degrees':{'$contains':[2]}},
                   {'isogeny_degrees':{'$contains':[3]}},
                   {'isogeny_degrees':{'$contains':[4]}},
                   {'isogeny_degrees':{'$contains':[5]}},
                   {'cm':0},
                   {'cm':{'$ne':0}},
                   {'number':1}]
    db.ec_curves.stats._clear_stats_counts()
    db.ec_curves.stats.add_stats_auto(cols, constraints, max_depth=3)

def ec_nfcurves_addstats():
    cols = ['torsion_order', 'torsion_structure', 'field_label']
    constraints = [None,
                   {'isogeny_degrees':{'$contains':[2]}},
                   {'isogeny_degrees':{'$contains':[3]}},
                   {'isogeny_degrees':{'$contains':[4]}},
                   {'isogeny_degrees':{'$contains':[5]}},
                   {'cm':0},
                   {'cm':{'$ne':0}},
                   {'number':1},
                   {'base_change':[]}]
    db.ec_nfcurves.stats._clear_stats_counts()
    db.ec_nfcurves.stats.add_stats_auto(cols, constraints)

def nf_fields_addstats():
    cols = ['degree', 'r2', 'galt', 'class_number', 'class_group', 'disc_rad']
    constraints = [{'ramps':{'$contains':[2]}},
                   {'ramps':{'$contains':[3]}},
                   {'ramps':{'$contains':[5]}},
                   {'ramps':{'$contains':[7]}},
                   {'ramps':{'$contains':[11]}},
                   {'ramps':{'$contains':[13]}},
                   {'ramps':{'$contains':[2,3]}},
                   {'ramps':{'$contains':[2,5]}},
                   {'ramps':{'$contains':[2,7]}},
                   {'ramps':{'$contains':[3,5]}},
                   {'ramps':{'$contains':[3,7]}},
                   {'ramps':{'$contains':[5,7]}},
                   {'ramps':{'$contains':[2,3,5]}}]
    db.nf_fields.stats._clear_stats_counts()
    db.nf_fields.stats.add_stats_auto(cols, [None], max_depth=3)
    db.nf_fields.stats.add_stats_auto(cols[:-1], constraints, max_depth=3)

def av_fqisog_addstats():
    cols = ['g', 'q', 'p_rank', 'ang_rank', ['galois_n', 'galois_t'], 'is_simp', 'is_prim', 'is_jac', 'is_pp']
    db.av_fqisog.stats._clear_stats_counts()
    db.av_fqisog.stats.add_stats_auto(cols, max_depth=3)

def hmf_forms_addstats():
    cols = ['field_label', 'deg', 'weight', 'level_norm', 'dimension', 'is_CM', 'is_base_change']
    db.hmf_forms.stats._clear_stats_counts()
    db.hmf_forms.stats.add_stats_auto(cols, max_depth=3)

def lat_lattices_addstats():
    cols = ['dim', 'det', 'level', 'minimum', 'class_number', 'aut']
    db.lat_lattices.stats._clear_stats_counts()
    db.lat_lattices.stats.add_stats_auto(cols, max_depth=3)

def hecke_algebras_addstats():
    cols = ['level', 'weight', 'num_orbits', 'ell']
    db.hecke_algebras.stats._clear_stats_counts()
    db.hecke_algebras.stats.add_stats_auto(cols)

def hgcwa_passports_addstats():
    cols = ['group', 'genus', 'signature', 'dim', 'hyperelliptic', 'cyclic_trigonal']
    constraints = [None,
                   {'full_auto':{'$exists':True}},
                   {'full_auto':{'$exists':False}}]
    db.hgcwa_passports.stats._clear_stats_counts()
    db.hgcwa_passports.stats.add_stats_auto(cols, constraints, max_depth=3)

def g2c_curves_addstats():
    cols = ['num_rat_pts', 'num_rat_wpts', 'aut_grp_id', 'geom_aut_grp_id', 'analytic_rank', 'two_selmer_rank', 'has_square_sha', 'locally_solvable', 'is_gl2_type', 'real_geom_end_alg', 'st_group', 'torsion_order']
    db.g2c_curves.stats._clear_stats_counts()
    for col in cols:
        db.g2c_curves.stats.add_stats([col])
    cols.extend(['is_simple_geom','torsion_subgroup'])
    db.g2c_curves.stats.add_stats_auto(cols, max_depth=3)

def lf_fields_addstats():
    cols = ['n', 'p', 'c', 'e', 'galT', 'topslope']
    db.lf_fields.stats._clear_stats_counts()
    db.lf_fields.stats.add_stats_auto(cols, max_depth=3)

def bmf_forms_addstats():
    cols = ['field_label', 'weight', 'level_norm', 'dimension']
    db.bmf_forms.stats._clear_stats_counts()
    db.bmf_forms.stats.add_stats_auto(cols, max_depth=3)

def artin_reps_addstats():
    cols = ['Indicator', 'Container', ['Galn', 'Galt'], 'Dim']
    constraints = [None,
                   {'GalConjSigns':{'$contains':1}},
                   {'GalConjSigns':{'$contains':-1}},
                   {'BadPrimes':{'$contains':[2]}},
                   {'BadPrimes':{'$contains':[3]}},
                   {'BadPrimes':{'$contains':[5]}},
                   {'BadPrimes':{'$contains':[7]}},
                   {'BadPrimes':{'$notcontains':[2]}},
                   {'BadPrimes':{'$notcontains':[3]}},
                   {'BadPrimes':{'$notcontains':[5]}},
                   {'BadPrimes':{'$notcontains':[7]}}]
    db.artin_reps.stats._clear_stats_counts()
    db.artin_reps.stats.add_stats_auto(cols, contstaints, max_depth=3)

def gps_transitive_addstats():
    pass

def gps_sato_tate_addstats():
    pass

def hgm_motives_addstats():
    pass

def hgm_families_addstats():
    cols = ['degree', 'weight']
    db.hgm_families.stats._clear_stats_counts()
    db.hgm_families.stats.add_stats_auto(cols)

def fq_fields_addstats():
    pass

def bmf_dims_addstats():
    pass

def belyi_galmaps_addstats():
    pass

def belyi_passports_addstats():
    pass

def g2c_endomorphisms_addstats():
    pass
def g2c_ratpts_addstats():
    pass
def g2c_tamagawa_addstats():
    pass
def halfmf_forms_addstats():
    pass
def gps_small_addstats():
    pass
def hecke_orbits_addstats():
    pass
def hmf_fields_addstats():
    pass
def hmf_hecke_addstats():
    pass
def artin_field_data_addstats():
    pass
def ec_padic_addstats():
    pass
def smf_dims_addstats():
    pass
def hecke_ladic_addstats():
    pass
def gps_gmodules_addstats():
    pass
def ec_iqf_labels_addstats():
    pass
def smf_families_addstats():
    pass
def mwf_coeffs_addstats():
    pass
def lfunc_instances_addstats():
    pass
def mwf_plots_addstats():
    pass
def mwf_tables_addstats():
    pass
def mwf_forms_addstats():
    pass
def smf_samples_addstats():
    pass
def smf_ev_addstats():
    pass
def smf_fc_addstats():
    pass
def gps_sato_tate0_addstats():
    pass
def lfunc_zeros_addstats():
    pass
def modlgal_reps_addstats():
    pass
def modlmf_forms_addstats():
    pass
def mwfp_forms_addstats():
    pass
def sl2z_subgroups_addstats():
    pass
def lfunc_lfunctions_addstats():
    pass
