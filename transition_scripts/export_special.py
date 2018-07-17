from collections import defaultdict
import traceback
import lmfdb
from datetime import datetime, timedelta, date, time
from sage.libs.pari import pari
import json
import gridfs
from gridfs.errors import CorruptGridFile
import re
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
from sage.structure.sage_object import loads
from lmfdb.modular_forms.maass_forms.maass_waveforms.backend.maass_forms_db import MaassDB
stripzero_re = re.compile(r'^([0-9-]+)(\.0*)?$')
all_digits = re.compile(r'^-?[0-9]+$')
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
def binarywrap(n):
    return r'\N' if n is None else r'\\x' + ''.join(c.encode('hex') for c in n)
def nbrwrap(n):
    return r'\N' if n is None else str(n).replace(" ","")
def asiswrap(n):
    return r'\N' if n is None else n
def strwrap(n):
    return ur'\N' if n is None else n.replace(u'\\', u'\\\\').replace(u"\r", ur"\r").replace(u"\n", ur"\n").replace(u"\t", ur"\t").replace(u'"',ur'\"')
def datetimewrap(n):
    return r'\N' if not n else str(n)
def boolwrap(b):
    if b is None:
        return r'\N'
    elif b:
        return 't'
    else:
        return 'f'
def jsonwrap(n):
    return r'\N' if n is None else Json.dumps(n)
def ijsonwrap(n):
    return r'\N' if n is None else Json.dumps(n, make_int=True)
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
def prod(iterable):
    return reduce(operator.mul, iterable, 1)
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
def report_time(i, total, name, t0, interval=10000):
    if i and i % interval == 0:
        t = datetime.now()
        print "%s: %s/%s in %s"%(name, i, total, t-t0)
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
            for i, rec in enumerate(coll.find({},sort,no_cursor_timeout=True)):
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
def integer_size(n):
    if n < 2**7:
        return 'smallint'
    elif n < 2**15:
        return 'integer'
    elif n < 2**31:
        return 'bigint'
    else:
        return 'numeric'
class Json(object):
    @classmethod
    def dumps(cls, obj, make_int=False):
        return json.dumps(cls.prep(obj, make_int=make_int))

    @classmethod
    def prep(cls, obj, make_int=False, escape_backslashes=True):
        """
        Returns a version of the object that is parsable by the standard json dumps function.
        For example, replace Integers with ints, encode various Sage types using dictionaries....
        """
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
            if make_int and all_digits.match(obj):
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
            if make_int and all_digits.match(obj):
                return int(obj)
            return obj
        elif isinstance(obj, (basestring, int, long, bool, float)):
            return obj
        else:
            raise ValueError("Unsupported type: %s"%(type(obj)))

def export_mwfp_forms():
    conn = lmfdb.base.getDBConnection()
    mwfp_forms = conn.HTPicard.picard
    maxvals = defaultdict(int)
    try:
        ordered_cols = ["maass_id", "ev", "prec", "sym", "coef"]
        with open('exports/mwfp_forms.txt', 'w') as Fout:
            for i, rec in sort_collection(mwfp_forms, None, 'HTPicard.picard'):
                flds = {}

                flds["maass_id"] = str(rec["_id"])
                for fld in ["ev", "prec"]: # real columns
                    flds[fld] = realwrap(rec.get(fld))
                for fld in ["sym", "coef"]: # json columns
                    flds[fld] = jsonwrap(rec.get(fld))
                if i is not None:
                    flds["id"] = str(i)
                Fout.write("\t".join(flds[fld] for fld in ordered_cols) + "\n")
        types = defaultdict(list)
        for fld in ["ev", "prec"]:
            types['numeric'].append(fld)
        for fld in ["sym", "coef"]:
            types['jsonb'].append(fld)
        types['text'] = ["maass_id"]
        with open('import_special.py', 'a') as Fimp:
            Fimp.write('''
def import_mwfp_forms():
    try:
        db.create_table('mwfp_forms', {0}, 'maass_id', ['ev'], False, search_order={1})
    except Exception:
        print "Failure in creating mwfp_forms"
        traceback.print_exc()
        return
    try:
                db.mwfp_forms.copy_from('/scratch/importing/mwfp_forms.txt', search_cols={1}, includes_ids=True, resort=False)
    except Exception:
        print "Failure loading data into mwfp_forms"
        traceback.print_exc()
        return
    print "Successfully imported mwfp_forms"
    index_mwfp_forms()
def index_mwfp_forms():
    try:
        print "Indexing mwfp_forms"
        db.mwfp_forms.restore_pkeys()
        db.mwfp_forms.create_index([u'ev'], 'btree')
    except Exception:
        print "Failure in indexing mwfp_forms"
        traceback.print_exc()
    else:
        print "Successfully indexed mwfp_forms"
'''.format(dict(types), ordered_cols))
    except Exception:
        print "Failure in exporting mwfp_forms"
        traceback.print_exc()
        with open("export_failures.txt",'a') as F:
            F.write("mwfp_forms\n")
    else:
        print "Successfully exported mwfp_forms"

def export_g2c_instances():
    conn = lmfdb.base.getDBConnection()
    instances = conn.genus2_curves.Lfunctions.instances
    try:
        with open("exports/g2c_instances.txt", 'w') as Fout:
            ordered_cols = ["Lhash","type","url"]
            for instance in instances.find():
                flds = {}
                for fld in ordered_cols:
                    flds[fld] = asiswrap(instance.get(fld))
                if "urls" in instance:
                    assert "url" not in instance
                    flds["url"] = u"|".join(instance.get("urls"))
                Fout.write((u"\t".join(flds[fld] for fld in ordered_cols) + u"\n").encode("utf-8"))
        with open("import_special.py", 'a') as Fimp:
            Fimp.write('''
def import_g2c_instances():
    try:
        db.lfunc_instances.copy_from('/scratch/importing/g2c_instances.txt', search_cols=['Lhash', 'type', 'url'])
        conn.commit()
    except DatabaseError as err:
        conn.rollback()
        print "Failure in importing g2c_instances"
        print err
    else:
        print "Successfully imported g2c_instances"
''')
    except Exception as err:
        print "Failure in exporting g2c_instances"
        traceback.print_exc()
    else:
        print "Successfully exported g2c_instances"

def export_users():
    conn = lmfdb.base.getDBConnection()
    users = conn.userdb.users
    try:
        with open("exports/users.txt", 'w') as Fout:
            ordered_cols = ["username","password","bcpassword","admin","color_scheme","full_name","email","url","about","created"]
            for user in users.find():
                flds = {}
                flds["username"] = user["_id"]
                for fld in ["about","full_name","password","bcpassword","url"]:
                    if user.get(fld):
                        flds[fld] = strwrap(user.get(fld))
                    else:
                        flds[fld] = ur"\N"
                flds["email"] = r'\N' # nobody has e-mail set correctly
                flds["admin"] = boolwrap(user.get("admin",False))
                flds["color_scheme"] = nbrwrap(user.get("color_scheme"))
                flds["created"] = datetimewrap(user.get("created"))
                Fout.write((u"\t".join(flds[fld] for fld in ordered_cols) + u"\n").encode("utf-8"))
        with open("import_special.py", 'a') as Fimp:
            Fimp.write('''
def import_users():
    try:
        conn = db.conn
        cur = conn.cursor()
        cur.execute('CREATE SCHEMA userdb;')
        cur.execute('CREATE TABLE userdb.users (username text COLLATE "C", password text COLLATE "C", bcpassword text COLLATE "C", admin boolean, color_scheme smallint, full_name text COLLATE "C", email text COLLATE "C", url text COLLATE "C", about text COLLATE "C", created timestamp);')
        with open('/scratch/importing/users.txt') as F:
            cur.copy_from(F, 'userdb.users', columns=["username", "password", "bcpassword", "admin", "color_scheme", "full_name", "email", "url", "about", "created"])
        conn.commit()
    except DatabaseError as err:
        conn.rollback()
        print "Failure in importing users"
        print err
    else:
        print "Successfully imported users"
''')
    except Exception as err:
        print "Failure in exporting users"
        traceback.print_exc()
    else:
        print "Successfully exported users"

def export_knowls():
    conn = lmfdb.base.getDBConnection()
    knowls = conn.knowledge.knowls
    with open('exports/kwl_knowls.txt', 'w') as F:
        for knowl in knowls.find():
            flds = {}
            # unchanged string fields
            for fld in ["_id", "cat", "title", "content", "last_author", "quality"]:
                flds[fld] = strwrap(knowl.get(fld))
            for fld in ["_keywords", "authors"]:
                flds[fld] = jsonwrap(knowl.get(fld))
            history = knowl.get("history")
            if history is None:
                flds["history"] = r'\N'
            else:
                history = [hist for hist in history if hist is not None]
                for i, hist in enumerate(history):
                    hist['id'] = hist.pop('_id')
                    if 'timestamp' in hist:
                        hist['timestamp'] = str(hist['timestamp'])
                    for key, value in hist.iteritems():
                        if isinstance(value, basestring):
                            hist[key] = value.replace('\\', '\\\\\\\\').replace("\r", r"\r").replace("\n", r"\n").replace("\t", r"\t").replace('"',r'\"')
                    history[i] = hist
                flds["history"] = json.dumps(history)
            flds["timestamp"] = str(knowl.get("timestamp",r"\N"))
            F.write((u"\t".join(flds[fld] for fld in ["_id", "cat", "title", "content", "authors", "last_author", "quality", "timestamp", "_keywords", "history"]) + u"\n").encode('utf-8'))
    history_db = conn.knowledge.history
    with open('kwl_history.txt', 'w') as F:
        for history in history_db.find():
            flds = {}
            for fld in ["_id", "title", "who", "state"]:
                flds[fld] = asiswrap(history.get(fld))
            flds["time"] = str(history.get("time", r"\N"))
            F.write("\t".join(flds[fld] for fld in ["_id", "title", "time", "who", "state"]) + "\n")
    with open('import_special.py', 'a') as F:
        F.write('''
def import_knowls():
    cur = db.conn.cursor()
    try:
        cur.execute("CREATE TABLE kwl_knowls (id text, cat text, title text, content text, authors jsonb, last_author text, quality text, timestamp timestamp, _keywords jsonb, history jsonb)")
        cur.execute("CREATE TABLE kwl_deleted (id text, cat text, title text, content text, authors jsonb, last_author text, quality text, timestamp timestamp, _keywords jsonb, history jsonb)")
        cur.execute("CREATE TABLE kwl_history (id text, title text, time timestamp, who text, state text)")
        for tbl in ["kwl_knowls", "kwl_deleted", "kwl_history"]:
            for action in ["INSERT", "UPDATE", "DELETE"]:
                db._grant(action, tbl, ['webserver'], False)
            db.grant_select(tbl, commit=False)
        with open('/scratch/importing/kwl_knowls.txt') as F:
            cur.copy_from(F, 'kwl_knowls', columns=["id", "cat", "title", "content", "authors", "last_author", "quality", "timestamp", "_keywords", "history"])
        with open('/scratch/importing/kwl_history.txt') as F:
            cur.copy_from(F, 'kwl_history', columns=["id", "title", "time", "who", "state"])
        cur.execute("ALTER TABLE kwl_knowls ADD CONSTRAINT kwl_knowls_pkey PRIMARY KEY (id)")
	# no primary key on deleted
        #cur.execute("ALTER TABLE kwl_deleted ADD CONSTRAINT kwl_deleted_pkey PRIMARY KEY (id)")
        cur.execute("ALTER TABLE kwl_history ADD CONSTRAINT kwl_history_pkey PRIMARY KEY (id)")
    except Exception:
        print "Failure in importing knowls"
        traceback.print_exc()
        db.conn.rollback()
    else:
        db.conn.commit()
        print "Succeeded in importing knowls"
''')

def export_g2c_curves():
    conn = lmfdb.base.getDBConnection()
    g2c_curves = conn.genus2_curves.curves
    maxvals = defaultdict(int)
    try:
        ordered_cols = ["geom_aut_grp_id", "igusa_inv", "num_rat_wpts", "is_simple_base", "is_simple_geom", "torsion_order", "g2_inv", "cond", "Lhash", "abs_disc", "has_square_sha", "two_selmer_rank", "analytic_rank", "label", "st_group", "aut_grp_id", "eqn", "num_rat_pts", "locally_solvable", "bad_lfactors", "is_gl2_type", "non_solvable_places", "disc_sign", "class", "igusa_clebsch_inv", "globally_solvable", "real_geom_end_alg", "torsion_subgroup", "root_number", "two_torsion_field"]
        with open('exports/g2c_curves.txt', 'w') as Fout:
            for i, rec in sort_collection(g2c_curves, ['cond', 'class', 'abs_disc', 'disc_sign', 'label'], 'genus2_curves.curves'):
                flds = {}
                for fld in ["num_rat_wpts", "torsion_order", "cond", "abs_disc", "two_selmer_rank", "analytic_rank", "num_rat_pts", "disc_sign", "globally_solvable", "root_number"]: # integer columns
                    maxval_update(flds, maxvals, fld, rec)
                for fld in ["non_solvable_places", "two_torsion_field"]: # json columns
                    flds[fld] = jsonwrap(rec.get(fld))
                for fld in ["is_simple_base", "is_simple_geom", "has_square_sha", "locally_solvable", "is_gl2_type"]: # bool columns
                    flds[fld] = boolwrap(rec.get(fld))
                for fld in ["geom_aut_grp_id", "igusa_inv", "g2_inv", "Lhash", "label", "st_group", "aut_grp_id", "bad_lfactors", "class", "igusa_clebsch_inv", "real_geom_end_alg", "torsion_subgroup"]: # str columns
                    flds[fld] = strwrap(rec.get(fld))
                for fld in ["eqn"]: # mixlist_as_str columns
                    flds[fld] = asiswrap(rec.get(fld))
                flds["id"] = str(i)
                Fout.write("\t".join(flds[fld] for fld in ["id"] + ordered_cols) + "\n")
        types = defaultdict(list)

        for fld in ["num_rat_wpts", "torsion_order", "cond", "abs_disc", "two_selmer_rank", "analytic_rank", "num_rat_pts", "disc_sign", "globally_solvable", "root_number"]:
            types[integer_size(maxvals[fld])].append(fld)
        for fld in ["non_solvable_places", "two_torsion_field"]:
            types['jsonb'].append(fld)
        for fld in ["is_simple_base", "is_simple_geom", "has_square_sha", "locally_solvable", "is_gl2_type"]:
            types['boolean'].append(fld)
        for fld in ["geom_aut_grp_id", "igusa_inv", "g2_inv", "Lhash", "label", "st_group", "aut_grp_id", "bad_lfactors", "class", "igusa_clebsch_inv", "real_geom_end_alg", "torsion_subgroup"]:
            types['text'].append(fld)
        for fld in ["eqn"]:
            types['text'].append(fld)
        with open('import_special.py', 'a') as Fimp:
            Fimp.write('''
def import_g2c_curves():
    try:
        db.create_table('g2c_curves', {0}, 'label', ['cond', 'class', 'abs_disc', 'disc_sign', 'label'], False, search_order={1})
    except Exception:
        print "Failure in creating g2c_curves"
        traceback.print_exc()
        return
    try:
        db.g2c_curves.copy_from('/scratch/importing/g2c_curves.txt', search_cols={1}, includes_ids=True, resort=False)
    except Exception:
        print "Failure loading data into g2c_curves"
        traceback.print_exc()
        return
    print "Successfully imported g2c_curves"
    index_g2c_curves()
def index_g2c_curves():
    try:
        print "Indexing g2c_curves"
        db.g2c_curves.restore_pkeys()
        db.g2c_curves.create_index([u'abs_disc'], 'btree')
        db.g2c_curves.create_index([u'analytic_rank'], 'btree')
        db.g2c_curves.create_index([u'aut_grp_id'], 'btree')
        db.g2c_curves.create_index([u'class'], 'btree')
        db.g2c_curves.create_index([u'cond'], 'btree')
        db.g2c_curves.create_index(['cond', 'class', 'abs_disc', 'label'], 'btree')
        db.g2c_curves.create_index([u'g2_inv'], 'btree')
        db.g2c_curves.create_index([u'geom_aut_grp_id'], 'btree')
        db.g2c_curves.create_index([u'has_square_sha'], 'btree')
        db.g2c_curves.create_index([u'is_gl2_type'], 'btree')
        db.g2c_curves.create_index([u'is_simple_geom'], 'btree')
        db.g2c_curves.create_index([u'label'], 'btree')
        db.g2c_curves.create_index([u'locally_solvable'], 'btree')
        db.g2c_curves.create_index([u'num_rat_pts'], 'btree')
        db.g2c_curves.create_index([u'real_geom_end_alg'], 'btree')
        db.g2c_curves.create_index([u'st_group'], 'btree')
        db.g2c_curves.create_index([u'torsion_order'], 'btree')
        db.g2c_curves.create_index([u'torsion_subgroup'], 'btree')
    except Exception:
        print "Failure in indexing g2c_curves"
        traceback.print_exc()
    else:
        print "Successfully indexed g2c_curves"
'''.format(dict(types), ordered_cols))
    except Exception:
        print "Failure in exporting g2c_curves"
        traceback.print_exc()
        with open("export_failures.txt",'a') as F:
            F.write("g2c_curves\n")
    else:
        print "Successfully exported g2c_curves"

def export_ec_curves():
    conn = lmfdb.base.getDBConnection()
    ec = conn.elliptic_curves.curves
    total = ec.count()
    t0 = datetime.now()
    maxvals = defaultdict(int)
    try:
        search_cols = ['label', 'lmfdb_label', 'iso', 'lmfdb_iso', 'iso_nlabel', 'number', 'lmfdb_number', 'ainvs', 'jinv', 'conductor', 'torsion', 'rank', 'sha', 'torsion_structure', 'cm', 'isogeny_degrees', 'nonmax_primes', 'nonmax_rad']
        extra_cols = ['equation', 'signD', 'torsion_generators', 'xcoord_integral_points', 'gens', 'heights', 'regulator', 'tamagawa_product', 'special_value', 'real_period', 'degree', 'modp_images', '2adic_label', '2adic_index', '2adic_log_level', '2adic_gens', 'isogeny_matrix', 'class_deg', 'class_size', 'sha_an', 'sha_primes', 'torsion_primes', 'tor_degs', 'tor_fields', 'tor_gro', 'local_data', 'min_quad_twist', 'aplist', 'anlist', 'iwdata', 'iwp0', 'galois_images']
        with open("exports/ec_curves.txt", 'w') as Fsmall:
            with open("exports/ec_curves_extras.txt", 'w') as Fextras:
                for i, rec in enumerate(ec.find().sort([('conductor', 1), ('iso_nlabel', 1), ('lmfdb_number', 1)])):
                    if i and i % 10000 == 0:
                        t = datetime.now()
                        print "%s/%s in %s"%(i, total, t-t0)
                    flds = {}
                    flds["id"] = str(i)
                    # unchanged string fields
                    for fld in ["label", "lmfdb_label", "iso", "lmfdb_iso", "jinv", "2adic_label", "equation"]:
                        flds[fld] = strwrap(rec.get(fld))
                    for fld in ["real_period", "special_value", "regulator", "sha_an"]: # real columns
                        flds[fld] = realwrap(rec.get(fld))
                    # integers
                    for fld in ["conductor", "iso_nlabel", "number", "lmfdb_number", "signD", "cm", "rank", "torsion", "tamagawa_product", "degree", "2adic_index", "2adic_log_level", "class_deg", "class_size", "sha", "iwp0"]:
                        maxval_update(flds, maxvals, fld, rec)
                    # lists of integers currently stored as strings
                    for new,old in [("ainvs", "xainvs"), ("xcoord_integral_points", "x-coordinates_of_integral_points")]:
                        if old in rec:
                            flds[new] = splitwrap(rec[old][1:-1])
                        else:
                            flds[new] = r'\N'
                    nonmax = rec.get("non-maximal_primes")
                    flds["nonmax_primes"] = jsonwrap(nonmax)
                    if nonmax is None:
                        flds["nonmax_rad"] = r'\N'
                    else:
                        flds["nonmax_rad"] = str(prod(nonmax))
                    for fld in ["torsion_structure", "isogeny_degrees", "sha_primes", "galois_images", "aplist", "iwdata", "min_quad_twist", "local_data", "tor_gro", "heights", "torsion_generators", "2adic_gens", "torsion_primes", "tor_degs", "isogeny_matrix", "tor_fields", "gens", "anlist"]: # json columns
                        flds[fld] = jsonwrap(rec.get(fld))
                    flds["modp_images"] = jsonwrap(rec.get("mod-p_images"))
                    Fsmall.write("\t".join(flds[fld] for fld in ["id"] + search_cols) + "\n")
                    Fextras.write("\t".join(flds[fld] for fld in ["id"] + extra_cols) + "\n")
        types = defaultdict(list)
        for fld in ["real_period", "special_value", "regulator", "sha_an"]:
            types['numeric'].append(fld)
        #for fld in ["degree", "cm", "number", "rank", "torsion", "iso_nlabel", "2adic_log_level", "class_deg", "lmfdb_number", "conductor", "iwp0", "signD", "tamagawa_product", "2adic_index", "sha", "class_size"]:
        #    types[integer_size(maxvals[fld])].append(fld)
        for fld in ["signD", "cm", "2adic_index", "2adic_log_level", "class_deg", "class_size", "number", "lmfdb_number", "iso_nlabel", "rank", "iwp0", "torsion"]:
            types['smallint'].append(fld)
        for fld in ["tamagawa_product", "sha", "nonmax_rad"]:
            types['integer'].append(fld)
        for fld in ["degree", "conductor"]:
            types['numeric'].append(fld)
        for fld in ["xcoord_integral_points", "nonmax_primes", "torsion_structure", "ainvs", "isogeny_degrees", "sha_primes", "modp_images", "galois_images", "aplist", "iwdata", "min_quad_twist", "local_data", "tor_gro", "heights", "torsion_generators", "2adic_gens", "torsion_primes", "tor_degs", "isogeny_matrix", "tor_fields", "gens", "anlist"]:
            types['jsonb'].append(fld)
        for fld in ["label", "lmfdb_iso", "2adic_label", "jinv", "equation", "iso", "lmfdb_label"]:
            types['text'].append(fld)
        search_types = defaultdict(list)
        extra_types = defaultdict(list)
        for typ in types:
            for col in types[typ]:
                if col in search_cols:
                    search_types[typ].append(col)
                else:
                    extra_types[typ].append(col)
        with open('import_special.py', 'a') as Fimp:
            Fimp.write('''
def import_ec_curves():
    try:
        db.create_table('ec_curves', {0}, 'lmfdb_label', ['conductor', 'iso_nlabel', 'lmfdb_number'], True, {1}, search_order={2}, extra_order={3})
    except Exception:
        print "Failure in creating ec_curves"
        traceback.print_exc()
        return
    try:
        db.ec_curves.copy_from('/scratch/importing/ec_curves.txt', '/scratch/importing/ec_curves_extras.txt', search_cols={2}, extra_cols={3}, includes_ids=True, resort=False)
    except Exception:
        print "Failure loading data into ec_curves"
        traceback.print_exc()
        return
    print "Successfully imported ec_curves"
    index_ec_curves()
def index_ec_curves():
    try:
        print "Indexing ec_curves"
        db.ec_curves.restore_pkeys()
        db.ec_curves.create_index(['ainvs'], 'btree')
        db.ec_curves.create_index(['cm'], 'btree')
        db.ec_curves.create_index(['conductor', 'iso_nlabel', 'lmfdb_number'], 'btree')
        db.ec_curves.create_index(['iso'], 'btree')
        db.ec_curves.create_index(['isogeny_degrees'], 'gin')
        db.ec_curves.create_index(['jinv', 'id'], 'btree')
        db.ec_curves.create_index(['label'], 'btree')
        db.ec_curves.create_index(['label', 'number'], 'btree')
        db.ec_curves.create_index(['lmfdb_label'], 'btree')
        db.ec_curves.create_index(['lmfdb_label', 'number'], 'btree')
        db.ec_curves.create_index(['lmfdb_iso'], 'btree')
        db.ec_curves.create_index(['lmfdb_number'], 'btree')
        db.ec_curves.create_index(['nonmax_primes'], 'gin')
        db.ec_curves.create_index(['number'], 'btree')
        db.ec_curves.create_index(['rank'], 'btree')
        db.ec_curves.create_index(['rank', 'number'], 'btree')
        db.ec_curves.create_index(['sha', 'id'], 'btree')
        db.ec_curves.create_index(['sha', 'rank', 'id'], 'btree')
        db.ec_curves.create_index(['sha', 'rank', 'torsion', 'id'], 'btree')
        db.ec_curves.create_index(['torsion'], 'btree')
        db.ec_curves.create_index(['torsion_structure'], 'btree')
    except Exception:
        print "Failure in indexing ec_curves"
        traceback.print_exc()
    else:
        print "Successfully indexed ec_curves"
'''.format(dict(search_types), dict(extra_types), search_cols, extra_cols))
    except Exception:
        print "Failure in exporting ec_curves"
        traceback.print_exc()
        with open("export_failures.txt",'a') as F:
            F.write("ec_curves\n")
    else:
        print "Successfully exported ec_curves"

def export_nf_fields():
    conn = lmfdb.base.getDBConnection()
    nf = conn.numberfields.fields
    total = nf.count()
    t0 = datetime.now()
    maxvals = defaultdict(int)
    pentuples = []
    for i, field in enumerate(nf.find({},['degree', 'disc_abs_key', 'disc_sign', 'label'])):
        if i and i%10000 == 0:
            t = datetime.now()
            print "nf_fields sort %s/%s in %s"%(i, total, t-t0)
        label = field['label']
        degree, sig, disc, num = map(int, label.split('.'))
        pentuples.append((degree, abs(disc), field['disc_sign'], num, label))
    pentuples.sort()
    id_dict = {}
    for i, (degree, disc_abs, disc_sign, num, label) in enumerate(pentuples):
        id_dict[label] = str(i)
    try:
        search_cols = ['label', 'coeffs', 'degree', 'r2', 'cm', 'iso_number', 'disc_abs', 'disc_sign', 'disc_rad', 'ramps', 'galt', 'class_number', 'class_group', 'used_grh']
        extra_cols = ['zk', 'units', 'reg', 'subs', 'unitsGmodule', 'unitsType', 'res', 'loc_algebras']
        with open('exports/nf_fields.txt', 'w') as Fsmall:
            #with open('exports/nf_fields_extras.txt', 'w') as Fextras:
            with open('exports/nf_fields_extras.txt', 'a') as Fextras:
                for i, rec in enumerate(nf.find()):
                    if i and i%10000 == 0:
                        t = datetime.now()
                        print "nf_fields %s/%s in %s"%(i, total, t-t0)
                    flds = {}
                    flds["id"] = id_dict[rec["label"]]
                    flds["disc_abs"] = rec["disc_abs_key"][3:]
                    for fld in ["class_number", "degree", "disc_sign", "r2"]: # integer columns
                        maxval_update(flds, maxvals, fld, rec)
                    flds["reg"] = realwrap(rec.get("reg"))
                    for fld in ["class_group", "coeffs"]:
                        flds[fld] = splitwrap(rec.get(fld))
                    flds["ramps"] = ijsonwrap(rec.get("ramps"))
                    if "ramps" in rec:
                        flds["disc_rad"] = str(prod(int(p) for p in rec["ramps"]))
                    else:
                        flds["disc_rad"] = r'\N'
                    #for fld in ["subs", "zk", "unitsGmodule", "units", "res", "loc_algebras"]: # json columns
                    #    flds[fld] = jsonwrap(rec.get(fld))
                    if "galois" in rec:
                        flds["galt"] = intwrap(rec["galois"]["t"])
                    else:
                        flds["galt"] = r'\N'
                    for fld in ["used_grh"]: # bool columns
                        flds[fld] = boolwrap(rec.get(fld))
                    for fld in ["label", "unitsType"]:
                        flds[fld] = strwrap(rec.get(fld))
                    # iso_number
                    flds["iso_number"] = intwrap(rec["label"].split(".")[-1])
                    # cm
                    flds["subs"] = jsonwrap(rec.get("subs"))
                    subs = flds["subs"]
                    if rec["r1"] != '0':
                        cm = 'f'
                    elif rec["r2"] == '1':
                        cm = 't'
                    elif subs == r'\N':
                        cm = r'\N'
                    elif len(subs) == 2: # []
                        cm = 'f'
                    else:
                        subfields = [s[1:s.rfind('"')] for s in subs[2:-2].split('],[')]
                        for subfield in subfields:
                            if subfield.count(',') == int(flds["r2"]):
                                subfield = pari("[" + subfield + "]").Polrev("x")
                                nr_real_rts = subfield.polsturm()
                                if nr_real_rts == int(flds["r2"]):
                                    cm = 't'
                                    break
                        else:
                            cm = 'f'
                    flds["cm"] = cm
                    Fsmall.write('\t'.join(flds[fld] for fld in ["id"] + search_cols) + '\n')
                    #Fextras.write('\t'.join(flds[fld] for fld in ["id"] + extra_cols) + '\n')
        search_types = defaultdict(list)
        extra_types = defaultdict(list)
        for fld in ["class_number", "degree", "disc_abs", "disc_sign", "r2"]:
            search_types[integer_size(maxvals[fld])].append(fld)
        search_types["smallint"].append("iso_number")
        search_types["integer"].append("galt")
        search_types["numeric"].extend(["disc_abs", "disc_rad"])
        extra_types["numeric"].append("reg")
        search_types["jsonb"].extend(["class_group", "coeffs", "ramps"])
        extra_types["jsonb"].extend(["subs", "zk", "unitsGmodule", "units", "loc_algebras", "res"])
        search_types["boolean"].extend(["cm", "used_grh"])
        search_types["text"].append("label")
        extra_types["text"].append("unitsType")
        with open('import_special.py', 'a') as Fimp:
            Fimp.write('''
def import_nf_fields():
    try:
        db.create_table('nf_fields', {0}, 'label', ['degree', 'disc_abs', 'disc_sign', 'iso_number'], True, {1}, search_order={2}, extra_order={3})
    except Exception:
        print "Failure in creating nf_fields"
        traceback.print_exc()
        return
    try:
        db.nf_fields.copy_from('/scratch/importing/nf_fields.txt', '/scratch/importing/nf_fields_extras.txt', search_cols={2}, extra_cols={3}, includes_ids=True, resort=False)
    except Exception:
        print "Failure loading data into nf_fields"
        traceback.print_exc()
        return
    print "Successfully imported nf_fields"
    index_nf_fields()
def index_nf_fields():
    try:
        print "Indexing nf_fields"
        db.nf_fields.restore_pkeys()
        db.nf_fields.create_index(['ramps'], 'gin')
        db.nf_fields.create_index(['label'])
        db.nf_fields.create_index(['class_group', 'id'])
        db.nf_fields.create_index(['class_number', 'id'])
        db.nf_fields.create_index(['coeffs'])
        db.nf_fields.create_index(['degree', 'disc_abs', 'disc_sign', 'iso_number', 'r2'])
        db.nf_fields.create_index(['disc_sign', 'disc_abs', 'id'])
        db.nf_fields.create_index(['disc_rad', 'id'])
        db.nf_fields.create_index(['degree', 'galt', 'disc_abs', 'disc_sign'])
        db.nf_fields.create_index(['degree', 'galt', 'id'])
        db.nf_fields.create_index(['degree', 'r2', 'id'])
    except Exception:
        print "Failure in indexing nf_fields"
        traceback.print_exc()
    else:
        print "Successfully indexed nf_fields"
'''.format(dict(search_types), dict(extra_types), search_cols, extra_cols))
    except Exception:
        print "Failure in exporting nf_fields"
        traceback.print_exc()
        with open("export_failures.txt",'a') as F:
            F.write("nf_fields\n")
    else:
        print "Successfully exported nf_fields"

def export_lfunc_zeros():
    conn = lmfdb.base.getDBConnection()
    coll = conn.test.Lfunctions_test2
    ordered_cols = ['degree', 'level', 'signature', 'first_zero', 'description', 'coeffs', 'eta', 'mu', 'special']
    icols = ['degree', 'level']
    jcols = ['signature', 'coeffs', 'eta', 'mu', 'special']
    rcols = ['first_zero']
    tcols = ['description']
    try:
        with open('exports/lfunc_zeros.txt', 'w') as Fout:
            for rec in coll.find():
                flds = {}
                for fld in icols:
                    flds[fld] = intwrap(rec.get(fld))
                for fld in jcols:
                    flds[fld] = jsonwrap(rec.get(fld))
                for fld in rcols:
                    flds[fld] = realwrap(rec.get(fld))
                for fld in tcols:
                    flds[fld] = strwrap(rec.get(fld))
                Fout.write("\t".join(flds[fld] for fld in ordered_cols) + "\n")
        types = defaultdict(list)
        types['smallint'] = icols
        types['jsonb'] = jcols
        types['numeric'] = rcols
        types['text'] = tcols
        with open('/home/roed/import_special.py', 'a') as Fimp:
            Fimp.write('''
def import_lfunc_zeros():
    try:
        db.create_table('lfunc_zeros', {0}, None, ['first_zero'], search_order={1})
        db.lfunc_zeros.copy_from('/scratch/importing/lfunc_zeros.txt', search_cols={1}, includes_ids=True, resort=False)
    except Exception:
        print "Failure in importing lfunc_zeros"
        traceback.print_exc()
    else:
        print "Successfully imported lfunc_zeros"
'''.format(dict(types), ordered_cols))
    except Exception:
        print "Failure in exporting lfunc_zeros"
        traceback.print_exc()
        with open("export_failures.txt",'a') as F:
            F.write("lfunc_zeros\n")
    else:
        print "Successfully exported lfunc_zeros"

def export_mwfp_forms():
    conn = lmfdb.base.getDBConnection()
    mwfp_forms = conn.HTPicard.picard
    maxvals = defaultdict(int)
    try:
        ordered_cols = ["maass_id", "ev", "prec", "sym", "coef"]
        with open('exports/mwfp_forms.txt', 'w') as Fout:
            for i, rec in sort_collection(mwfp_forms, ['ev'], 'HTPicard.picard'):
                flds = {}
                flds["maass_id"] = str(rec["_id"])
                for fld in ["ev", "prec"]: # real columns
                    flds[fld] = realwrap(rec.get(fld))
                for fld in ["sym", "coef"]: # json columns
                    flds[fld] = jsonwrap(rec.get(fld))
                flds["id"] = str(i)
                Fout.write("\t".join(flds[fld] for fld in ["id"] + ordered_cols) + "\n")
        types = defaultdict(list)

        for fld in ["ev", "prec"]:
            types['numeric'].append(fld)
        for fld in ["sym", "coef"]:
            types['jsonb'].append(fld)
        types['text'] = ["maass_id"]
        with open('import.py', 'a') as Fimp:
            Fimp.write('''
def import_mwfp_forms():
    try:
        db.create_table('mwfp_forms', {0}, 'maass_id', ['ev'], False, search_order={1})
    except Exception:
        print "Failure in creating mwfp_forms"
        traceback.print_exc()
        return
    try:
        db.mwfp_forms.copy_from('/scratch/importing/mwfp_forms.txt', search_cols={1}, includes_ids=True, resort=False)
    except Exception:
        print "Failure loading data into mwfp_forms"
        traceback.print_exc()
        return
    print "Successfully imported mwfp_forms"
    index_mwfp_forms()
def index_mwfp_forms():
    try:
        print "Indexing mwfp_forms"
        db.mwfp_forms.restore_pkeys()
        db.mwfp_forms.create_index([u'ev'], 'btree')
    except Exception:
        print "Failure in indexing mwfp_forms"
        traceback.print_exc()
    else:
        print "Successfully indexed mwfp_forms"
'''.format(dict(types), ordered_cols))
    except Exception:
        print "Failure in exporting mwfp_forms"
        traceback.print_exc()
        with open("export_failures.txt",'a') as F:
            F.write("mwfp_forms\n")
    else:
        print "Successfully exported mwfp_forms"

def export_mwf_coeffs():
    conn = lmfdb.base.getDBConnection()
    mwf_coeffs = conn.MaassWaveForms.Coefficients
    mwf_gridfs = gridfs.GridFS(conn.MaassWaveForms, 'Coefficients')
    maxvals = defaultdict(int)
    Numc = {}
    try:
        with open('exports/mwf_coeffs.txt', 'w') as Fout:
            for i, rec in sort_collection(mwf_coeffs, None, 'MaassWaveForms.Coefficients'):
                cid = rec.get('coeff_id')
                if cid is not None and mwf_gridfs.exists(cid):
                    gfile = mwf_gridfs.get(cid)
                    Numc[gfile.filename] = rec.get('Numc')
            for filename in mwf_gridfs.list():
                flds = {}
                gfile = mwf_gridfs.find_one({'filename':filename})
                flds['coefficients'] = binarywrap(gfile.read())
                flds['label'] = filename
                flds['Numc'] = intwrap(Numc.get(filename))
                Fout.write("\t".join(flds[fld] for fld in ["label", "Numc", "coefficients"]) + "\n")
        with open('/home/roed/import_special.py', 'a') as Fimp:
            Fimp.write('''
def import_mwf_coeffs():
    try:
        db.create_table('mwf_coeffs', {'int':'Numc', 'bytea':'coefficients', 'text':'label'}, 'label', search_order=['label', 'Numc', 'coefficients'])
        db.mwf_coeffs.copy_from('/scratch/importing/mwf_coeffs.txt', search_cols=['label', 'Numc', 'coefficients'], includes_ids=True, resort=False)
    except DatabaseError as err:
        print "Failure in importing mwf_coeffs"
        traceback.print_exc()
    else:
        print "Successfully imported mwf_coeffs"
''')
    except Exception as err:
        print "Failure in exporting mwf_coeffs"
        traceback.print_exc()
        with open("export_failures.txt",'a') as F:
            F.write("mwf_coeffs\n")
    else:
        print "Successfully exported mwf_coeffs"

def export_mwf_forms():
    conn = lmfdb.base.getDBConnection()
    mwf_forms = conn.MaassWaveForms.FS
    mwf_gridfs = gridfs.GridFS(conn.MaassWaveForms, 'Coefficients')
    maxvals = defaultdict(int)
    ordered_cols = ["maass_id", "Eigenvalue", "Character", "date", "M0", "Norm", "Coefficient", "ObjectUrl", "Precision", "Numc", "Cusp_evs", "Dim", "Level", "Newform", "Error", "Y", "Dimension", "dim", "Symmetry", "Weight", "Fricke", "ObjectLabel", "Sign", "coeff_label", "Contributor", "software"]
    try:
        with open('exports/mwf_forms.txt', 'w') as Fout:
            for i, rec in sort_collection(mwf_forms, ['Weight', 'Level', 'Character', 'Eigenvalue'], 'MaassWaveForms.FS'):
                flds = {}

                for fld in ["Eigenvalue", "Precision", "Error", "Y"]: # real columns
                    flds[fld] = realwrap(rec.get(fld))
                if rec.get("Conrey") != 1:
                    N = rec['Level']
                    chi = rec['Character']
                    new_chi = MaassDB().getDircharConrey(N, chi)
                    rec['Character'] = new_chi
                for fld in ["Character", "M0", "Norm", "Numc", "Dim", "Level", "Dimension", "dim", "Symmetry", "Weight", "Fricke", "Sign"]: # integer columns
                    maxval_update(flds, maxvals, fld, rec)
                for fld in ["date"]: # datetime columns
                    flds[fld] = datetimewrap(rec.get(fld))
                for fld in ["Coefficient", "Cusp_evs"]: # json columns
                    flds[fld] = jsonwrap(rec.get(fld))
                for fld in ["Newform"]: # bool columns
                    flds[fld] = boolwrap(rec.get(fld))
                for fld in ["ObjectUrl", "ObjectLabel", "Contributor", "software"]: # str columns
                    flds[fld] = strwrap(rec.get(fld))
                flds['maass_id'] = str(rec['_id'])
                cid = rec.get('coeff_id')
                if cid is None or not mwf_gridfs.exists(cid):
                    flds['coeff_label'] = r'\N'
                else:
                    flds['coeff_label'] = mwf_gridfs.get(cid).filename
                flds["id"] = str(i)
                Fout.write("\t".join(flds[fld] for fld in ["id"] + ordered_cols) + "\n")
        types = defaultdict(list)
        types['numeric'] = ["Eigenvalue", "Precision", "Error", "Y"]
        types['timestamp'] = ["date"]
        types['jsonb'] = ["Coefficient", "Cusp_evs"]
        types['boolean'] = ["Newform"]
        types['text'] = ["ObjectUrl", "ObjectLabel", "Contributor", "software", "coeff_label", "maass_id"]
        for fld in ["Character", "M0", "Norm", "Numc", "Dim", "Level", "Dimension", "dim", "Symmetry", "Weight", "Fricke", "Sign"]:
            types[integer_size(maxvals[fld])].append(fld)
        with open('/home/roed/import_special.py', 'a') as Fimp:
            Fimp.write('''
def import_mwf_forms():
    try:
        db.create_table('mwf_forms', {0}, 'maass_id', sort=['Weight', 'Level', 'Character', 'Eigenvalue'], id_ordered=True, search_order={1})
        db.mwf_forms.copy_from('/scratch/importing/mwf_forms.txt', search_cols={1}, includes_ids=True, resort=False)
        db.mwf_forms._set_ordered()
        db.mwf_forms.restore_pkeys()
        db.mwf_forms.create_index(["Character"])
        db.mwf_forms.create_index(["Eigenvalue", "id"])
        db.mwf_forms.create_index(["Level", "Weight", "id"])
        #db.mwf_forms.create_index(["Coefficient"])
        #db.mwf_forms.create_index(["Level"])
        #db.mwf_forms.create_index(["Numc"])
        #db.mwf_forms.create_index(["Weight"])
    except DatabaseError:
        print "Failure in importing mwf_forms"
        traceback.print_exc()
    else:
        print "Successfully imported mwf_forms"
'''.format(dict(types), ordered_cols))
    except Exception as err:
        #print "Failure in exporting mwf_forms"
        #traceback.print_exc()
        raise
    else:
        print "Successfully exported mwf_forms"

def export_mwf_plots():
    conn = lmfdb.base.getDBConnection()
    mwf_plots = conn.MaassWaveForms.maassform_plots
    maxvals = defaultdict(int)
    ordered_cols = ["eigenvalue", "plot", "level", "maass_id", "num_pts", "dpi"]
    try:
        with open('exports/mwf_plots.txt', 'w') as Fout:
            for i, rec in sort_collection(mwf_plots, None, 'MaassWaveForms.maassform_plots'):
                flds = {}
                flds["eigenvalue"] = realwrap(rec.get("eigenvalue"))
                flds["plot"] = binarywrap(rec.get("plot"))
                for fld in ["level", "num_pts", "dpi"]: # integer columns
                    maxval_update(flds, maxvals, fld, rec)
                flds["maass_id"] = str(rec["maass_id"]) if rec.get("maass_id") else r"\N"
                Fout.write("\t".join(flds[fld] for fld in ordered_cols) + "\n")
        types = defaultdict(list)
        types['numeric'] = ["eigenvalue"]
        types['text'] = ["maass_id"]
        types['bytea'] = ["plot"]
        for fld in ["level", "num_pts", "dpi"]:
            types[integer_size(maxvals[fld])].append(fld)
        with open('/home/roed/import_special.py', 'a') as Fimp:
            Fimp.write('''
def import_mwf_plots():
    try:
        db.create_table('mwf_plots', {0}, 'maass_id', search_order={1})
        db.mwf_plots.copy_from('/scratch/importing/mwf_plots.txt', search_cols={1}, includes_ids=True, resort=False)
        db.mwf_plots.create_index(["maass_id"])
    except DatabaseError as err:
        print "Failure in importing mwf_plots"
        print err
    else:
        print "Successfully imported mwf_plots"
'''.format(dict(types), ordered_cols))
    except Exception as err:
        #print "Failure in exporting mwf_plots"
        #traceback.print_exc()
        raise
        with open("export_failures.txt",'a') as F:
            F.write("mwf_plots\n")
    else:
        print "Successfully exported mwf_plots"

def export_mwf_tables():
    conn = lmfdb.base.getDBConnection()
    fs = gridfs.GridFS(conn.MaassWaveForms, 'Table')
    ff = fs.get_version(filename=fs.list()[0])
    table = loads(ff.read())
    ordered_cols = ['ncols', 'nrows', 'levels', 'characters', 'weights', 'keylist', 'data']
    icols = ['ncols', 'nrows']
    jcols = ['levels', 'characters', 'weights', 'keylist']
    with open('exports/mwf_tables.txt', 'w') as Fout:
        flds = {}
        for fld in icols:
            flds[fld] = intwrap(table[fld])
        for fld in jcols:
            flds[fld] = jsonwrap(table[fld])
        flds['data'] = jsonwrap({','.join(str(x) for x in k):v for k,v in table['data'].iteritems()})
        Fout.write("\t".join(flds[fld] for fld in ordered_cols) + "\n")
    types = {'smallint':icols, 'jsonb':jcols+['data']}
    with open('/home/roed/import_special.py', 'a') as Fimp:
        Fimp.write('''
def import_mwf_tables():
    db.create_table('mwf_tables', {0}, None, search_order={1})
    db.mwf_tables.copy_from('/scratch/importing/mwf_tables.txt', search_cols={1}, includes_ids=True, resort=False)
    print "Successfully imported mwf_tables"
'''.format(types, ordered_cols))

maxid = 0
def export_smf_samples():
    conn = lmfdb.base.getDBConnection()
    smf_samples = conn.siegel_modular_forms.samples
    maxvals = defaultdict(int)
    icols = ["weight", "degree", "fdeg", "representation"]
    jcols = ["Fourier_coefficients", "collection", "eigenvalues"]
    bcols = ["is_integral", "is_eigenform"]
    tcols = ["field_poly", "field", "type", "name", "explicit_formula", "courtesy_of"]

    sordered_cols = ["id_link", "name", "weight", "degree", "field_poly", "fdeg", "is_integral", "Fourier_coefficients", "field", "type", "collection", "courtesy_of", "eigenvalues", "is_eigenform", "representation", "explicit_formula"]
    fordered_cols = ["owner_id", "det", "data"]
    eordered_cols = ["owner_id", "index", "data"]
    id_links = {}
    def get_id_link(_id):
        global maxid
        _id = str(_id)
        if _id in id_links:
            return id_links[_id]
        else:
            maxid += 1
            id_links[_id] = maxid
            return maxid
    try:
        with open('exports/smf_samples.txt', 'w') as Fsamp:
            with open('exports/smf_ev.txt', 'w') as Fev:
                with open('exports/smf_fc.txt', 'w') as Ffc:
                    for i, rec in sort_collection(smf_samples, None, 'siegel_modular_forms.samples'):
                        flds = {}
                        if rec['data_type'] == 'sample':
                            ordered_cols = sordered_cols
                            flds["id_link"] = intwrap(get_id_link(rec["_id"]))
                            for fld in icols:
                                flds[fld] = intwrap(rec.get(fld))
                            for fld in jcols:
                                flds[fld] = jsonwrap(rec.get(fld))
                            for fld in bcols:
                                flds[fld] = boolwrap(rec.get(fld))
                            for fld in tcols:
                                flds[fld] = strwrap(rec.get(fld))
                            F = Fsamp
                        elif rec['data_type'] == 'fc':
                            ordered_cols = fordered_cols
                            flds["det"] = intwrap(rec.get("det"))
                            flds["data"] = jsonwrap(rec.get("data"))
                            flds["owner_id"] = intwrap(get_id_link(rec["owner_id"]))
                            F = Ffc
                        elif rec['data_type'] == 'ev':
                            ordered_cols = eordered_cols
                            flds["index"] = intwrap(rec.get("index"))
                            flds["data"] = jsonwrap(rec.get("data"))
                            flds["owner_id"] = intwrap(get_id_link(rec["owner_id"]))
                            F = Fev
                        else:
                            raise RuntimeError
                        F.write("\t".join(flds[fld] for fld in ordered_cols) + "\n")
        stypes = {'smallint':icols, 'integer':['id_link'], 'jsonb':jcols, 'boolean':bcols, 'text':tcols}
        ftypes = {'smallint':'det', 'jsonb':'data', 'integer':'owner_id'}
        etypes = {'smallint':'index', 'jsonb':'data', 'integer':'owner_id'}
        with open('/home/roed/import_special.py', 'a') as Fimp:
            Fimp.write('''
def import_smf_samples():
    try:
        db.create_table('smf_samples', {0}, None, ['name'], search_order={1})
        db.create_table('smf_ev', {2}, None, ['index'], search_order={3})
        db.create_table('smf_fc', {4}, None, ['det'], search_order={5})
        db.smf_samples.copy_from('/scratch/importing/smf_samples.txt', search_cols={1}, includes_ids=True, resort=False)
        db.smf_ev.copy_from('/scratch/importing/smf_ev.txt', search_cols={3}, includes_ids=True, resort=False)
        db.smf_fc.copy_from('/scratch/importing/smf_fc.txt', search_cols={5}, includes_ids=True, resort=False)
        db.smf_samples.create_index(['collection'], 'gin')
        db.smf_samples.create_index(['weight'])
        db.smf_ev.create_index(['owner_id', 'index'])
        db.smf_fc.create_index(['owner_id', 'det'])
    except DatabaseError as err:
        print "Failure in importing smf_samples"
        print err
    else:
        print "Successfully imported smf_samples"
'''.format(stypes, sordered_cols, etypes, eordered_cols, ftypes, fordered_cols))
    except Exception as err:
        #print "Failure in exporting smf_samples"
        #traceback.print_exc()
        raise
        with open("export_failures.txt",'a') as F:
            F.write("smf_samples\n")
    else:
        print "Successfully exported smf_samples"

def export_artin_reps():
    conn = lmfdb.base.getDBConnection()
    artin_reps = conn.artin.representations
    maxvals = defaultdict(int)
    try:
        ordered_cols = ["Baselabel", "Dim", "Conductor", "Galn", "Galt", "Indicator", "BadPrimes", "HardPrimes", "GaloisConjugates", "GalConjSigns", "CharacterField", "NFGal", "Hide"]
        with open('exports/artin_reps.txt', 'w') as Fout:
            for i, rec in sort_collection(artin_reps, ['Dim', 'Conductor_key'], 'artin.representations'):
                flds = {}

                for fld in ["Dim", "Indicator", "CharacterField", "Hide"]: # integer columns
                    maxval_update(flds, maxvals, fld, rec)
                for fld in ["Conductor"]: # int_as_str columns
                    flds[fld] = intwrap(rec.get(fld))
                for fld in ["NFGal"]: # list_as_str columns
                    flds[fld] = splitwrap(rec.get(fld))
                if "GaloisConjugates" in rec:
                    galconj = rec["GaloisConjugates"]
                    flds["GaloisConjugates"] = jsonwrap(galconj)
                    flds["GalConjSigns"] = jsonwrap([X['Sign'] for X in galconj])
                else:
                    flds["GaloisConjugates"] = r'\N'
                    flds["GalConjSigns"] = r'\N'
                for fld in ["BadPrimes", "HardPrimes"]: # json columns
                    flds[fld] = ijsonwrap(rec.get(fld))
                galnt = rec.get("Galois_nt")
                if galnt is None:
                    flds["Galt"] = r'\N'
                    flds["Galn"] = r'\N'
                else:
                    flds["Galn"] = intwrap(galnt[0])
                    flds["Galt"] = intwrap(galnt[1])
                for fld in ["Baselabel"]: # str columns
                    flds[fld] = strwrap(rec.get(fld))
                if i is not None:
                    flds["id"] = str(i)
                Fout.write("\t".join(flds[fld] for fld in ["id"] + ordered_cols) + "\n")
        types = defaultdict(list)
        for fld in ["Dim", "Indicator", "CharacterField", "Hide"]:
            types[integer_size(maxvals[fld])].append(fld)
        for fld in ["Conductor"]:
            types['numeric'].append(fld)
        for fld in ["NFGal"]:
            types['jsonb'].append(fld)
        for fld in ["BadPrimes", "GaloisConjugates", "HardPrimes", "GalConjSigns"]:
            types['jsonb'].append(fld)
        for fld in ["Baselabel"]:
            types['text'].append(fld)
        types['smallint'].append("Galn")
        types['integer'].append("Galt")
        with open('import_special.py', 'a') as Fimp:
            Fimp.write('''
def import_artin_reps():
    try:
        db.create_table('artin_reps', {0}, 'Baselabel', ['Dim', 'Conductor'], True, search_order={1})
    except Exception:
        print "Failure in creating artin_reps"
        traceback.print_exc()
        return
    try:
        db.artin_reps.copy_from('/scratch/importing/artin_reps.txt', search_cols={1}, includes_ids=True, resort=False)
    except Exception:
        print "Failure loading data into artin_reps"
        traceback.print_exc()
        return
    print "Successfully imported artin_reps"
    index_artin_reps()
def index_artin_reps():
    try:
        print "Indexing artin_reps"
        db.artin_reps.restore_pkeys()
        db.artin_reps.create_index([u'BadPrimes'], 'gin')
        db.artin_reps.create_index([u'Baselabel'], 'btree')
        db.artin_reps.create_index(['Dim', 'Hide', 'Conductor'], 'btree')
        db.artin_reps.create_index(['Hide', 'Conductor'], 'btree')
        db.artin_reps.create_index([u'Hide', u'Galn', u'Galt'], 'btree')
        db.artin_reps.create_index([u'NFGal'], 'btree')
    except Exception:
        print "Failure in indexing artin_reps"
        traceback.print_exc()
    else:
        print "Successfully indexed artin_reps"
'''.format(dict(types), ordered_cols))
    except Exception:
        print "Failure in exporting artin_reps"
        traceback.print_exc()
        with open("export_failures.txt",'a') as F:
            F.write("artin_reps\n")
    else:
        print "Successfully exported artin_reps"

def export_artin_field_data():
    conn = lmfdb.base.getDBConnection()
    artin_field_data = conn.artin.field_data
    maxvals = defaultdict(int)
    try:
        ordered_cols = ["ConjClasses", "G-Name", "QpRts-minpoly", "ComplexConjugation", "FrobResolvents", "QpRts-prec", "ArtinReps", "G-Gens", "Polynomial", "QpRts", "Frobs", "QpRts-p", "TransitiveDegree", "Size"]
        with open('exports/artin_field_data.txt', 'w') as Fout:
            for i, rec in sort_collection(artin_field_data, None, 'artin.field_data'):
                flds = {}

                for fld in ["ComplexConjugation", "QpRts-prec", "QpRts-p", "TransitiveDegree"]: # integer columns
                    maxval_update(flds, maxvals, fld, rec)
                for fld in ["Size"]: # int_as_str columns
                    flds[fld] = intwrap(rec.get(fld))
                if "Polynomial" in rec: # list_as_str columns
                    flds["Polynomial"] = splitwrap(rec.get("Polynomial"))
                else:
                    flds["Polynomial"] = splitwrap(rec.get("NFGal"))
                for fld in ["ConjClasses", "QpRts-minpoly", "FrobResolvents", "ArtinReps", "G-Gens", "QpRts", "Frobs"]: # json columns
                    flds[fld] = jsonwrap(rec.get(fld))
                for fld in ["G-Name"]: # str columns
                    flds[fld] = strwrap(rec.get(fld))
                if i is not None:
                    flds["id"] = str(i)
                Fout.write("\t".join(flds[fld] for fld in ordered_cols) + "\n")
        types = defaultdict(list)
        for fld in ["ComplexConjugation", "QpRts-prec", "QpRts-p", "TransitiveDegree"]:
            types[integer_size(maxvals[fld])].append(fld)
        for fld in ["Size"]:
            types['numeric'].append(fld)
        for fld in ["Polynomial"]:
            types['jsonb'].append(fld)
        for fld in ["ConjClasses", "QpRts-minpoly", "FrobResolvents", "ArtinReps", "G-Gens", "QpRts", "Frobs"]:
            types['jsonb'].append(fld)
        for fld in ["G-Name"]:
            types['text'].append(fld)
        with open('import_special.py', 'a') as Fimp:
            Fimp.write('''
def import_artin_field_data():
    try:
        db.create_table('artin_field_data', {0}, None, None, False, search_order={1})
    except Exception:
        print "Failure in creating artin_field_data"
        traceback.print_exc()
        return
    try:
            db.artin_field_data.copy_from('/scratch/importing/artin_field_data.txt', search_cols={1}, includes_ids=True, resort=False)
    except Exception:
        print "Failure loading data into artin_field_data"
        traceback.print_exc()
        return
    print "Successfully imported artin_field_data"
    index_artin_field_data()
def index_artin_field_data():
    try:
        print "Indexing artin_field_data"
        db.artin_field_data.create_index([u'Polynomial'], 'btree')
    except Exception:
        print "Failure in indexing artin_field_data"
        traceback.print_exc()
    else:
        print "Successfully indexed artin_field_data"
'''.format(dict(types), ordered_cols))
    except Exception:
        print "Failure in exporting artin_field_data"
        traceback.print_exc()
        with open("export_failures.txt",'a') as F:
            F.write("artin_field_data\n")
    else:
        print "Successfully exported artin_field_data"

def export_oldstats(tables):
    with open('/home/roed/import_special.py', 'a') as Fimp:
        Fimp.write("def import_oldstats():\n")
        for collection, search_table in tables:
            filename = "exports/" + search_table + "_oldstats.txt"
            with open(filename, 'w') as F:
                for D in collection.find():
                    if search_table == 'av_fqisog':
                        Did = D.pop('label')
                        D = D['counts']
                    else:
                        Did = D.pop('_id')
                    F.write(Did + "\t" + jsonwrap(D) + "\n")
            Fimp.write("    db.{0}.stats.create_oldstats('/scratch/importing/{0}_oldstats.txt')\n".format(search_table))

def export_special():
    with open('import_special.py', 'w') as Fimp:
        Fimp.write('''
from psycopg2 import DatabaseError, connect
from lmfdb.db_backend import db
import traceback
def import_special():
    import_mwfp_forms()
    import_g2c_instances()
    import_users()
    import_lfunc_zeros
    import_mwf_coeffs()
    import_mwf_forms()
    import_mwf_plots()
    import_mwf_tables()
    import_smf_samples()
    import_artin_reps()
    import_artin_field_data()
    import_lfunc_zeros()
    import_oldstats()
''')
    export_mwfp_forms()
    export_g2c_instances()
    export_users()
    export_lfunc_zeros
    export_mwf_coeffs()
    export_mwf_forms()
    export_mwf_plots()
    export_mwf_tables()
    export_smf_samples()
    export_artin_reps()
    export_artin_field_data()
    export_lfunc_zeros()
    conn = lmfdb.base.getDBConnection()
    export_oldstats([(conn.elliptic_curves.nfcurves.stats, 'ec_nfcurves'),
                     (conn.Lattices.lat.stats, 'lat_lattices'),
                     (conn.elliptic_curves.curves.stats, 'ec_curves'),
                     (conn.hmfs.forms.stats, 'hmf_forms'),
                     (conn.numberfields.stats, 'nf_fields'),
                     (conn.hecke_algebras.hecke_algebras.stats, 'hecke_algebras'),
                     (conn.abvar.fq_isog.stats, 'av_fqisog'),
                     (conn.curve_automorphisms.passports.stats, 'hgcwa_passports')])
