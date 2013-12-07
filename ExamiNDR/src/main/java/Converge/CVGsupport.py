# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.
import _CVGsupport
def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class IntVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, IntVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, IntVector, name)
    def __init__(self,*args):
        _swig_setattr(self, IntVector, 'this', apply(_CVGsupport.new_IntVector,args))
        _swig_setattr(self, IntVector, 'thisown', 1)
    def __len__(*args): return apply(_CVGsupport.IntVector___len__,args)
    def __nonzero__(*args): return apply(_CVGsupport.IntVector___nonzero__,args)
    def clear(*args): return apply(_CVGsupport.IntVector_clear,args)
    def append(*args): return apply(_CVGsupport.IntVector_append,args)
    def pop(*args): return apply(_CVGsupport.IntVector_pop,args)
    def __getitem__(*args): return apply(_CVGsupport.IntVector___getitem__,args)
    def __getslice__(*args): return apply(_CVGsupport.IntVector___getslice__,args)
    def __setitem__(*args): return apply(_CVGsupport.IntVector___setitem__,args)
    def __setslice__(*args): return apply(_CVGsupport.IntVector___setslice__,args)
    def __delitem__(*args): return apply(_CVGsupport.IntVector___delitem__,args)
    def __delslice__(*args): return apply(_CVGsupport.IntVector___delslice__,args)
    def __del__(self, destroy= _CVGsupport.delete_IntVector):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C IntVector instance at %s>" % (self.this,)

class IntVectorPtr(IntVector):
    def __init__(self,this):
        _swig_setattr(self, IntVector, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, IntVector, 'thisown', 0)
        _swig_setattr(self, IntVector,self.__class__,IntVector)
_CVGsupport.IntVector_swigregister(IntVectorPtr)

class DoubleVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DoubleVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DoubleVector, name)
    def __init__(self,*args):
        _swig_setattr(self, DoubleVector, 'this', apply(_CVGsupport.new_DoubleVector,args))
        _swig_setattr(self, DoubleVector, 'thisown', 1)
    def __len__(*args): return apply(_CVGsupport.DoubleVector___len__,args)
    def __nonzero__(*args): return apply(_CVGsupport.DoubleVector___nonzero__,args)
    def clear(*args): return apply(_CVGsupport.DoubleVector_clear,args)
    def append(*args): return apply(_CVGsupport.DoubleVector_append,args)
    def pop(*args): return apply(_CVGsupport.DoubleVector_pop,args)
    def __getitem__(*args): return apply(_CVGsupport.DoubleVector___getitem__,args)
    def __getslice__(*args): return apply(_CVGsupport.DoubleVector___getslice__,args)
    def __setitem__(*args): return apply(_CVGsupport.DoubleVector___setitem__,args)
    def __setslice__(*args): return apply(_CVGsupport.DoubleVector___setslice__,args)
    def __delitem__(*args): return apply(_CVGsupport.DoubleVector___delitem__,args)
    def __delslice__(*args): return apply(_CVGsupport.DoubleVector___delslice__,args)
    def __del__(self, destroy= _CVGsupport.delete_DoubleVector):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C DoubleVector instance at %s>" % (self.this,)

class DoubleVectorPtr(DoubleVector):
    def __init__(self,this):
        _swig_setattr(self, DoubleVector, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, DoubleVector, 'thisown', 0)
        _swig_setattr(self, DoubleVector,self.__class__,DoubleVector)
_CVGsupport.DoubleVector_swigregister(DoubleVectorPtr)

def Motif2c_PSSM(motif):
    c_PSSM = SeqMat(motif.width)
    LjT = zip(['A','C','G','T'],[0,1,2,3])
    for i in range(motif.width):
        for L,j in LjT:
            c_PSSM.set(i,j, motif.ll[i][L])
    return c_PSSM


seq2int = _CVGsupport.seq2int

print_seq = _CVGsupport.print_seq

list2double = _CVGsupport.list2double

printdouble = _CVGsupport.printdouble

log2_sum = _CVGsupport.log2_sum

rc = _CVGsupport.rc

class Probe(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Probe, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Probe, name)
    __swig_setmethods__["iseq"] = _CVGsupport.Probe_iseq_set
    __swig_getmethods__["iseq"] = _CVGsupport.Probe_iseq_get
    if _newclass:iseq = property(_CVGsupport.Probe_iseq_get,_CVGsupport.Probe_iseq_set)
    __swig_setmethods__["gaps"] = _CVGsupport.Probe_gaps_set
    __swig_getmethods__["gaps"] = _CVGsupport.Probe_gaps_get
    if _newclass:gaps = property(_CVGsupport.Probe_gaps_get,_CVGsupport.Probe_gaps_set)
    __swig_setmethods__["len"] = _CVGsupport.Probe_len_set
    __swig_getmethods__["len"] = _CVGsupport.Probe_len_get
    if _newclass:len = property(_CVGsupport.Probe_len_get,_CVGsupport.Probe_len_set)
    __swig_setmethods__["genomes"] = _CVGsupport.Probe_genomes_set
    __swig_getmethods__["genomes"] = _CVGsupport.Probe_genomes_get
    if _newclass:genomes = property(_CVGsupport.Probe_genomes_get,_CVGsupport.Probe_genomes_set)
    __swig_setmethods__["probebg"] = _CVGsupport.Probe_probebg_set
    __swig_getmethods__["probebg"] = _CVGsupport.Probe_probebg_get
    if _newclass:probebg = property(_CVGsupport.Probe_probebg_get,_CVGsupport.Probe_probebg_set)
    __swig_setmethods__["seqbgs"] = _CVGsupport.Probe_seqbgs_set
    __swig_getmethods__["seqbgs"] = _CVGsupport.Probe_seqbgs_get
    if _newclass:seqbgs = property(_CVGsupport.Probe_seqbgs_get,_CVGsupport.Probe_seqbgs_set)
    __swig_setmethods__["wmerbgs"] = _CVGsupport.Probe_wmerbgs_set
    __swig_getmethods__["wmerbgs"] = _CVGsupport.Probe_wmerbgs_get
    if _newclass:wmerbgs = property(_CVGsupport.Probe_wmerbgs_get,_CVGsupport.Probe_wmerbgs_set)
    __swig_setmethods__["Zs"] = _CVGsupport.Probe_Zs_set
    __swig_getmethods__["Zs"] = _CVGsupport.Probe_Zs_get
    if _newclass:Zs = property(_CVGsupport.Probe_Zs_get,_CVGsupport.Probe_Zs_set)
    __swig_setmethods__["Hs"] = _CVGsupport.Probe_Hs_set
    __swig_getmethods__["Hs"] = _CVGsupport.Probe_Hs_get
    if _newclass:Hs = property(_CVGsupport.Probe_Hs_get,_CVGsupport.Probe_Hs_set)
    def __init__(self,*args):
        _swig_setattr(self, Probe, 'this', apply(_CVGsupport.new_Probe,args))
        _swig_setattr(self, Probe, 'thisown', 1)
    def calc_gaps(*args): return apply(_CVGsupport.Probe_calc_gaps,args)
    def __del__(self, destroy= _CVGsupport.delete_Probe):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C Probe instance at %s>" % (self.this,)

class ProbePtr(Probe):
    def __init__(self,this):
        _swig_setattr(self, Probe, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Probe, 'thisown', 0)
        _swig_setattr(self, Probe,self.__class__,Probe)
_CVGsupport.Probe_swigregister(ProbePtr)

class Probelist(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Probelist, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Probelist, name)
    __swig_setmethods__["count"] = _CVGsupport.Probelist_count_set
    __swig_getmethods__["count"] = _CVGsupport.Probelist_count_get
    if _newclass:count = property(_CVGsupport.Probelist_count_get,_CVGsupport.Probelist_count_set)
    __swig_setmethods__["m_probes"] = _CVGsupport.Probelist_m_probes_set
    __swig_getmethods__["m_probes"] = _CVGsupport.Probelist_m_probes_get
    if _newclass:m_probes = property(_CVGsupport.Probelist_m_probes_get,_CVGsupport.Probelist_m_probes_set)
    def __init__(self,*args):
        _swig_setattr(self, Probelist, 'this', apply(_CVGsupport.new_Probelist,args))
        _swig_setattr(self, Probelist, 'thisown', 1)
    def append(*args): return apply(_CVGsupport.Probelist_append,args)
    def __del__(self, destroy= _CVGsupport.delete_Probelist):
        try:
            if self.thisown: destroy(self)
        except: pass
    def get_Z(*args): return apply(_CVGsupport.Probelist_get_Z,args)
    def __repr__(self):
        return "<C Probelist instance at %s>" % (self.this,)

class ProbelistPtr(Probelist):
    def __init__(self,this):
        _swig_setattr(self, Probelist, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Probelist, 'thisown', 0)
        _swig_setattr(self, Probelist,self.__class__,Probelist)
_CVGsupport.Probelist_swigregister(ProbelistPtr)

class SeqMat(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SeqMat, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SeqMat, name)
    __swig_setmethods__["width"] = _CVGsupport.SeqMat_width_set
    __swig_getmethods__["width"] = _CVGsupport.SeqMat_width_get
    if _newclass:width = property(_CVGsupport.SeqMat_width_get,_CVGsupport.SeqMat_width_set)
    __swig_setmethods__["numgenomes"] = _CVGsupport.SeqMat_numgenomes_set
    __swig_getmethods__["numgenomes"] = _CVGsupport.SeqMat_numgenomes_get
    if _newclass:numgenomes = property(_CVGsupport.SeqMat_numgenomes_get,_CVGsupport.SeqMat_numgenomes_set)
    __swig_setmethods__["gamma"] = _CVGsupport.SeqMat_gamma_set
    __swig_getmethods__["gamma"] = _CVGsupport.SeqMat_gamma_get
    if _newclass:gamma = property(_CVGsupport.SeqMat_gamma_get,_CVGsupport.SeqMat_gamma_set)
    __swig_setmethods__["gamma_wt"] = _CVGsupport.SeqMat_gamma_wt_set
    __swig_getmethods__["gamma_wt"] = _CVGsupport.SeqMat_gamma_wt_get
    if _newclass:gamma_wt = property(_CVGsupport.SeqMat_gamma_wt_get,_CVGsupport.SeqMat_gamma_wt_set)
    __swig_setmethods__["theta_wt"] = _CVGsupport.SeqMat_theta_wt_set
    __swig_getmethods__["theta_wt"] = _CVGsupport.SeqMat_theta_wt_get
    if _newclass:theta_wt = property(_CVGsupport.SeqMat_theta_wt_get,_CVGsupport.SeqMat_theta_wt_set)
    __swig_setmethods__["zeta_wt"] = _CVGsupport.SeqMat_zeta_wt_set
    __swig_getmethods__["zeta_wt"] = _CVGsupport.SeqMat_zeta_wt_get
    if _newclass:zeta_wt = property(_CVGsupport.SeqMat_zeta_wt_get,_CVGsupport.SeqMat_zeta_wt_set)
    __swig_setmethods__["deltamin"] = _CVGsupport.SeqMat_deltamin_set
    __swig_getmethods__["deltamin"] = _CVGsupport.SeqMat_deltamin_get
    if _newclass:deltamin = property(_CVGsupport.SeqMat_deltamin_get,_CVGsupport.SeqMat_deltamin_set)
    __swig_setmethods__["tdeltamin"] = _CVGsupport.SeqMat_tdeltamin_set
    __swig_getmethods__["tdeltamin"] = _CVGsupport.SeqMat_tdeltamin_get
    if _newclass:tdeltamin = property(_CVGsupport.SeqMat_tdeltamin_get,_CVGsupport.SeqMat_tdeltamin_set)
    __swig_setmethods__["beta"] = _CVGsupport.SeqMat_beta_set
    __swig_getmethods__["beta"] = _CVGsupport.SeqMat_beta_get
    if _newclass:beta = property(_CVGsupport.SeqMat_beta_get,_CVGsupport.SeqMat_beta_set)
    __swig_setmethods__["m_bg"] = _CVGsupport.SeqMat_m_bg_set
    __swig_getmethods__["m_bg"] = _CVGsupport.SeqMat_m_bg_get
    if _newclass:m_bg = property(_CVGsupport.SeqMat_m_bg_get,_CVGsupport.SeqMat_m_bg_set)
    __swig_setmethods__["mask"] = _CVGsupport.SeqMat_mask_set
    __swig_getmethods__["mask"] = _CVGsupport.SeqMat_mask_get
    if _newclass:mask = property(_CVGsupport.SeqMat_mask_get,_CVGsupport.SeqMat_mask_set)
    __swig_setmethods__["theta"] = _CVGsupport.SeqMat_theta_set
    __swig_getmethods__["theta"] = _CVGsupport.SeqMat_theta_get
    if _newclass:theta = property(_CVGsupport.SeqMat_theta_get,_CVGsupport.SeqMat_theta_set)
    __swig_setmethods__["ntheta"] = _CVGsupport.SeqMat_ntheta_set
    __swig_getmethods__["ntheta"] = _CVGsupport.SeqMat_ntheta_get
    if _newclass:ntheta = property(_CVGsupport.SeqMat_ntheta_get,_CVGsupport.SeqMat_ntheta_set)
    __swig_setmethods__["ptheta"] = _CVGsupport.SeqMat_ptheta_set
    __swig_getmethods__["ptheta"] = _CVGsupport.SeqMat_ptheta_get
    if _newclass:ptheta = property(_CVGsupport.SeqMat_ptheta_get,_CVGsupport.SeqMat_ptheta_set)
    __swig_setmethods__["priortheta"] = _CVGsupport.SeqMat_priortheta_set
    __swig_getmethods__["priortheta"] = _CVGsupport.SeqMat_priortheta_get
    if _newclass:priortheta = property(_CVGsupport.SeqMat_priortheta_get,_CVGsupport.SeqMat_priortheta_set)
    __swig_setmethods__["zeta_1"] = _CVGsupport.SeqMat_zeta_1_set
    __swig_getmethods__["zeta_1"] = _CVGsupport.SeqMat_zeta_1_get
    if _newclass:zeta_1 = property(_CVGsupport.SeqMat_zeta_1_get,_CVGsupport.SeqMat_zeta_1_set)
    __swig_setmethods__["nzeta_1"] = _CVGsupport.SeqMat_nzeta_1_set
    __swig_getmethods__["nzeta_1"] = _CVGsupport.SeqMat_nzeta_1_get
    if _newclass:nzeta_1 = property(_CVGsupport.SeqMat_nzeta_1_get,_CVGsupport.SeqMat_nzeta_1_set)
    __swig_setmethods__["pzeta_1"] = _CVGsupport.SeqMat_pzeta_1_set
    __swig_getmethods__["pzeta_1"] = _CVGsupport.SeqMat_pzeta_1_get
    if _newclass:pzeta_1 = property(_CVGsupport.SeqMat_pzeta_1_get,_CVGsupport.SeqMat_pzeta_1_set)
    __swig_setmethods__["priorzeta_1"] = _CVGsupport.SeqMat_priorzeta_1_set
    __swig_getmethods__["priorzeta_1"] = _CVGsupport.SeqMat_priorzeta_1_get
    if _newclass:priorzeta_1 = property(_CVGsupport.SeqMat_priorzeta_1_get,_CVGsupport.SeqMat_priorzeta_1_set)
    __swig_setmethods__["zeta_2"] = _CVGsupport.SeqMat_zeta_2_set
    __swig_getmethods__["zeta_2"] = _CVGsupport.SeqMat_zeta_2_get
    if _newclass:zeta_2 = property(_CVGsupport.SeqMat_zeta_2_get,_CVGsupport.SeqMat_zeta_2_set)
    __swig_setmethods__["nzeta_2"] = _CVGsupport.SeqMat_nzeta_2_set
    __swig_getmethods__["nzeta_2"] = _CVGsupport.SeqMat_nzeta_2_get
    if _newclass:nzeta_2 = property(_CVGsupport.SeqMat_nzeta_2_get,_CVGsupport.SeqMat_nzeta_2_set)
    __swig_setmethods__["pzeta_2"] = _CVGsupport.SeqMat_pzeta_2_set
    __swig_getmethods__["pzeta_2"] = _CVGsupport.SeqMat_pzeta_2_get
    if _newclass:pzeta_2 = property(_CVGsupport.SeqMat_pzeta_2_get,_CVGsupport.SeqMat_pzeta_2_set)
    __swig_setmethods__["priorzeta_2"] = _CVGsupport.SeqMat_priorzeta_2_set
    __swig_getmethods__["priorzeta_2"] = _CVGsupport.SeqMat_priorzeta_2_get
    if _newclass:priorzeta_2 = property(_CVGsupport.SeqMat_priorzeta_2_get,_CVGsupport.SeqMat_priorzeta_2_set)
    __swig_setmethods__["bg"] = _CVGsupport.SeqMat_bg_set
    __swig_getmethods__["bg"] = _CVGsupport.SeqMat_bg_get
    if _newclass:bg = property(_CVGsupport.SeqMat_bg_get,_CVGsupport.SeqMat_bg_set)
    __swig_setmethods__["bg_ng"] = _CVGsupport.SeqMat_bg_ng_set
    __swig_getmethods__["bg_ng"] = _CVGsupport.SeqMat_bg_ng_get
    if _newclass:bg_ng = property(_CVGsupport.SeqMat_bg_ng_get,_CVGsupport.SeqMat_bg_ng_set)
    __swig_setmethods__["joint"] = _CVGsupport.SeqMat_joint_set
    __swig_getmethods__["joint"] = _CVGsupport.SeqMat_joint_get
    if _newclass:joint = property(_CVGsupport.SeqMat_joint_get,_CVGsupport.SeqMat_joint_set)
    def __init__(self,*args):
        _swig_setattr(self, SeqMat, 'this', apply(_CVGsupport.new_SeqMat,args))
        _swig_setattr(self, SeqMat, 'thisown', 1)
    def get_best(*args): return apply(_CVGsupport.SeqMat_get_best,args)
    def find_matching_probes(*args): return apply(_CVGsupport.SeqMat_find_matching_probes,args)
    def calcdist(*args): return apply(_CVGsupport.SeqMat_calcdist,args)
    def calcKL(*args): return apply(_CVGsupport.SeqMat_calcKL,args)
    def scanbest(*args): return apply(_CVGsupport.SeqMat_scanbest,args)
    def scanbest2(*args): return apply(_CVGsupport.SeqMat_scanbest2,args)
    def score(*args): return apply(_CVGsupport.SeqMat_score,args)
    def set(*args): return apply(_CVGsupport.SeqMat_set,args)
    def get(*args): return apply(_CVGsupport.SeqMat_get,args)
    def get_c(*args): return apply(_CVGsupport.SeqMat_get_c,args)
    def get_theta(*args): return apply(_CVGsupport.SeqMat_get_theta,args)
    def get_gamma(*args): return apply(_CVGsupport.SeqMat_get_gamma,args)
    def getKL(*args): return apply(_CVGsupport.SeqMat_getKL,args)
    def matchstarts(*args): return apply(_CVGsupport.SeqMat_matchstarts,args)
    def matchstarts2(*args): return apply(_CVGsupport.SeqMat_matchstarts2,args)
    def score_probe(*args): return apply(_CVGsupport.SeqMat_score_probe,args)
    def EMstep(*args): return apply(_CVGsupport.SeqMat_EMstep,args)
    def setBg(*args): return apply(_CVGsupport.SeqMat_setBg,args)
    def setTheta(*args): return apply(_CVGsupport.SeqMat_setTheta,args)
    def setZeta(*args): return apply(_CVGsupport.SeqMat_setZeta,args)
    def setmask(*args): return apply(_CVGsupport.SeqMat_setmask,args)
    def loglikelihood(*args): return apply(_CVGsupport.SeqMat_loglikelihood,args)
    def __del__(self, destroy= _CVGsupport.delete_SeqMat):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C SeqMat instance at %s>" % (self.this,)

class SeqMatPtr(SeqMat):
    def __init__(self,this):
        _swig_setattr(self, SeqMat, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, SeqMat, 'thisown', 0)
        _swig_setattr(self, SeqMat,self.__class__,SeqMat)
_CVGsupport.SeqMat_swigregister(SeqMatPtr)


