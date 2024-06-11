from localcider.sequenceParameters import SequenceParameters


def get_local_cider_features(sequence_str):
    Seqob = SequenceParameters(sequence_str)
    d = {
        "localcider - PPII_propensity": Seqob.get_PPII_propensity(),
        "localcider - uversky_hydropathy": Seqob.get_uversky_hydropathy(),
        "localcider - mean_hydropathy": Seqob.get_mean_hydropathy(),
        "localcider - mean_net_charge": Seqob.get_mean_net_charge(),
        "localcider - Omega": Seqob.get_Omega(),
        "localcider - kappa": Seqob.get_kappa(),
        "localcider - fraction_expanding": Seqob.get_fraction_expanding(),
        "localcider - fraction_positive": Seqob.get_fraction_positive(),
        "localcider - fraction_negative": Seqob.get_fraction_negative(),
        "localcider - countNeut": Seqob.get_countNeut(),
        "localcider - countPos": Seqob.get_countPos(),
        "localcider - countNeg": Seqob.get_countNeg(),
        "localcider - isoelectric_point": Seqob.get_isoelectric_point(),
        "localcider - NCPR": Seqob.get_NCPR(),
    }
    return d
