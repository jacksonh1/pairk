# from enum import Enum
# class MyEnum(Enum):
#     """Docstring for MyEnum."""
#     FIRST_ENUM = "some_value"
#     SECOND_ENUM = "some_other_value"

import copy
from typing import Literal, Optional, Union

from attrs import define, field, validators

import local_env_variables.env_variables as env


@define
class FilterConf:
    """sequence filtering parameters
    
    Attributes:
    `min_fraction_shorter_than_query`: float,
        the minimum fraction of the query sequence length that each orthogroup sequence must be.
        Default: 0.5
    """
    min_fraction_shorter_than_query: float = field(
        default=0.5, validator=validators.and_(validators.le(1), validators.ge(0))
    )


@define
class OGSelectConf:
    '''
    orthogroup selection parameters

    Attributes:
    `OG_selection_method`: Union[str, Literal["level_name"]]
        The method used to select the orthogroup.
        `level_name`: select the orthogroup with the given level name
        Default: "level_name"
    `OG_level_name`: str,
        the level name to use if OG_selection_method is "level_name".
        Default: "Vertebrata"
    '''
    OG_selection_method: Union[str, Literal["level_name"]] = field(
        default="level_name", validator=validators.in_(["level_name"])
    )
    OG_level_name: str = field(default="Vertebrata")


@define
class LDOSelectConf:
    """
    least divergent ortholog selection parameters

    Attributes:
    `LDO_selection_method`: Union[str, Literal["msa_by_organism", "alfpy_google_distance", "pairwise", "msa"]]
        The method used to select the least divergent orthologs. For each method, 
        the paralog in each organism with the highest percent identity (PID) to 
        the query sequence is chosen as the LDO. The only exception is `alfpy_google_distance`,
        where the similarity between the sequences is used instead.
        Methods:
        `msa`: Aligns all of the sequences in the ortholog group. 
        `msa_by_organism`: performs a separate alignment of the paralogs in each organism with the query sequence.
        `alfpy_google_distance`: Uses an alignment free word-based method to calculate 
            the similarity (1-distance) between the query sequence and each paralog.
        `pairwise`: performs a pairwise alignment between the query sequence and each ortholog using BioPython.
        Default: "alfpy_google_distance"
    `_LDO_msa_exe`: str,
        path to the msa executable. only used if LDO_selection_method is "msa" or "msa_by_organism". intended to be specified in the .env file
        Default: path specified in the .env file
    `_LDO_msa_additional_args`: str,
        additional arguments to pass to the msa executable. only used if LDO_selection_method is "msa" or "msa_by_organism". intended to be specified in the .env file
        Default: arguments specified in the .env file
    `LDO_msa_threads`: int,
        number of threads to use for msa. only used if LDO_selection_method is "msa" or "msa_by_organism".
        Defaults: 8

    if LDO_selection_method doesn't require an msa (i.e. `alfpy_google_distance` or `pairwise`),
     then LDO_msa_exe and LDO_msa_threads are ignored
    """
    LDO_selection_method: Union[
        str, Literal["msa_by_organism", "alfpy_google_distance", "pairwise", "msa"]
    ] = field(
        default="alfpy_google_distance",
        validator=validators.in_(
            ["msa_by_organism", "alfpy_google_distance", "pairwise", "msa"]
        ),
    )
    _LDO_mafft_exe: str = field(default=env.MAFFT_EXECUTABLE)
    _LDO_mafft_additional_args: str = field(default=env.MAFFT_ADDITIONAL_ARGUMENTS)
    LDO_mafft_threads: int = field(default=8, converter=int)


@define
class AlignConf:
    align: bool = field(default=False, converter=bool)
    n_align_threads: int = field(default=8, converter=int)
    _mafft_exe: str = field(default=env.MAFFT_EXECUTABLE)
    _mafft_additional_args: str = field(default=env.MAFFT_ADDITIONAL_ARGUMENTS)


@define
class PipelineParams:
    filter_params: FilterConf = field(default=FilterConf())
    og_select_params: OGSelectConf = field(default=OGSelectConf())
    ldo_select_params: LDOSelectConf = field(default=LDOSelectConf())
    align_params: AlignConf = field(default=AlignConf())
    _cd_hit_exe: str = field(default=env.CD_HIT_EXECUTABLE)
    _cd_hit_additional_args: str = field(default=env.CD_HIT_ADDITIONAL_ARGUMENTS)
    main_output_folder: str = field(default="./processed_odb_groups_output")
    write_files: bool = field(default=True)
    overwrite: bool = field(default=False)

    @classmethod
    def from_dict(cls, d):
        d = copy.deepcopy(d)
        return cls(
            filter_params=FilterConf(**d.pop("filter_params", {})),
            og_select_params=OGSelectConf(**d.pop("og_select_params", {})),
            ldo_select_params=LDOSelectConf(**d.pop("ldo_select_params", {})),
            align_params=AlignConf(**d.pop("align_params", {})),
            **d,
        )


# ==============================================================================
# // test
# ==============================================================================

# DEFAULT_PARAM_DICT = {
#     "filter_params": {
#         "min_fraction_shorter_than_query": 0.5,
#     },
#     "og_select_params": {
#         "OG_selection_method": "level_name",
#         "OG_level_name": "Eukaryota",
#     },
#     "ldo_select_params": {
#         "LDO_selection_method": "alfpy_google_distance",
#     },
#     "align_params": {
#         "align": False,
#         "n_align_threads": 8,
#     },
#     "write_files": True,
#     "main_output_folder": "./orthoDB_analysis",
# }
# config = PipelineParams.from_dict(DEFAULT_PARAM_DICT)
# print(config)