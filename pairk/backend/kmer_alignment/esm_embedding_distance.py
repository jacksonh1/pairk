import copy
import pairk.backend.tools.sequence_utils as tools
import pairk.backend.tools.pairwise_tools as pairwise_tools
import pairk.backend.exceptions as _exceptions
import pairk.backend.tools.esm_tools as esm_tools
import torch
from collections import defaultdict
import torch


def get_subsequences(indices, input_string, k):
    """Get subsequences of ortho string based on best indices"""
    char_list = list(input_string)
    start_indices = indices.view(-1, 1)
    end_indices = start_indices + k
    max_length = len(char_list)
    subsequences = [
        "".join(char_list[start:end]) for start, end in zip(start_indices, end_indices)
    ]
    return subsequences


def run_pairwise_kmer_emb_aln(
    query_id: str,
    embedding_dict: dict,
    k: int,
):
    """
    This function was originally too slow but was vectorized/optimized by
    Foster Birnbaum. Thanks Foster!
    """
    # get the query sequence and remove it from the embedding_dict
    ref_seq_str, ref_seq_embedding = embedding_dict.pop(query_id)
    kmers = tools.gen_kmers(ref_seq_str, k)
    positions = list(range(len(ref_seq_str) - (k - 1)))
    score_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(embedding_dict.keys())
    )
    orthokmer_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(embedding_dict.keys())
    )
    pos_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(embedding_dict.keys())
    )
    score_df.loc[positions, "query_kmer"] = kmers
    orthokmer_df.loc[positions, "query_kmer"] = kmers
    pos_df.loc[positions, "query_kmer"] = kmers

    # Make expanded ref tensor
    expand_inds_ref = torch.arange(k).view(1, -1) + torch.arange(
        ref_seq_embedding.shape[0]
    ).view(-1, 1)
    expand_inds_ref[ref_seq_embedding.shape[0] - (k - 1) :] = 0
    expand_inds_ref = (
        expand_inds_ref.unsqueeze(-1)
        .expand(-1, -1, ref_seq_embedding.shape[1])
        .to(dtype=torch.int64)
    )
    expand_ref = ref_seq_embedding.unsqueeze(1).expand(-1, expand_inds_ref.shape[1], -1)
    expand_ref = torch.gather(expand_ref, 0, expand_inds_ref)
    expand_ref = expand_ref[: ref_seq_embedding.shape[0] - (k - 1)].reshape(
        -1, k * expand_ref.shape[2]
    )

    # for each ortholog sequence
    for ortholog_id, v in embedding_dict.items():
        ortholog_seq = v[0]
        ortholog_embedding = v[1]
        if ortholog_seq == "no idr" or ortholog_embedding == "no idr":
            orthokmer_df.loc[positions, ortholog_id] = "-" * k
            continue
        if len(ortholog_seq) < k:
            orthokmer_df.loc[positions, ortholog_id] = "-" * k
            continue

        # Make expanded ortholog tensor
        expand_inds_ortho = torch.arange(k).view(1, -1) + torch.arange(
            ortholog_embedding.shape[0]
        ).view(-1, 1)
        expand_inds_ortho[ortholog_embedding.shape[0] - (k - 1) :] = 0
        expand_inds_ortho = (
            expand_inds_ortho.unsqueeze(-1)
            .expand(-1, -1, ortholog_embedding.shape[1])
            .to(dtype=torch.int64)
        )
        expand_ortho = ortholog_embedding.unsqueeze(1).expand(
            -1, expand_inds_ortho.shape[1], -1
        )
        expand_ortho = torch.gather(expand_ortho, 0, expand_inds_ortho)
        expand_ortho = expand_ortho[: ortholog_embedding.shape[0] - (k - 1)].reshape(
            -1, k * expand_ortho.shape[2]
        )

        # Calculate pairwise distances and get stats
        pairwise_dists = torch.cdist(
            expand_ref, expand_ortho, p=2
        )  # Optional: compute_mode='donot_use_mm_for_euclid_dist'
        min_dists, min_dists_pos = torch.min(pairwise_dists, dim=-1)
        score_df.loc[positions, ortholog_id] = min_dists.cpu().numpy()
        pos_df.loc[positions, ortholog_id] = min_dists_pos.cpu().numpy()
        orthokmer_df.loc[positions, ortholog_id] = get_subsequences(
            min_dists_pos, ortholog_seq, k
        )
    return score_df, orthokmer_df, pos_df


def get_idr_embeddings(
    seq_str: str,
    idr_start: int,
    idr_end: int,
    mod: esm_tools.ESM_Model,
    device="cuda",
):
    idr_str = seq_str[idr_start : idr_end + 1]
    if len(idr_str) == 0:
        return "no idr", "no idr"
        # print("no idr")
        # return None, None
    orth_tensor = mod.encode(seq_str, device=device)
    idr_ortho_tensor = orth_tensor[idr_start + 1 : idr_end + 2, :]  # type: ignore # +1 to account for the start token
    return idr_str, idr_ortho_tensor


# def get_idr_embedding_dict(
#     full_length_sequence_dict: dict[str, str],
#     idr_position_map: dict[str, list[int]],
#     mod: esm_tools.ESM_Model,
#     device="cuda",
# ):
#     embedding_dict = defaultdict(list)
#     # get dictionary keys and values as 2 separate lists
#     # fl_seqs = [(i, seq) for i, seq in full_length_sequence_dict.items()]
#     seq_ids, seqs = zip(*full_length_sequence_dict.items())
#     seq_embeddings = mod.encode_multiple_seqs(seqs, device=device)
#     for i, seq, seq_embedding in zip(seq_ids, seqs, seq_embeddings):
#         idrst, idrend = idr_position_map[i][0], idr_position_map[i][1]
#         idr_str = seq[idrst : idrend + 1]
#         if len(idr_str) == 0:
#             embedding_dict[i].append("no idr")
#             embedding_dict[i].append("no idr")
#             continue
#         idr_ortho_tensor = seq_embedding[idrst + 1 : idrend + 2, :]  # type: ignore # +1 to account for the start token
#         embedding_dict[i].append(idr_str)
#         embedding_dict[i].append(idr_ortho_tensor)
#     return embedding_dict


def get_idr_embedding_dict(
    full_length_sequence_dict: dict[str, str],
    idr_position_map: dict[str, list[int]],
    mod: esm_tools.ESM_Model,
    device="cuda",
):
    embedding_dict = defaultdict(list)
    for i, seq in full_length_sequence_dict.items():
        idr_str, idr_tensor = get_idr_embeddings(
            seq,
            idr_position_map[i][0],
            idr_position_map[i][1],
            mod,
            device,
        )
        embedding_dict[i].append(idr_str)
        embedding_dict[i].append(idr_tensor)
    return embedding_dict


def slice_idr_embedding_dict(
    full_length_sequence_dict: dict[str, str],
    idr_position_map: dict[str, list[int]],
    precomputed_embeddings: dict[str, torch.Tensor],
):
    embedding_dict = defaultdict(list)
    for i, seq in full_length_sequence_dict.items():
        idr_start = idr_position_map[i][0]
        idr_end = idr_position_map[i][1]
        idr_str = seq[idr_start : idr_end + 1]
        if len(idr_str) == 0:
            embedding_dict[i].append("no idr")
            embedding_dict[i].append("no idr")
            continue
        idr_ortho_tensor = precomputed_embeddings[i][idr_start + 1 : idr_end + 2, :]  # type: ignore # +1 to account for the start token
        embedding_dict[i].append(idr_str)
        embedding_dict[i].append(idr_ortho_tensor)
    return embedding_dict


def pairk_alignment_embedding_distance(
    full_length_sequence_dict: dict[str, str],
    idr_position_map: dict[str, list[int]],
    query_id: str,
    k: int,
    mod: esm_tools.ESM_Model,
    device: str = "cuda",
    precomputed_embeddings: None | dict[str, torch.Tensor] = None,
):
    """run pairwise k-mer alignment method using sequence embeddings from the
    ESM2 protein large language model to find the best k-mer matches from each
    homolog. If a ortholog IDR is shorter than the k-mer, a string of "-"
    characters ("-"\\*k) is assigned as the best matching ortholog k-mer for that
    ortholog

    **Note**: if there are multiple top-scoring matches, only one is returned.

    Sequence embeddings are calculated for each full length sequence in the
    input dictionary. The `idr_position_map` dictionary is used to extract the
    IDR and the IDR embeddings from each sequence. The Euclidean distance is
    calculated between each query k-mer embedding slice and each ortholog k-mer
    embedding slice to find the best matching ortholog k-mer from each ortholog.


    Parameters
    ----------
    full_length_sequence_dict : dict[str, str]
        input sequences in dictionary format with the key being the sequence id and
        the value being the sequence as a string
    idr_position_map : dict[str, list[int]]
        a dictionary where the keys are the sequence ids in `full_length_sequence_dict`
        and the values are the start and end positions of the IDR in the sequence
        (using python indexing). This is used to slice out the IDR
        embeddings/sequences from the full-length embeddings/sequences.
    query_id : str
        the id of the query sequence within the `full_length_sequence_dict` dictionary
        and the `idr_position_map` dictionary. The query id must be present in
        both dictionaries.
    k : int
        the length of the k-mers to use for the alignment
    mod : esm_tools.ESM_Model
        ESM2 model used to generate the embeddings
    device : str, optional
        whether to use cuda or cpu for pytorch, must be either "cpu" or "cuda",
        by default "cuda". If "cuda" fails, it will default to "cpu". This
        argument is passed to the `esm_tools.ESM_Model.encode` method.
    precomputed_embeddings : None | dict[str, torch.Tensor], optional
        a dictionary where the keys are the sequence ids in `full_length_sequence_dict`
        and the values are the precomputed embeddings for each sequence. If this
        is provided, the function will use these embeddings instead of computing
        them.

    Returns
    -------
    pairwise_tools.PairkAln
        an object containing the alignment results. See the `pairk.PairkAln`
        class for more information.
    """
    _exceptions.check_queryid_in_idr_dict(full_length_sequence_dict, query_id)
    _exceptions.check_queryid_in_idr_dict(idr_position_map, query_id)
    assert set(full_length_sequence_dict.keys()) == set(
        idr_position_map.keys()
    ), "Keys in full_length_dict and idr_position_map must be the same"
    full_length_dict = copy.deepcopy(full_length_sequence_dict)
    if precomputed_embeddings is None:
        embedding_dict = get_idr_embedding_dict(
            full_length_dict, idr_position_map, mod, device
        )
    else:
        _exceptions.check_queryid_in_idr_dict(precomputed_embeddings, query_id)
        assert set(full_length_sequence_dict.keys()) == set(
            precomputed_embeddings.keys()
        ), "Keys in full_length_dict and precomputed_embeddings must be the same"
        embedding_dict = slice_idr_embedding_dict(
            full_length_sequence_dict=full_length_dict,
            idr_position_map=idr_position_map,
            precomputed_embeddings=precomputed_embeddings,
        )
    score_df, orthokmer_df, pos_df = run_pairwise_kmer_emb_aln(
        query_id,
        embedding_dict,
        k,
    )
    return pairwise_tools.PairkAln(
        orthokmer_df=orthokmer_df,
        pos_df=pos_df,
        score_df=score_df,
    )
