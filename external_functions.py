from typing import Iterable, List, Tuple, Optional

import pandas as pd



import bisect

        
def make_arrow_from_SE_coords(start, end, base_ycoord, top_ycoord, arrow_direction='forward'):
    """"""
    if arrow_direction=='forward':
        br = start + 0.85*(end-start)
        mid=base_ycoord+0.5
        xvals=[start,            br,             end,      br,           start,         start        ]
        yvals=[base_ycoord,      base_ycoord,    mid,     top_ycoord,    top_ycoord,    base_ycoord]
        return xvals, yvals
    elif arrow_direction=='reverse':
        return make_arrow_from_SE_coords(end, start, base_ycoord, top_ycoord)



def map_snps_to_genes(snps_df, annot_df): 
    
    colnames  = [
        "protein_ID",
        "common_name",
        "product",
        "start",
        "end",
        "strand",
        "contig",
        'location', 
        'snp_id',
        'snp_position'
        ]
    cols_d = {c:[] for c in colnames}
    starts = list(annot_df['start'])
    ends = list(annot_df['end'])
    
    
    
    for ind, row in snps_df.iterrows():
        start_index = bisect.bisect_left(starts, row['ref_pos'])
        end_index   = bisect.bisect_left(ends, row['ref_pos'])
        s = annot_df.iloc[end_index]
        cols_d['snp_id'].append(ind)
        cols_d['snp_position'].append(row['ref_pos'])
        for col in colnames:
            if col in s.index:
                cols_d[col].append(s[col])
            elif col == 'location':
                if start_index-1 == end_index: 
                    location = 'within_gene'
                    cols_d[col].append(location)
                else:
                    location = 'intergenic'
                    cols_d[col].append(location)
    return snps_df.join(pd.DataFrame(cols_d).set_index('snp_id'))  


def make_snp_df(snps_fp):
    snps_df = pd.read_csv(
        snps_fp,
        sep="\t",
        header=None
    )
    
    snps_df.columns = [
        "ref_pos",
        "ref_base",
        "query_base",
        "query_pos",
        "ref_pos_aln",
        "query_pos_aln",
        "ref_contig_len",
        "query_contig_len",
        "ref_strand",
        "query_strand",
        "ref_contig",
        "query_contig",
    ]
    return snps_df


def unaligned_segments(
    starts: Iterable[int],
    ends: Iterable[int],
    contig_len: int,
    *,
    assume_half_open: bool = True,
    merge_adjacent: bool = True,
    min_len: int = 1
) -> List[Tuple[int, int]]:
    """
    Return reference segments with no alignments.

    Parameters
    ----------
    starts, ends
        Coordinates for aligned segments on the reference.
        By default these are treated as 0-based, half-open intervals [start, end).
    contig_len
        Total length of the contig.
    assume_half_open
        If False, treat input as closed intervals [start, end] and convert to half-open.
    merge_adjacent
        If True, merge intervals that touch (e.g., [10,20) and [20,30)) as continuous.
    min_len
        Only return unaligned segments with length >= min_len.

    Returns
    -------
    List of (start, end) half-open intervals [start, end) that are unaligned.
    """
    if contig_len < 0:
        raise ValueError("contig_len must be >= 0")

    # Build and sanitize intervals
    intervals: List[Tuple[int, int]] = []
    for s, e in zip(starts, ends):
        if s is None or e is None:
            continue
        s = int(s)
        e = int(e)
        if not assume_half_open:
            # convert closed [s, e] to half-open [s, e+1]
            e += 1

        # normalize if reversed
        if e < s:
            s, e = e, s

        # clamp to reference bounds
        s = max(0, min(s, contig_len))
        e = max(0, min(e, contig_len))

        # skip empty
        if e <= s:
            continue

        intervals.append((s, e))

    if not intervals:
        return [(0, contig_len)] if contig_len >= min_len else []

    # Sort and merge
    intervals.sort()
    merged: List[Tuple[int, int]] = []
    cur_s, cur_e = intervals[0]

    for s, e in intervals[1:]:
        touch = (s <= cur_e) if merge_adjacent else (s < cur_e)
        if touch:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))

    # Compute complement
    gaps: List[Tuple[int, int]] = []
    prev_end = 0
    for s, e in merged:
        if s > prev_end:
            if (s - prev_end) >= min_len:
                gaps.append((prev_end, s))
        prev_end = max(prev_end, e)

    if contig_len > prev_end and (contig_len - prev_end) >= min_len:
        gaps.append((prev_end, contig_len))

    return gaps



def parse_show_coords(show_coords_fp):
    df = pd.read_csv(show_coords_fp,skiprows=2).drop(0, axis=0)#.str.split('|')
    df.columns = ['merged']
    t_df = df['merged'].str.split('|')
    t_df = pd.DataFrame(t_df.to_dict()).transpose()

    parsed_df = pd.DataFrame()
    parsed_df['S1,E1'] = t_df[0].apply(lambda x: [x for x in x.split(' ') if x] )
    parsed_df['S1'] = parsed_df.apply(lambda x: x['S1,E1'][0], axis=1)
    parsed_df['E1'] = parsed_df.apply(lambda x: x['S1,E1'][1], axis=1)
    parsed_df = parsed_df.drop('S1,E1', axis=1)
    
    parsed_df['S2,E2'] = t_df[1].apply(lambda x: [x for x in x.split(' ') if x] )
    parsed_df['S2'] = parsed_df.apply(lambda x: x['S2,E2'][0], axis=1)
    parsed_df['E2'] = parsed_df.apply(lambda x: x['S2,E2'][1], axis=1)
    parsed_df = parsed_df.drop('S2,E2', axis=1)
    
    parsed_df['len1len2'] = t_df[2].apply(lambda x: [x for x in x.split(' ') if x] )
    parsed_df['len1'] = parsed_df.apply(lambda x: x['len1len2'][0], axis=1)
    parsed_df['len2'] = parsed_df.apply(lambda x: x['len1len2'][1], axis=1)
    parsed_df = parsed_df.drop('len1len2', axis=1)
    
    parsed_df['Identity'] = t_df[3].replace(' ', '')
    parsed_df['Tag1'] = t_df[4].apply(lambda x: x.split('\t')[0])
    parsed_df['Tag2'] = t_df[4].apply(lambda x: x.split('\t')[1])
    for col in ['S1', 'E1', 'S2', 'E2', 'len1', 'len2']:
        parsed_df[col] = parsed_df[col].astype(int)
    parsed_df['same_orientation'] = (parsed_df['S2'] < parsed_df['E2'] )
    return parsed_df