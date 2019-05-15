from collections import defaultdict


def get_size_bins(bins_strs):
    result = [-3000000000]
    for str_value in bins_strs:
        result.append(int(float(str_value)))
    result.append(3000000000)
    return result


def get_adj_size(adjacency, size_extra_field="svlen", size_extra_field_abs=True, size_extra_seq_field=None):
    adj_size = None
    try:
        adj_size = int(float(adjacency.extra[size_extra_field]))
        if size_extra_field_abs:
            adj_size = abs(adj_size)
    except (KeyError, ValueError):
        pass
    if adj_size is None:
        try:
            adj_size = len(adjacency.extra[size_extra_seq_field])
        except (KeyError, ValueError):
            pass
    if adj_size is None:
        adj_size = adjacency.distance_non_hap
    return adj_size


def merged_source_tally(adjacencies, bins=None, extra_sources_field="supporting_sources", size_extra_field="svlen", size_extra_field_abs=True, size_extra_seq_field=None):
    if bins is None:
        bins = [-3000000000, 3000000000]
    result = defaultdict(lambda: defaultdict(int))
    for adj in adjacencies:
        adj_size = get_adj_size(adjacency=adj, size_extra_field=size_extra_field, size_extra_field_abs=size_extra_field_abs, size_extra_seq_field=size_extra_seq_field)
        sources_string = adj.extra.get(extra_sources_field, None)
        target_bin = None
        for bin in bins:
            if adj_size < bin:
                target_bin = bin
                break
        source = None if sources_string is None else tuple(sorted(sources_string.split(",")))
        result[source][target_bin] += 1
    return result
