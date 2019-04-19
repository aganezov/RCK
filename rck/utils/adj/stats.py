from collections import defaultdict


def merged_source_tally(adjacencies, extra_sources_field="supporting_sources"):
    counter = defaultdict(int)
    for adj in adjacencies:
        sources_string = adj.extra.get(extra_sources_field, None)
        if sources_string is None:
            counter[sources_string] += 1
        else:
            sources = tuple(sorted(sources_string.split(",")))
            counter[sources] += 1
    return counter
