"""SPRITE and GAM experiment sets of the 4D Nucleome Data Portal by submitting laboratory (Introduction of the
manuscript: 17 of 18 experiment sets from the laboratories that developed the methods, queried on 2026-10-28).

Usage: python 4dn_multiway_labs.py
"""

import collections
import json
import time
import urllib.parse
import urllib.request

TYPES = ['DNA SPRITE', 'RNA-DNA SPRITE', 'scSPRITE', 'GAM', 'multiplex-GAM']


def main():
    for experiment_type in TYPES:
        query = urllib.parse.urlencode({'type': 'ExperimentSetReplicate', 'limit': 'all', 'field': 'lab.display_title',
                                        'experiments_in_set.experiment_type.display_title': experiment_type})
        request = urllib.request.Request('https://data.4dnucleome.org/search/?' + query,
                                         headers={'Accept': 'application/json'})
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                graph = json.load(response).get('@graph', [])
        except urllib.error.HTTPError as error:
            graph = [] if error.code == 404 else None  # the portal answers 404 for an empty result
        labs = collections.Counter(x.get('lab', {}).get('display_title', '?') for x in graph or [])
        print('{}\t{}\t{}'.format(experiment_type, len(graph or []), '; '.join(
            '{} ({})'.format(lab, n) for lab, n in labs.most_common())))
        time.sleep(1)


if __name__ == '__main__':
    main()
