#!/usr/bin/env python3
"""The last job of a campaign whose fills were queued as plain jobs: merges and the validated collection.

usage: campaign_finish.py BINARY --work WORK --out OUT [--parallel N] [--quiet-minutes M]
(called as `tree_to_hists.sh campaign-finish --work WORK --out OUT`, queued together with the fills)

A campaign is a graph (WORK/campaign_graph.json): fills, merges whose inputs are fills or merges, and the collection.
A fill is a group of small tree parts (fill-group G<n>) or one range of a large part (fill F<n>: a part with more
events than Campaign.MaxEventsPerJob is filled range by range and merged by a node of its own; DATA has thousands).
Where the graph cannot be handed to DAGMan, every fill and this job go into one submission:
  1. waits until every fill has its receipt (done.json beside a group's package, or in a range's folder);
  2. when at least half are done and no new one has arrived for M minutes (default 30), runs the open fills here,
     N at a time. A fill publishes whole, so a job that finishes later only finds the result;
  3. runs the merges in the order of the graph, N at a time (a merge that fails is tried three times), then the
     collection. A merge whose package and sidecar exist is not run again, so this job can be started again after an
     interruption.
Prints one line per step; exit status 0 only when the collection is there.
"""
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor


def option(name, default=None):
    if name in sys.argv:
        return sys.argv[sys.argv.index(name) + 1]
    if default is None:
        sys.exit(f'campaign_finish.py: missing {name}')
    return default


def say(*words):
    print(time.strftime('%H:%M:%S'), *words, flush=True)


def run(binary, arguments, tries=1):
    for attempt in range(tries):
        code = subprocess.call([binary] + arguments)
        if code == 0:
            return True
        say('failed', *arguments[:1], *arguments[-2:], f'(attempt {attempt + 1} of {tries}, exit {code})')
        time.sleep(20)
    return False


def main():
    binary = sys.argv[1]
    work, out = option('--work'), option('--out')
    parallel, quiet = int(option('--parallel', '4')), 60 * int(option('--quiet-minutes', '30'))
    nodes = {node['id']: node for node in json.load(open(f'{work}/campaign_graph.json'))['nodes']}
    fills = [n for n in nodes.values() if n['operation'] in ('fill', 'fill-group')]
    merges = [n for n in nodes.values() if n['operation'] == 'merge-group']
    if not fills:
        sys.exit('campaign_finish.py: the graph has no fills; plan the campaign first')

    def receipt(node):  # the file that says a fill is complete: beside a group's packages, or in a range's folder
        if node['operation'] == 'fill-group':
            return f"{work}/chunks/{int(node['id'][1:]) // 1000}/{node['id']}/done.json"
        return node['output'] + '/done.json'

    remaining, last_count, last_change = list(fills), -1, time.time()
    while True:
        remaining = [node for node in remaining if not os.path.exists(receipt(node))]
        done = len(fills) - len(remaining)
        if done != last_count:
            say(f'{done} of {len(fills)} fills done')
            last_count, last_change = done, time.time()
        if not remaining:
            break
        if 2 * done >= len(fills) and time.time() - last_change > quiet:
            say(f'running the {len(remaining)} open fills here')
            with ThreadPoolExecutor(parallel) as pool:
                list(pool.map(lambda node: os.path.exists(receipt(node)) or run(binary, node['argv']), remaining))
            last_change = time.time()
            continue
        time.sleep(60)

    def merged(node):
        return os.path.exists(node['output']) and os.path.exists(node['output'][:-5] + '.json')

    finished = {node['id'] for node in fills} | {node['id'] for node in merges if merged(node)}
    if len(finished) > len(fills):
        say(f'{len(finished) - len(fills)} of {len(merges)} merges are there already')
    while len(finished) < len(fills) + len(merges):
        level = [node for node in merges if node['id'] not in finished and all(p in finished for p in node['parents'])]
        if not level:
            sys.exit('campaign_finish.py: the merge graph does not close')
        say(f'merging {len(level)} nodes')
        with ThreadPoolExecutor(parallel) as pool:
            results = list(pool.map(lambda node: run(binary, node['argv'], tries=3), level))
        if not all(results):
            sys.exit('campaign_finish.py: a merge failed three times: ' + ' '.join(n['id'] for n, ok in zip(level, results) if not ok)[:2000])
        finished.update(node['id'] for node in level)
    if not os.path.exists(out + '/collection.yaml') and not run(binary, ['collection', '--work', work, '--out', out]):
        sys.exit('campaign_finish.py: the collection failed')
    say('CAMPAIGN_DONE', out)


if __name__ == '__main__':
    main()
