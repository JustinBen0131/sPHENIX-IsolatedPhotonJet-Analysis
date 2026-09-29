"""Read-only runtime integrity check: python -m eventtree check FILE.root."""
import argparse,json
from .check import check

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['check']);p.add_argument('file')
    a=p.parse_args();result=check(a.file)
    print(json.dumps(result,indent=2));return 0 if result['ok'] else 1

if __name__=='__main__':raise SystemExit(main())
