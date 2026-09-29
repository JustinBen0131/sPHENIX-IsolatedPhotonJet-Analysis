"""Compact saved-source certificates for complete-group publication.

Full payload/inverse/nominal validation runs in each bounded source worker.
Publication rehashes those immutable files and certificates, then checks a
compact DATA identity stream; it does not rerun the full analysis per file.
"""
import json
from pathlib import Path
import sqlite3

import numpy as np
import uproot

from .integrity import sha256
from .tables import require


def create(path,proof,index_path):
    require(proof.get('status')=='PASS' and proof.get('proof',{}).get('status')=='PASS','source conservation proof missing')
    require(proof.get('nominal',{}).get('status')=='PASS','source nominal proof missing')
    require(proof.get('schema') in ('RJLogicalSourceReceiptV1','RJAugmentationReceiptV1'),'unqualified source proof')
    with uproot.open(path) as f:
        source=json.loads(str(f['meta/source']));schema=json.loads(str(f['meta/schema']));tree=f['eventTree']
        require(tree.num_entries==proof['entries'] and source.get('logical_source_id')==proof['logical_source_id'],'certificate source differs')
        identity={k:source.get(k) for k in ('dataset','run','sample','period','si_di_role')}
        data=source.get('si_di_role')=='DATA';index=None
        if data:
            count=0;p=Path(index_path)
            with p.open('xb') as stream:
                for a in tree.iterate(['evt_run','evt_phys_seq','evt_phys_seq_valid'],step_size='8 MB',library='np'):
                    require(bool((a['evt_phys_seq_valid']==1).all()),'DATA source lacks physical event identity')
                    require(bool((a['evt_run']==source['run']).all()),'DATA source event run differs')
                    np.column_stack((a['evt_run'],a['evt_phys_seq'])).astype('<i8').tofile(stream)
                    count+=len(a['evt_run'])
            require(count==tree.num_entries,'physical identity index coverage differs')
            index=dict(path=str(p),sha256=sha256(p),entries=count,encoding='little_endian_int64_run_seq')
        return dict(schema='RJSavedSourceCertificateV1',status='PASS',logical_source_id=proof['logical_source_id'],
            output_sha256=proof['output_sha256'],entries=tree.num_entries,source=identity,
            schema_sha256=schema['schema_sha256'],identity_index=index)


def population_check(certificates,index_path):
    """Bounded SQLite uniqueness over 16-byte DATA keys; SIM needs no ROOT read."""
    db=sqlite3.connect(index_path);db.execute('PRAGMA cache_size=-8192')
    db.execute('CREATE TABLE events(run INTEGER,seq INTEGER,PRIMARY KEY(run,seq)) WITHOUT ROWID')
    reference=None;entries=0;physical=0
    try:
        for c in certificates:
            require(c.get('schema')=='RJSavedSourceCertificateV1' and c.get('status')=='PASS','missing saved-source certificate')
            s=c['source'];data=s['si_di_role']=='DATA'
            identity=tuple(s.get(k) for k in ('dataset','run' if data else 'sample','period','si_di_role'))
            if reference is None:reference=identity
            require(reference==identity,'group mixes DATA runs or SIM samples')
            entries+=c['entries'];index=c['identity_index']
            if data:
                require(index and index['entries']==c['entries'] and index['encoding']=='little_endian_int64_run_seq','DATA identity certificate incomplete')
                require(sha256(index['path'])==index['sha256'],'DATA identity index changed')
                require(Path(index['path']).stat().st_size==16*c['entries'],'DATA identity index size differs')
                with Path(index['path']).open('rb') as stream:
                    while True:
                        block=np.fromfile(stream,dtype='<i8',count=16384)
                        if not len(block):break
                        require(len(block)%2==0,'truncated DATA identity index')
                        rows=block.reshape(-1,2)
                        require(bool((rows[:,0]==s['run']).all()),'identity index names another run')
                        try:db.executemany('INSERT INTO events VALUES(?,?)',((int(r),int(e)) for r,e in rows))
                        except sqlite3.IntegrityError as exc:raise ValueError('duplicate physical DATA event across logical sources') from exc
                        physical+=len(rows)
                db.commit()
                require(sha256(index['path'])==index['sha256'],'identity index changed during publication')
            else:require(index is None,'SIM identity must remain source-scoped')
        return dict(status='PASS',entries=entries,unique_data_events=physical,logical_sources=len(certificates))
    finally:db.close()
