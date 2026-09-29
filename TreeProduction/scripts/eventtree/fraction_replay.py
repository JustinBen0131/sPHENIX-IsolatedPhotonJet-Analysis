"""Bounded, exact join of native replay observations to a stage-2 eventTree.

SQLite holds the native jet pool on scratch. Retained jets must match one and
only one native raw-jet fingerprint, with the same physical event and radius.
Extra native jets are allowed because production's stored collection has cuts;
their number is reported. Corrected pT is never reconstructed or overwritten.
"""
import json
from pathlib import Path
import sqlite3

import numpy as np
import uproot

from .augmentation import FINGERPRINT, RADII
from .integrity import sha256
from .tables import require


def fingerprint(row):
    a=np.asarray([row[k] for k in FINGERPRINT],dtype=np.float64)
    require(np.isfinite(a).all(),'nonfinite native/retained jet fingerprint')
    return a.tobytes()


def join(input_path, native_paths, output, *, index_path, progress=None):
    require(not Path(index_path).exists(),'replay scratch index already exists')
    before=sha256(input_path); native_pins=[(str(p),sha256(p)) for p in native_paths]
    db=sqlite3.connect(index_path)
    db.executescript('''PRAGMA journal_mode=DELETE; PRAGMA cache_size=-8192;
      CREATE TABLE events(run INTEGER,seq INTEGER,radius TEXT,PRIMARY KEY(run,seq,radius));
      CREATE TABLE jets(run INTEGER,seq INTEGER,radius TEXT,jet_index INTEGER,fp BLOB,body TEXT,used INTEGER DEFAULT 0,
                        PRIMARY KEY(run,seq,radius,jet_index));''')
    try:
        native_count=0
        for path,pin in native_pins:
            with open(path) as f:
                for line in f:
                    r=json.loads(line);key=(r['run'],r['physical_event'],r['radius'])
                    require(r['radius'] in RADII,'native replay radius outside recipe')
                    try:
                        if r['kind']=='event':db.execute('INSERT INTO events VALUES(?,?,?)',key)
                        elif r['kind']=='jet':
                            require(db.execute('SELECT 1 FROM events WHERE run=? AND seq=? AND radius=?',key).fetchone(),
                                    'native jet precedes event marker')
                            db.execute('INSERT INTO jets(run,seq,radius,jet_index,fp,body) VALUES(?,?,?,?,?,?)',
                                       (*key,int(r['jet_index']),fingerprint(r),json.dumps(r)))
                            native_count+=1
                            if native_count%10000==0:db.commit()
                        else:raise ValueError('unknown native replay record')
                    except sqlite3.IntegrityError as exc:
                        raise ValueError('duplicate/ambiguous replay event or jet') from exc
            require(sha256(path)==pin,'native replay changed while indexing')
        db.commit();matched=0;events=0
        with uproot.open(input_path) as source,Path(output).open('x') as out:
            names=['evt_run','evt_phys_seq','evt_phys_seq_valid','evt_src_entry']
            names+=['jet'+r+'_'+k for r in RADII for k in FINGERPRINT]
            for arrays in source['eventTree'].iterate(names,step_size='8 MB',library='ak'):
                for e in range(len(arrays)):
                    for radius in RADII:
                        for j in range(len(arrays['jet'+radius+'_rawpt'][e])):
                            require(int(arrays['evt_phys_seq_valid'][e])==1,'retained jet has no physical identity')
                            key=(int(arrays['evt_run'][e]),int(arrays['evt_phys_seq'][e]),radius)
                            values={k:float(arrays['jet'+radius+'_'+k][e][j]) for k in FINGERPRINT}
                            fp=fingerprint(values)
                            found=db.execute('SELECT body,used FROM jets WHERE run=? AND seq=? AND radius=? AND jet_index=? AND fp=?',(*key,j,fp)).fetchone()
                            require(found is not None and found[1]==0,'missing, mismatched or multiply used replay jet')
                            row=json.loads(found[0]);row.update(source_entry=int(arrays['evt_src_entry'][e]),jet_index=j)
                            row['jet_energy']=float('nan') if row['jet_energy'] is None else row['jet_energy']
                            row['calo_energy']=[float('nan') if v is None else v for v in row['calo_energy']]
                            out.write(json.dumps(row,sort_keys=True)+'\n')
                            db.execute('UPDATE jets SET used=1 WHERE run=? AND seq=? AND radius=? AND jet_index=? AND fp=?',(*key,j,fp));matched+=1
                events+=len(arrays);db.commit()
                if progress:progress(events,source['eventTree'].num_entries)
            out.flush()
            import os
            os.fsync(out.fileno())
        require(sha256(input_path)==before,'stage-2 file changed during replay join')
        return dict(schema='RJFractionReplayJoinV1',status='PASS',input_sha256=before,events=events,
                    matched_jets=matched,unretained_native_jets=native_count-matched,native_pins=native_pins,
                    output=str(output),output_sha256=sha256(output))
    finally:db.close()
