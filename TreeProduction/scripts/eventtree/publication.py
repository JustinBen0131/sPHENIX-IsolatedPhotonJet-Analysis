"""Immutable complete-group versions and a compare-and-swap publication pointer.

Readers resolve CURRENT exactly once before collecting their input files. A
failed augmentation leaves the previous version accessible and unchanged.
"""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import uuid
import sqlite3
import uproot

from .integrity import sha256
from .tables import require


def population_check(files,index_path):
    """Run-wide DATA physical identity uniqueness; SIM identity stays source-scoped."""
    db=sqlite3.connect(index_path);db.execute('PRAGMA cache_size=-8192')
    db.execute('CREATE TABLE events(run INTEGER,seq INTEGER,PRIMARY KEY(run,seq))')
    reference=None;entries=0;physical=0
    try:
        for item in files:
            with uproot.open(item['path']) as f:
                source=json.loads(str(f['meta/source']));role=source.get('si_di_role')
                identity=tuple(source.get(k) for k in ('dataset','run' if role=='DATA' else 'sample','period','si_di_role'))
                if reference is None:reference=identity
                require(identity==reference,'group mixes DATA runs or SIM sample identity')
                require(source.get('logical_source_id')==item['logical_id'],'file logical source differs from group manifest')
                for a in f['eventTree'].iterate(['evt_run','evt_phys_seq','evt_phys_seq_valid'],step_size='8 MB',library='np'):
                    entries+=len(a['evt_run'])
                    if role=='DATA':
                        require(bool((a['evt_phys_seq_valid']==1).all()),'DATA group has unverifiable physical event identity')
                        try:db.executemany('INSERT INTO events VALUES(?,?)',((int(r),int(s)) for r,s in zip(a['evt_run'],a['evt_phys_seq'])))
                        except sqlite3.IntegrityError as exc:raise ValueError('duplicate physical DATA event across logical sources') from exc
                        physical+=len(a['evt_run']);db.commit()
        return dict(status='PASS',entries=entries,unique_data_events=physical,logical_sources=len(files))
    finally:db.close()


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',', ':'),allow_nan=False)


def safe_relative(name):
    p=Path(name)
    require(not p.is_absolute() and '..' not in p.parts and p.suffix=='.root','unsafe logical output name')
    return p


def resolve(group):
    """Pin a version once; callers must use this returned directory for the run."""
    group=Path(group).resolve(); current=group/'CURRENT'
    require(current.is_symlink(),'group is not published')
    version=current.resolve(strict=True)
    require(version.parent==group/'versions','publication points outside immutable version root')
    manifest=json.loads((version/'MANIFEST.json').read_text())
    body={k:v for k,v in manifest.items() if k!='version'}
    require(hashlib.sha256(canonical(body).encode()).hexdigest()==manifest['version']==version.name,'group manifest changed')
    return version,manifest


def publish(group, files, *, group_id, expected_sources, recipe, base_version=None, stage=2, validator=None, exclusions=()):
    """Files contain logical_id, relative_path, path, sha256 and a PASS receipt pin.

    Rereads every saved file. A supplied maintained validator runs before any
    publication; missing validation refuses the operation. Pointer update holds
    a local/farm advisory lock and compares the exact expected base version.
    """
    require(stage in (2,3),'invalid publication stage')
    expected=set(expected_sources)
    require(len(files)==len(expected) and {f['logical_id'] for f in files}==expected,'group logical-source roster differs')
    require(len({f['relative_path'] for f in files})==len(files),'colliding logical file names')
    group=Path(group).resolve(); group.mkdir(parents=True,exist_ok=True)
    records=[];certificates=[];identities={}
    def file_identity(p):
        s=Path(p).stat();return (s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns)
    for f in sorted(files,key=lambda x:x['logical_id']):
        safe_relative(f['relative_path'])
        before=file_identity(f['path'])
        require(before[0]==group.stat().st_dev,'publication candidate and version root must share a filesystem')
        require(sha256(f['path'])==f['sha256'],'group input hash differs')
        require(file_identity(f['path'])==before,'group input changed during hashing')
        identities[f['logical_id']]=before
        require(sha256(f['receipt_path'])==f['receipt_sha256'],'group validation receipt changed')
        r=json.loads(Path(f['receipt_path']).read_text())
        require(r['status']=='PASS' and r['logical_source_id']==f['logical_id'] and r['output_sha256']==f['sha256'],
                'source validation receipt does not bind file')
        if validator is not None:require(validator(f['path'])['ok'],'saved-file validation failed')
        else:
            c=r.get('saved_source_certificate')
            require(c and c.get('logical_source_id')==f['logical_id'] and c.get('output_sha256')==f['sha256']
                    and c.get('entries')==r['entries'] and r.get('proof',{}).get('status')=='PASS'
                    and r.get('nominal',{}).get('status')=='PASS','bounded publication requires complete source validation certificate')
            certificates.append(c)
        records.append({k:f[k] for k in ('logical_id','relative_path','sha256','receipt_sha256')})
    checks=group/'.validation';checks.mkdir(exist_ok=True)
    if validator is None:
        from .certificate import population_check as certified_population
        population=certified_population(certificates,checks/(uuid.uuid4().hex+'.sqlite'))
    else:population=population_check(files,checks/(uuid.uuid4().hex+'.sqlite'))
    body=dict(schema='RJDatasetGroupV1',group_id=group_id,stage=stage,recipe=recipe,base_version=base_version,files=records,population=population,exclusions=list(exclusions))
    version=hashlib.sha256(canonical(body).encode()).hexdigest()
    manifest=dict(body,version=version)
    folder=group/'versions'/version
    folder.mkdir(parents=True,exist_ok=True)
    for f in files:
        target=folder/safe_relative(f['relative_path']); target.parent.mkdir(parents=True,exist_ok=True)
        try: os.link(f['path'],target)
        except FileExistsError:require(file_identity(target)==identities[f['logical_id']],'version file collision')
        require(file_identity(target)==identities[f['logical_id']],'version file changed while publishing')
    raw=canonical(manifest)+'\n'; m=folder/'MANIFEST.json'
    candidate=folder/('.MANIFEST.'+uuid.uuid4().hex)
    with candidate.open('x') as stream:
        stream.write(raw);stream.flush();os.fsync(stream.fileno())
    try:
        os.link(candidate,m)
    except FileExistsError: require(m.read_text()==raw,'manifest collision')
    # Lock and pointer live under this one group, so independent samples proceed.
    with (group/'.publication.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        old,previous=(resolve(group) if (group/'CURRENT').is_symlink() else (None,None))
        actual=previous['version'] if previous else None
        if actual==version: return manifest
        require(actual==base_version,'publication race: base version changed')
        if stage==3:
            require(previous is not None,'augmentation requires a published base')
            require([(r['logical_id'],r['relative_path']) for r in previous['files']]==
                    [(r['logical_id'],r['relative_path']) for r in records],'augmentation changes logical file set')
        # Recheck the linked saved version while holding the publication lock.
        # Candidates are create-once and may never be edited after validation.
        for f in files:
            target=folder/safe_relative(f['relative_path'])
            require(file_identity(target)==identities[f['logical_id']],'candidate changed before pointer publication')
            fd=os.open(target,os.O_RDONLY)
            try:os.fsync(fd)
            finally:os.close(fd)
        tmp=group/('.CURRENT.'+uuid.uuid4().hex)
        tmp.symlink_to(Path('versions')/version,target_is_directory=True)
        os.replace(tmp,group/'CURRENT')
        fd=os.open(group,os.O_RDONLY)
        try: os.fsync(fd)
        finally: os.close(fd)
    return manifest
