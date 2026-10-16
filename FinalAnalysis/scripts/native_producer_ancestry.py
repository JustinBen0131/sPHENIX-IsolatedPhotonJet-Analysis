"""Verify recorded producer ancestry against retained qualified source releases.

Unknown input ancestry cannot be replaced by a newer executable witness. This
module reads no files and never treats provenance hashes as scientific values.
"""

def ancestry_errors(witness, ancestry, inventories):
    errors=[]; remaining=100000
    def visit(w, a, label, depth):
        nonlocal remaining
        remaining-=1
        if depth>128 or remaining<0:
            errors.append(label+": ancestry exceeds validation bounds"); return
        if (not isinstance(w,dict) or w.get('schema')!='NativeProducerSourceWitnessV1' or
            w.get('state')!='recorded' or not w.get('source_inventory') or
            w.get('source_inventory') not in inventories):
            errors.append(label+": producing code is not a retained qualified release")
        if (not isinstance(a,dict) or a.get('schema')!='NativeProducerAncestryV1' or
            a.get('state')!='recorded' or not isinstance(a.get('inputs'),list)):
            errors.append(label+": input ancestry was not recorded"); return
        for index,row in enumerate(a['inputs']):
            if not isinstance(row,dict) or not isinstance(row.get('source_path'),str) or not row['source_path']:
                errors.append(label+": malformed parent record"); continue
            visit(row.get('witness'),row.get('ancestry'),label+'/'+str(index),depth+1)
    visit(witness,ancestry,'producer',0)
    return errors
