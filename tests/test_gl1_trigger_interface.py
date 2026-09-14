"""Direct GL1 transport and failure semantics, including both package layouts."""
from pathlib import Path
import tempfile
import unittest
from unittest import mock
import numpy as np
import uproot
import build_photonjet_collaboration_tree as builder
import validate_photonjet_collaboration_tree as validator
import gl1_trigger_interface as interface
import test_collaborator_interface as fixtures


def decisions(events, simulation=False):
    n=len(events['event_id_hi'])
    events.update({key:np.zeros(n,dtype=dtype) for key,dtype in interface.EVENT_SCHEMA.items()})
    events['trigger_capture_version'][:]=1
    for key in interface.WORDS:events[key]=np.zeros(n,dtype='uint64')
    events['run'][:]=69577
    if simulation:
        events['trigger_decision_state'][:]=5
    else:
        events['trigger_packet_version'][:]=3
        events['trigger_decision_available'][:]=7
        events['live_trigger_bits'][:]=np.uint64(2**63+2**30+1)
        events['trigger_bits'][:]=events['live_trigger_bits']
        # Genuine zero, MB-only, photon-trigger and bit 63 events all remain.
        events['scaled_trigger_bits'][:]=np.resize(np.asarray([0,1,2**30,2**63],dtype='uint64'),n)
    events['scaled_bit30']=((events['scaled_trigger_bits']>>np.uint64(30))&np.uint64(1)).astype('int32')
    return events


def run_rows():
    rows={key:np.full(64,'a'*64,dtype=object) if dtype=='string' else np.zeros(64,dtype=dtype)
          for key,dtype in interface.RUN_SCHEMA.items() if key!='source_file_index'}
    rows['run'][:]=69577
    rows['bit'][:]=np.arange(64)
    rows['name'][:]=[f'bit {bit}' for bit in range(64)]
    rows['initial_prescale'][:]=np.arange(64)
    rows['prescale'][:]=np.arange(64)+0.5
    for key in ('run_raw','run_live','run_scaled'):
        rows[key][:]=np.arange(64,dtype='uint64')+np.uint64(2**63+7)
    return rows


class DirectGl1InterfaceTests(unittest.TestCase):
    def test_all_events_survive_both_layouts_and_direct_readback(self):
        original_fixture=fixtures.event_fixture
        for system in ('pp','auau'):
            for simulation in (False,True):
                for layout in ('expanded','normalized_v1'):
                    with self.subTest(system=system,simulation=simulation,layout=layout),tempfile.TemporaryDirectory() as tmp:
                        source,output=Path(tmp)/'input.root',Path(tmp)/'out.root'
                        def fixture(complete=True):return decisions(original_fixture(complete),simulation)
                        with mock.patch.object(fixtures,'event_fixture',fixture):
                            fixtures.write_source(source,system=system,simulation=simulation)
                        if not simulation:
                            with uproot.update(source) as root:
                                root.mktree(interface.RUN_TREE,{k:v for k,v in interface.RUN_SCHEMA.items()
                                                               if k!='source_file_index'}).extend(run_rows())
                        counts=builder.build([source],output,system,layout=layout)
                        self.assertEqual(counts['events'],8)
                        self.assertEqual(counts['triggerRunInfo'],0 if simulation else 64)
                        report=validator.validate([source],output,system)
                        self.assertIn('direct_gl1_decisions_and_run_configuration=PASS',report)
                        with uproot.open(output) as root:
                            rows=root['events'].arrays(list(interface.EVENT_SCHEMA)+[*interface.WORDS,'scaled_bit30'],library='np')
                            self.assertEqual(rows['scaled_trigger_bits'].dtype,np.dtype('uint64'))
                            np.testing.assert_array_equal(rows['scaled_bit30'],[-1]*8 if simulation else [0,0,1,0]*2)
                            if not simulation:self.assertEqual(int(rows['scaled_trigger_bits'][3]),2**63)

    def test_missing_is_not_valid_zero_and_simulation_has_no_requirement(self):
        rows=decisions(fixtures.event_fixture())
        self.assertEqual(interface.event_columns(rows)['trigger_decision_state'][0],0)
        rows['trigger_decision_state'][:]=1
        rows['trigger_packet_version'][:]=0
        rows['trigger_decision_available'][:]=0
        for key in interface.WORDS:rows[key][:]=0
        self.assertEqual(interface.event_columns(rows)['trigger_decision_state'][0],1)
        legacy=fixtures.event_fixture()
        self.assertTrue(np.all(interface.event_columns(legacy)['trigger_decision_state']==-1))

    def test_partial_lossy_and_contradictory_states_fail(self):
        for defect in ('partial','float_word','live_alias','missing_with_bits','false_bad','unknown_version'):
            with self.subTest(defect=defect):
                rows=decisions(fixtures.event_fixture())
                if defect=='partial':rows.pop('trigger_packet_version')
                elif defect=='float_word':rows['scaled_trigger_bits']=rows['scaled_trigger_bits'].astype('float64')
                elif defect=='live_alias':rows['trigger_bits'][0]=0
                elif defect=='missing_with_bits':rows['trigger_decision_state'][0]=1
                elif defect=='false_bad':rows['trigger_decision_state'][0]=4
                else:rows['trigger_capture_version'][0]=99
                with self.assertRaises(ValueError):interface.event_columns(rows)

    def test_data_cannot_silently_drop_run_configuration(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'input.root'
            with uproot.recreate(path):pass
            with uproot.open(path) as root:
                with self.assertRaisesRegex(ValueError,'omitted native run configuration'):
                    interface.copy_run_info(root,None,decisions(fixtures.event_fixture()),0)

    def test_invalid_source_or_bit_inventory_fails(self):
        for defect in ('run','bit','partial'):
            with self.subTest(defect=defect),tempfile.TemporaryDirectory() as tmp:
                values=run_rows()
                if defect=='run':values['run'][0]=68531
                elif defect=='bit':values['bit'][1]=0
                else:values={k:v[:63] for k,v in values.items()}
                source,output=Path(tmp)/'source.root',Path(tmp)/'out.root'
                with uproot.recreate(source) as root:
                    root.mktree(interface.RUN_TREE,{k:v for k,v in interface.RUN_SCHEMA.items()
                                                   if k!='source_file_index'}).extend(values)
                with uproot.open(source) as root,uproot.recreate(output) as out:
                    with self.assertRaises(ValueError):
                        interface.copy_run_info(root,out.mktree('triggerRunInfo',interface.RUN_SCHEMA),
                                                decisions(fixtures.event_fixture()),0)


if __name__=='__main__':unittest.main()
