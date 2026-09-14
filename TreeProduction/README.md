# TreeProduction: DST to collaborator trees

Two steps: a Fun4All producer that reads DSTs and writes a normalised event
record (the `ReplayFoundationV1` directory of `RJ*V1` record trees, plus the
direct histograms of the maintained analysis), and a Python builder that turns
that record into the eight collaborator trees described in
`contracts/photonjet_trees_v1_branches.json`.

## Where this source comes from

The producer here is a copy of the production source frozen on 2026-09-13,
before the current Au+Au production, passed through
`tools/scrub_producer_source.py`.  It has not been rebuilt from this
repository: there is no sPHENIX build environment on the machine that
assembled it.  When a production freezes, the rule is: copy the exact producer
source and tree builder that ran it into this directory verbatim, record their
hashes, rerun the scrub script, and commit.  `python tests/run_tests.py
test_producer_scrub.py` checks the result.

What the scrub changes is documented at the top of the script: study numbers
in comments, one user's absolute paths.  What it keeps, on purpose: C++
identifiers, environment-variable names, ROOT object names and log tags that
still carry study numbers (`the44`, `the134`, `the119`, `the221`).  Those are
runtime interfaces of the production wrappers and of the output files, and
they have no physics meaning.  `THE106Observation` is the name of an optional
instrumentation header in a patched PHOOL; the local
`THE106ObservationDisabled.h` compiles those call sites to no-ops when that
header is absent.

## Producer

* `src/` p+p analysis module (`RecoilJets`), `src_AuAu/` Au+Au module
  (`RecoilJets_AuAu`).  Both build with the usual sPHENIX autotools flow after
  `source /opt/sphenix/core/bin/sphenix_setup.sh -n` and
  `source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL`:

  ```bash
  cd TreeProduction/src && ./autogen.sh --prefix=$MYINSTALL && make -j4 install && cd ../..
  cd TreeProduction/src_AuAu && ./autogen.sh --prefix=$MYINSTALL && make -j4 install && cd ../..
  ```

  `src/` also carries the analysis copy of `PhotonClusterBuilder` (installed
  as part of `libcalo_reco`), which the macros load by name from the
  installation.
* `macros/Fun4All_recoilJets.C` and `macros/Fun4All_recoilJets_AuAu.C` are the
  entry points; both include `macros/Fun4All_recoilJets_unified_impl.C`.
  Arguments: number of events (0 = all), DST list file, output ROOT file,
  verbose flag, events to skip.
* `macros/analysis_config.yaml` is the producer configuration (photon and jet
  windows, jet radii 0.2/0.3/0.4, vertex windows, isolation cones, trigger
  selection, reweighting inputs).  The per-job YAML text is stamped into every
  output file.
* `macros/Calo_Calib.C` builds calibrated calorimeter towers and clusters;
  `macros/calo/` holds the tower-status and centrality helpers it uses.

### External inputs

These files are not in the repository.  They are analysis inputs provided by
collaborators and are referenced by their location on the sPHENIX file system:

| Input | Configured in | Provider |
| --- | --- | --- |
| p+p truth-vertex reweighting histograms (0 and 1.5 mrad periods) | `analysis_config.yaml` `vertex_reweight_file_pp`; period contracts in `src/RecoilJets.cc` | Shuhang Li (PPG12 efficiency tools) |
| Au+Au vertex and centrality reweighting histograms (off by default) | `analysis_config.yaml` `vertex_reweight_file_auau`, `centrality_reweight_file` | Blair Seidlitz |
| EMCal tower mask for the PPG12 yield replay | `src/RecoilJets.cc` `kPPG12YieldTowerMaskFile` | Shuhang Li |
| Scaled-trigger QA run list (optional) | environment `RJ_SCALED_TRIGGER_RUNLIST` | this analysis |

## Collaborator tree builder

```bash
python TreeProduction/build_photonjet_collaboration_tree.py --system auau \
    --input auau_producer_part1.root --input auau_producer_part2.root \
    --output trees/auau_data.root
python TreeProduction/validate_photonjet_collaboration_tree.py --system auau trees/auau_data.root
```

* Inputs are producer files; `--input-list` accepts a text file of paths.
* No model is needed: `bdt_score` is written unscored (`bdt_evaluation_state`
  0) and filled later by `PhotonID/score_trees.py`, which also regenerates the
  threshold and flag branches from `config/nominal.yaml`.  (`--model` is
  accepted for backward compatibility with the maintained analysis.)
* `--require-*` options enforce the presence of the trigger scaler, MBD and
  centrality-replay branches used by the Au+Au normalisation; use them for
  data production.
* The validator checks branch types, unique event and candidate identities,
  the event-array to flat-tree equivalence, pair kinematics and truth-link
  targets.

Helper modules used by the builder: `centrality_replay.py`,
`gl1_trigger_interface.py`, `trigger_scaler_interface.py`, `producer/`.
