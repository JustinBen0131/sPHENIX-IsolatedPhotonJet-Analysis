# TreeProduction: DST to collaborator trees

Two steps: a Fun4All producer that reads DSTs and writes a normalised event
record (`ReplayFoundationV1/*` trees plus the direct histograms of the
maintained analysis), and a Python builder that turns that record into the
eight collaborator trees described in `contracts/photonjet_trees_v1_branches.json`.

## Producer

* `src/` p+p analysis module (`RecoilJets`), `src_AuAu/` Au+Au module
  (`RecoilJets_AuAu`).  Both build with the usual sPHENIX autotools flow:
  `./autogen.sh --prefix=$MYINSTALL && make install`.
* `macros/Fun4All_recoilJets.C` and `macros/Fun4All_recoilJets_AuAu.C` are the
  entry points; both include `macros/Fun4All_recoilJets_unified_impl.C`.
  Arguments: number of events (0 = all), DST list file, output ROOT file,
  verbose flag, events to skip.
* `macros/analysis_config.yaml` is the producer configuration (photon and jet
  windows, jet radii 0.2/0.3/0.4, vertex windows, isolation cones, trigger
  selection).  The per-job YAML text is stamped into every output file.
* `macros/Calo_Calib.C` builds calibrated calorimeter towers and clusters;
  `macros/calo/` holds the tower-status and centrality helpers it uses.

### Site paths to edit before building elsewhere

The producer was last run from one user installation.  These lines contain
absolute paths that must point at your own install/scratch area:

* `macros/Fun4All_recoilJets.C`, `macros/Fun4All_recoilJets_AuAu.C`: the
  `#include` of the unified implementation macro.
* `macros/Fun4All_recoilJets_unified_impl.C`: `R__ADD_INCLUDE_PATH` /
  `pragma cling add_include_path` lines and the `PhotonClusterBuilder.h`
  include near the top.
* `macros/Calo_Calib.C`: one include path.
* `macros/analysis_config.yaml`: two reweighting-file paths.
* `src_AuAu/RecoilJets_AuAu.cc`: one include path.

`grep -n "/sphenix/u/\|/gpfs/\|/sphenix/tg/" -r TreeProduction` lists them.

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
* The validator checks branch types, unique event/candidate identities, the
  event-array/flat-tree equivalence, pair kinematics and truth-link targets.

Helper modules used by the builder: `centrality_replay.py`,
`gl1_trigger_interface.py`, `trigger_scaler_interface.py`, `producer/`.
